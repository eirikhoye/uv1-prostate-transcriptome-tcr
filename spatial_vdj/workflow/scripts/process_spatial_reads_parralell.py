import pandas as pd
import pysam
import numpy as np
import re
from numba import njit

# --- Optimized hamming distance using numba ---
@njit
def hamming_distance_numba(s1, s2):
    dist = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            dist += 1
    return dist

# --- Vectorized adapter finder ---
def find_adapter_vec(readstrings, adapter_seq='CGACGCTCTTCCGATCT', max_offset=7):
    adapter_len = len(adapter_seq)
    adapter_bytes = adapter_seq.encode()

    results = []
    for read in readstrings:
        best_pos = None
        best_dist = 1e9
        read_bytes = read.encode()
        for offset in range(max_offset):
            if offset + adapter_len > len(read_bytes):
                continue
            segment = read_bytes[offset:offset+adapter_len]
            dist = hamming_distance_numba(segment, adapter_bytes)
            if dist < best_dist:
                best_dist = dist
                best_pos = offset + adapter_len
        results.append((best_pos, best_dist))
    return results

def find_TSO_vec(readstrings, TSO_seq='TCTGCGTTGATACCACT', max_offset=7):
    T_len = len(TSO_seq)
    T_bytes = TSO_seq.encode()
    results = []
    for read in readstrings:
        best_pos = None
        best_dist = 1e9
        read_bytes = read.encode()
        read_len = len(read_bytes)
        for offset in range(max_offset):
            start = read_len - (offset + T_len)
            end = read_len - offset
            if start < 0 or end > read_len:
                continue
            segment = read_bytes[start:end]
            dist = hamming_distance_numba(segment, T_bytes)
            if dist < best_dist:
                best_dist = dist
                best_pos = -(offset + T_len + 10)  # as original script
        results.append((best_pos, best_dist))
    return results

# --- Reverse complement using numpy ---
_complement_table = np.array(
    [ord('N')] * 256, dtype=np.uint8
)
for a, b in zip(b'ACGTacgt', b'TGCAtgca'):
    _complement_table[a] = b
def reverse_complement(seqs):
    # seqs: list/array of strings
    rc_seqs = []
    for s in seqs:
        arr = np.frombuffer(s.encode(), dtype=np.uint8)
        comp = _complement_table[arr]
        rc = comp[::-1]
        rc_seqs.append(rc.tobytes().decode())
    return rc_seqs

# --- Load BAM file ---
path_bamfile = '/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/hifi_reads/m84212_250410_123503_s3.hifi_reads.bcAd1031T.bam'
bamfile = pysam.AlignmentFile(path_bamfile, 'rb', check_sq=False)

records = []
for read in bamfile.fetch(until_eof=True):
    records.append({
        "qname": read.query_name,
        "seq": read.query_sequence,
        "qual": ''.join(chr(q+33) for q in read.query_qualities),  # convert qual scores to ASCII
    })

df = pd.DataFrame(records)

# --- Reverse complement sequences ---
df['rcseq'] = reverse_complement(df['seq'].values)

# --- Find adapters (vectorized) ---
df['forward_distance'] = find_adapter_vec(df['seq'].values)
df['reverse_distance'] = find_adapter_vec(df['rcseq'].values)

# --- Choose direction vectorized ---
forward_dist = np.array([x[1] for x in df['forward_distance']])
reverse_dist = np.array([x[1] for x in df['reverse_distance']])

direction = np.where(reverse_dist < forward_dist, 'R', 'F')
df['direction'] = direction

min_distance = np.where(direction == 'R', reverse_dist, forward_dist)
df['min_distance'] = min_distance

barcode_pos = np.where(
    direction == 'R',
    [x[0] for x in df['reverse_distance']],
    [x[0] for x in df['forward_distance']]
)
df['barcode_pos'] = barcode_pos

# --- cseq depends on direction ---
df['cseq'] = np.where(df['direction'] == 'F', df['seq'], df['rcseq'])

# Drop unneeded cols early
df.drop(columns=['seq', 'rcseq', 'forward_distance', 'reverse_distance'], inplace=True)

# --- Find TSO vectorized ---
df['tso'] = find_TSO_vec(df['cseq'].values)

# --- Extract barcode, UMI, polyT ---
pos = df['barcode_pos'].values.astype(int)
cseqs = df['cseq'].values

st_barcode = [cseqs[i][p:p+16] for i, p in enumerate(pos)]
st_umi = [cseqs[i][p+16:p+28] for i, p in enumerate(pos)]
st_polyT = [cseqs[i][p+28:p+32] for i, p in enumerate(pos)]

df['st_barcode'] = st_barcode
df['st_umi'] = st_umi
df['st_polyT'] = st_polyT

# --- Filter UMI != all T and polyT == all T ---
df = df[(df['st_umi'] != 'TTTTTTTTTTTT') & (df['st_polyT'] == 'TTTT')]

# --- Precompile regex ---
pattern = re.compile('[^T]T{0,2}[^T]T{0,2}[^T]')

# --- Vectorized regex filter ---
def regex_filter(cseqs, barcode_pos, pattern):
    mask = []
    for cseq, pos in zip(cseqs, barcode_pos):
        substring = cseq[pos + 28:]
        mask.append(bool(pattern.search(substring)))
    return np.array(mask)

regex_mask = regex_filter(df['cseq'].values, df['barcode_pos'].values.astype(int), pattern)
df = df[regex_mask].reset_index(drop=True)

# --- Load spottab once and convert to dict for fast lookup ---
spottab = pd.read_csv(
    '/storage/mathelierarea/processed/eirikhoy/vdj_spatial/resources/spaceranger-4.0.1/lib/python/cellranger/barcodes/visium-v1_coordinates.txt',
    sep='\t', header=None, index_col=0, names=['barcode', 'x', 'y']
)
coords_dict = spottab.to_dict(orient='index')

def get_coord(barcode, axis):
    val = coords_dict.get(barcode)
    if val is None:
        return None
    return val[axis] - 1

df['y_coor'] = df['st_barcode'].map(lambda x: get_coord(x, 'y'))
df['x_coor'] = df['st_barcode'].map(lambda x: get_coord(x, 'x'))
df['optional'] = '+'

# --- Write FASTQ efficiently ---
def write_fastq(df, output_path):
    with open(output_path, 'w') as f:
        for _, row in df.iterrows():
            f.write(f"@{row.qname}\n{row.cseq}\n{row.optional}\n{row.qual}\n")

write_fastq(df, "/storage/mathelierarea/processed/eirikhoy/vdj_spatial/results/m84212_250410_123503_s3.hifi_reads.bcAd1031T_trimmed.fastq")
