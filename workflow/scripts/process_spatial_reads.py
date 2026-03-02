import argparse
import pandas as pd
import pysam
from Bio.Seq import Seq
from Bio import SeqIO
import re
import swifter  # pip install swifter


# ARGPARSE
# -----------------------------------------------------
parser = argparse.ArgumentParser(
    description="Process PacBio BAM files and produce filtered FASTQ"
)

parser.add_argument(
    "-i", "--input", required=True, help="Path to input BAM file"
)

parser.add_argument(
    "-o", "--output", required=True, help="Path to output FASTQ file"
)
args = parser.parse_args()

path_bamfile = args.input
output_fastq = args.output

def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def find_adapter(readstring):
    A_seq = 'CGACGCTCTTCCGATCT'
    A_len = len(A_seq)
    distances = [(offset+A_len,hamming_distance(readstring[offset:offset+A_len], A_seq)) for offset in range(7)]
    argmin   = min(distances, key=lambda x: x[1])
    return argmin

def find_TSO(readstring):
    T_seq = 'TCTGCGTTGATACCACT'
    T_len = len(T_seq)
    distances = [(-(offset+T_len+10),hamming_distance(readstring[-(offset+T_len):-offset], T_seq)) for offset in range(7)]
    argmin   = min(distances, key=lambda x: x[1])
    return argmin

# Load BAM and create DataFrame as before
bamfile = pysam.AlignmentFile(path_bamfile, 'rb', check_sq=False)

records = []
for read in bamfile.fetch(until_eof=True):
    
    
    records.append({
        "qname": read.query_name,
        "seq": read.query_sequence,
        "qual": ''.join(chr(q + 33) for q in read.query_qualities)#read.query_qualities,
    })

df = pd.DataFrame(records)

# Reverse complement with swifter parallel apply
df['rcseq'] = df['seq'].swifter.apply(lambda x: str(Seq(x).reverse_complement()))

# Find adapters forward and reverse in parallel
df['forward_distance'] = df['seq'].swifter.apply(find_adapter)
df['reverse_distance'] = df['rcseq'].swifter.apply(find_adapter)

# Direction, min_distance, barcode_pos
df['direction'] = df.swifter.apply(lambda x: 'R' if x['reverse_distance'][1] < x['forward_distance'][1] else 'F', axis=1)
df['min_distance'] = df.apply(lambda x: x['reverse_distance'][1] if x['direction'] == 'R' else x['forward_distance'][1], axis=1)
df['barcode_pos'] = df.apply(lambda x: x['reverse_distance'][0] if x['direction'] == 'R' else x['forward_distance'][0], axis=1)

# cseq column depending on direction
df['cseq'] = df.swifter.apply(lambda x: x['seq'] if x['direction'] == 'F' else x['rcseq'], axis=1)

# Drop unnecessary columns
df = df.drop(columns=['seq', 'rcseq', 'reverse_distance', 'forward_distance'])

# find TSO in parallel
df['tso'] = df['cseq'].swifter.apply(find_TSO)

# Extract barcode, UMI, polyT regions in parallel
df['st_barcode'] = df.swifter.apply(lambda x: x['cseq'][x['barcode_pos']:x['barcode_pos']+16], axis=1)
df['st_umi'] = df.swifter.apply(lambda x: x['cseq'][x['barcode_pos']+16:x['barcode_pos']+28], axis=1)
df['st_polyT'] = df.swifter.apply(lambda x: x['cseq'][x['barcode_pos']+28:x['barcode_pos']+32], axis=1)

# Filter rows where UMI is not all T and polyT is all T
df = df[(df.st_umi != 'TTTTTTTTTTTT') & (df.st_polyT == 'TTTT')]

# Precompile regex for filtering
pattern = re.compile('[^T]T{0,2}[^T]T{0,2}[^T]')

# Use swifter parallel apply for regex filtering
def pattern_match(row):
    substring = row.cseq[row.barcode_pos + 28:]
    return bool(pattern.search(substring))

df = df[df.swifter.apply(pattern_match, axis=1)].reset_index(drop=True)

def compute_polyT_end_pos(row):
    substring = row.cseq[row.barcode_pos + 28:]
    match = pattern.search(substring)
    return match.span()[0] + row.barcode_pos + 28 if match else -1  # Use -1 or np.nan for missing

df['polyT_end_pos'] = df.swifter.apply(compute_polyT_end_pos, axis=1)

df['filtered_seq'] = df.swifter.apply(
    lambda row: row.cseq[row.polyT_end_pos:-31] if row.polyT_end_pos >= 0 else '', axis=1
)

df['filtered_qual'] = df.swifter.apply(
    lambda row: row.qual[row.polyT_end_pos:-31] if row.polyT_end_pos >= 0 else '', axis=1
)


#df = df.assign(polyT_end_pos = [re.search('[^T]T{0,2}[^T]T{0,2}[^T]', x.cseq[x.barcode_pos+28:]).span()[0]+x.barcode_pos+28 for i,x in df.iterrows()])

#df = df.assign(filtered_seq = [x.cseq[x.polyT_end_pos:-31] for i, x in df.iterrows()])

# Load spottab once
spottab = pd.read_csv(
    '/storage/mathelierarea/processed/eirikhoy/vdj_spatial/resources/spaceranger-4.0.1/lib/python/cellranger/barcodes/visium-v1_coordinates.txt',
    sep='\t', header=None, index_col=0, names=['barcode', 'x', 'y']
)
bset = set(spottab.index)

# Add coordinates with vectorized logic
def get_y(x):
    if x in bset:
        return spottab.loc[x, 'y'] - 1
    else:
        return None

def get_x(x):
    if x in bset:
        return spottab.loc[x, 'x'] - 1
    else:
        return None

df['y_coor'] = df['st_barcode'].swifter.apply(get_y)
df['x_coor'] = df['st_barcode'].swifter.apply(get_x)
df['optional'] = '+'
df['read_id'] = df.swifter.apply(
    lambda row: f"{row.qname}:x_coor={row.x_coor}_y_coor={row.y_coor}", axis=1
)

def write_fastq(df, output_path):
    with open(output_path, 'w') as f:
        for read_id, seq, optional, qual in zip(
            df['read_id'], df['filtered_seq'], df['optional'], df['filtered_qual']
        ):
            f.write(f"@{read_id}\n{seq}\n{optional}\n{qual}\n")

write_fastq(df, output_fastq)


# TODO Concider using samtools or bedtools 
# TODO try bedder?
# TODO Try not aligning to J region, if there is an option
# TODO Look for alternative pipeline, nextflow or random repos published somewhere

