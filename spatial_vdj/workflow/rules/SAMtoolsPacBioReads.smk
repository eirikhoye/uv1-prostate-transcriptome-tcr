rule process_pacbio_reads:
    input:
        bamfiles = BAM_FILES
    
    output:
        fastqfiles = FASTQ_FILES
    
    shell:
        """
        {input.bamfiles} {output.fastqfiles}
        
        """        


# TODO Implement qc stats rule using the output from samtools

# module load SAMtools/1.21-GCC-13.3.0

#  1007  2026-02-10T11:12:04  samtools view /storage/mathelierarea/raw/urbanucci_vdj_spatial/data/BIDGEN_02_12_25/segmented.bam | head
# 1008  2026-02-10T11:14:12  samtools view /storage/mathelierarea/raw/urbanucci_vdj_spatial/data/BIDGEN_02_12_25/segmented.bam | head -n 1 | tr '\t' '\n'
# 1009  2026-02-10T11:19:05  samtools view -H /storage/mathelierarea/raw/urbanucci_vdj_spatial/data/BIDGEN_02_12_25/segmented.bam
# 1010  2026-02-10T11:20:59  samtools view /storage/mathelierarea/raw/urbanucci_vdj_spatial/data/BIDGEN_02_12_25/segmented.bam | awk '{print length($10)}' | head
# 1011  2026-02-10T11:21:36  samtools stats /storage/mathelierarea/raw/urbanucci_vdj_spatial/data/BIDGEN_02_12_25/segmented.bam | grep ^RL

"""
Perfect — this confirms exactly what we suspected. A few key takeaways from your `samtools view` output:

---

### 1️⃣ These are **unmapped, single-molecule CCS reads**

* `FLAG = 4` → unmapped
* `RNAME = *`, `POS = 0`, `CIGAR = *` → no alignment yet
* Single-end (`no paired reads`)

This is exactly what you get when a PacBio **CCS BAM** is stored post-segmentation but pre-alignment. Everything looks normal.

---

### 2️⃣ They are **segmented reads**

* `qs` / `qe` tags indicate start/end positions within the original polymerase read
* `ds:B:C,...` and `MM:Z`, `ML:B:C` contain **base modifications / consensus info**
* `zm` = ZMW identifier, `rq` = read quality → typical PacBio CCS metadata

This confirms that `segmented.bam` is **already processed into logical read units**.

---

### 3️⃣ The sequences themselves look sane

* Long, variable-length sequences (hundreds to thousands of bp) → consistent with PacBio HiFi CCS reads
* Poly-A tails / adapters present → typical for VDJ or transcriptomic libraries

---

### ✅ Conclusion

Nothing unusual is happening here. This is a **fully expected unaligned, unpaired PacBio BAM**. `samtools flagstat` isn’t informative for this kind of file; it just tells you “everything is unmapped and single-ended,” which is correct.

If your goal is to **analyze read quality, lengths, or VDJ segments**, you should instead focus on:

```bash
samtools stats segmented.bam    # read length distribution, coverage
samtools view segmented.bam | head -n 10 | less -S  # check tags / sequences
```

or feed directly into a **PacBio / VDJ-aware pipeline** for downstream analysis.

---

If you want, I can **show a quick one-liner** to get a histogram of read lengths directly from this BAM — it’s more meaningful than `flagstat` for these PacBio reads. Do you want me to do that?

"""