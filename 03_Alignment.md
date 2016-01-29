# Edinburgh Genomics - MinION Bioinformatics
## Read alignment

## Aligning reads to a reference

We can either use BWA or LAST.  In our experience BWA tends to be fast but less sensitive, and LAST tends to be slow but more sensitive.

BWA straight to sorted BAM:

```sh
# bwa index genome
bwa index Data/reference/Ecoli_MG1655/MG1655.fa

# create fasta index for samtools
samtools faidx Data/reference/Ecoli_MG1655/MG1655.fa

# run bwa and pipe straight to samtools to create BAM
bwa mem -x ont2d Data/reference/Ecoli_MG1655/MG1655.fa MAP006-1.2D.fastq | \
        samtools view -T Data/reference/Ecoli_MG1655/MG1655.fa -bS - | \
        samtools sort -T 2D_vs_MG1655.bwa -o 2D_vs_MG1655.bwa.bam -

# index bam
samtools index 2D_vs_MG1655.bwa.bam
```

LAST straight to sorted BAM:

```sh
# create a LAST db of the reference genome
lastdb MG1655 Data/reference/Ecoli_MG1655/MG1655.fa

# align high quaklity reads to reference genome with LAST
lastal -q 1 -a 1 -b 1 MG1655 MAP006-1.2D.fasta > MAP006-1.2D.maf

# convert the MAF to BAM with complete CIGAR (matches and mismatches)
maf-convert.py sam MAP006-1.2D.maf | \
    samtools view -T Data/reference/Ecoli_MG1655/MG1655.fa -bS - | \
    samtools sort -T 2D_vs_MG1655.last -o 2D_vs_MG1655.last.bam -

# index bam
samtools index 2D_vs_MG1655.last.bam
```

# Assessing the accuracy of alignment
Aaron Quinlan has written a few useful scripts to work with nanopore data available on [github](https://github.com/arq5x/nanopore-scripts/).  We can use these to generate some statistics on the quality of alignments.  Unfortunately, at present the indel counting only works with LAST not BWA, due to differences in the cigar strings:

```bash
# count errors
count-errors.py 2D_vs_MG1655.last.bam > 2D_vs_MG1655.last.profile.txt

# take a look
head 2D_vs_MG1655.last.profile.txt
```

We can visualise these in R:

```bash
R
```
```R
prof <- read.table("2D_vs_MG1655.last.profile.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# plot alignment length against read length
plot(prof$read_len, prof$align_len, xlab="Read Length", ylab="Aln length", pch=16)

# find the maximum aligned length
max(prof$align_len)
prof[ prof$align_len == max(prof$align_len) , ]
```

In their [Gigascience paper](http://www.gigasciencejournal.com/content/3/1/22), Quick et al defined:

 Read percentage identity is defined as 100 * matches/(matches + deletions + insertions + mismatches)

We can calculate this in R:

```R
pid <- 100 * prof$matches / (prof$matches + prof$deletions + prof$insertions + prof$mismatches)

# show a boxplot
boxplot(pid, ylab="% ID", col="skyblue2", xlab="2D")
```

# Calling variants with nanopolish

Nanopolish was written by Jared Simpson and represents a suite of tools for aligning events to a reference sequence and working with those (using HMMs) to improve e.g. SNP calling and consensus accuracy

Nanopolish eventalign aligns the raw squiggle data directly to a reference genome, but requires a sequence-based alignment to guide it.

Once we have events aligned, we can do SNP calling using Nanopolish variants

```bash
nanopolish eventalign --reads MAP006-1.2D.fasta -b 2D_vs_MG1655.bwa.bam -g Data/reference/Ecoli_MG1655/MG1655.fa --sam \|
           samtools view -bS - \|
           samtools sort - -o MAP006-1.2D.eventalign.bam

samtools index MAP006-1.2D.eventalign.bam
```

We can then SNP call these with the variants module
```bash

```
