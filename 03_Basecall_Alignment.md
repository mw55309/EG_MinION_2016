# Edinburgh Genomics - MinION Bioinformatics

## Alternative Basecallers

Up until recently, the only basecaller for MinION was the ONT cloud-based basecaller EPI2ME/Metrichor. However, two alternative basecallers were recently developed in the community: [nanocall](http://biorxiv.org/content/early/2016/03/28/046086) and [deepnano](http://arxiv.org/abs/1603.09195). They were released almost simultaneously. We are going to take a quick look at each of them here.

### Nanocall

[Nanocall](https://github.com/mateidavid/nanocall) was developed in Jared Simpson's lab, and performs basecalling using a hidden markov model (like the original Metrichor). For reasons of time, we have put just a couple of pre-basecalled files in ~/Data/pre_basecall. We can basecall them with nanocall like this:

```sh
nanocall/bin/nanocall -t 2 ~/Data/pre_basecall >~/Data/pre_basecall/nanocall_output.fa 2>~/Data/pre_basecall/nanocall.log
```

Nanocall just does 1D basecalling at the moment. In the ouptput file, you will find template and complement basecalls for both files.

### Deepnano

[Deepnano](https://bitbucket.org/vboza/deepnano) uses a deep recurrent neural network (RNN) model to call bases. ONT is also moving towards RNN basecalling and it is likely this will be the main method for calling R9 (and future data) going forward. Deepnano is a python package built on the Theano framework. We can run it on a directory of FAST5 files like this:

```sh
python ~/deepnano/basecall_no_metrichor.py --directory ~/Data/pre_basecall --output ~/Data/pre_basecall/deepnano.fasta
```

## Read alignment

## Aligning reads to a reference

A number of different alignment tools have been proposed for nanopore data. Three good ones are [BWA](http://bio-bwa.sourceforge.net/), [LAST](http://last.cbrc.jp/) and [Graphmap](https://github.com/isovic/graphmap). We will start with BWA.

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

Finally, Graphmap:

```sh
~/graphmap/bin/Linux-x64/graphmap align -t 2 -r ~/Data/reference/Ecoli_MG1655/MG1655.fa -d MAP006-1.2D.fastq -o gm.sam

samtools view -Sb gm.sam | \
samtools sort -o 2D_vs_MG1655.graphmap.bam - 

# index bam
samtools index 2D_vs_MG1655.graphmap.bam
```

## Viewing Alignments in IGV

You can read about viewing alignments in IGV [here](https://www.broadinstitute.org/software/igv/AlignmentData). We start by running 

```sh
igv.sh
```

If required, we load a genome from file ~/Data/reference/Ecoli_MG1655/MG1655.fa, and after that we can load one or more 2D_vs_MG1655.<aligner>.bam files to see how they compare.

## Assessing the accuracy of alignment
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

## Aligning events with Nanopolish EventAlign

Nanopolish was written by Jared Simpson and represents a suite of tools for aligning events to a reference sequence and working with those (using HMMs) to improve e.g. SNP calling and consensus accuracy

Nanopolish eventalign aligns the raw squiggle data directly to a reference genome, but requires a sequence-based alignment to guide it.

**If you run the command below it can take up to ~1hr on the VM we're working on, so I recommend you don't**

```bash
nanopolish eventalign --print-read-names --scale-events --reads MAP006-1.2D.fasta \
                      -b 2D_vs_MG1655.bwa.bam -g Data/reference/Ecoli_MG1655/MG1655.fa \
                      > Data/eventalign/MAP006-1.2D.eventalign.txt

# eventalign can also output SAM format with the --sam flag
```

We have a smaller dataset pre-run in Data/eventalign/singleread.eventalign.txt so let's take a look at that:

```bash
head Data/eventalign/singleread.eventalign.txt
```

Now we can visualise in R

```bash
R
```
```R
# Load the poRe library for the plot.as.squiggle function...
library(poRe)

# read in eventalign output as a data.frame
d <- read.table("Data/eventalign/singleread.eventalign.txt", header=TRUE, sep="\t")

# let's just look at the template strand
dt <- d[d$strand=="t",]

# plot.as.squiggle is a utility function in R that takes a vector of
# numbers and plots them as if they were an event-driven nanopore
# squiggle.  We need this because in real nanopore data, the events
# are differing times apart, so can be hard to visualise, whereas
# plot.as.squiggle simply uses the same false time interval

# get data for the first 100 event means
em <- plot.as.squiggle(dt$event_level_mean[1:100], plot=FALSE)

# get data for the first 100 model means they are aligned to
mm <- plot.as.squiggle(dt$model_mean[1:100], plot=FALSE)

# plot and visualise
plot(em$times, em$means, type="l", ylab="mean signal", xlab="fake time")
lines(mm$times, mm$means, col="red")
```
## Hairpins

Finally, we will take a look at hairpins. Hairpins divide the template and complement components of a MinION read. The hairpins form a characteristic current signal that is detected at the start of the basecalling process and accurate detection is crucial to accurate basecalling. 

As I will show in a demo, they look rather different in R7 and R9, and are not always easy to detect correctly.
