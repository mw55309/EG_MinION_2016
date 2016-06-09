# Edinburgh Genomics - MinION Bioinformatics
## Nanopore Data: FAST5, FASTA/Q, Events and Hairpins

## FAST5 Files 

MinION data files are FAST5 files, and FAST5 files are fundamentally [HDF5 files](https://www.hdfgroup.org/HDF5/whatishdf5.html). HDF5 supports storage of a wide variety of data types in the same file, and the format is portable, extensible and [widely supported] (https://www.hdfgroup.org/products/hdf5_tools/SWSummarybyName.htm). HDF5 files are binary files organized in a hierarchical, filesystem-like structure, with two primary object types: *groups* and *datasets*.

* Groups are container structures that can hold datasets and other groups
* Datasets are multidimensional arrays of data
* *Attributes* are user-defined data structures attached to groups and datasets, providing extra information about the HDF5 object to which they are atached.

In FAST5 files major data blocks like sequences are stored in datasets arranged in a tree-structure of groups. Metadata are stored in group and dataset attributes. The same file format is used pre- and post-basecalling, with groups, datasets and associated attributes added incrementally.

### A Brief History

Oxford Nanopore are very bad at releasing official definitions of file formats and there is consequently a significant amount of guess work and exploration required to work out where all the data is inside FAST5 files.

Most of the early ONT data was released from the SQK-MAP-005 kits - this includes 
[the MARC data](http://f1000research.com/articles/4-1075/v1), [Mick's B fragilis dataset](http://gigadb.org/dataset/100177) and [Nick Loman's first E coli dataset](http://gigadb.org/dataset/100102).  These data were encoded in what can be best described as FAST5 v.1.0 (ONT don't actually assign version numbers!)

Then SQK-MAP-006 came along, which was a major chemistry change that increased throughput.  A major change is that metrichor has switched from a 5mer model to a 6mer model. Late last year, [Nick released E coli SQK-MAP-006 data](http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/), and because he was very quick to do this, the files are still in FAST5 v.1.0. This is the main data set we will be working with today.

In November 2015, ONT released a new file format, which we can call FAST5 v1.1.  The major difference is that the template and complement FASTQ and events data have been moved to a new group within the FAST5 file, separate to the 2D data.  This actually makes logical sense, but moving things around can complicate data analysis.

In May 2016, ONT released the R9 pore and SQK-NSK007 kit (replacing the R7 pore and SQK-MAP-006). 

Inline-style: 
![alt text](https://2.bp.blogspot.com/-4Fb6ZFsqPkY/Vt-EBUVLsCI/AAAAAAAAhTs/kT0-tXx1s4A/s1600/2016-03-08-CsgG.PNG "Logo Title Text 1")

#### Major file format differences

The major difference is where the template and complement data are.  In version 1.0 they are all in a group called Basecall_2D_000; however, in v1.1 they have been moved to Basecall_1D_000

FAST5 v1.0
* /Analyses/Basecall_**2D**_000/BaseCalled_2D/
* /Analyses/Basecall_**2D**_000/BaseCalled_template/
* /Analyses/Basecall_**2D**_000/BaseCalled_complement/

FAST5 v1.1
* /Analyses/Basecall_**2D**_000/BaseCalled_2D/
* /Analyses/Basecall_**1D**_000/BaseCalled_template/
* /Analyses/Basecall_**1D**_000/BaseCalled_complement/


## File handling

The folder Data/read_data/MAP006-1_2100000-260000000_fast5/ contains some SQK-MAP-006 data from the Loman lab.  There are only a few thousand reads in there (to save time), but it is not uncommon to have folders containing 100,000s of FAST5 files, and therefore ls can be uninformative/not work.  However, we can use find to reveal directory structure (here using the top-level directory Data)

```sh
find Data -type d
```

## Using HDF5 command line tools

h5ls and h5dump can be quite useful.

h5ls reveals the structure of fast5 files.  Adding the -r flag makes this recursive:

```sh
h5ls Data/read_data/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch56_file159_strand.fast5
```
```
Analyses                 Group
Sequences                Group
UniqueGlobalKey          Group
```

Adding the -r flag makes this recursive
```sh
h5ls -r Data/read_data/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch56_file159_strand.fast5
```
```
/                        Group
/Analyses                Group
/Analyses/EventDetection_000 Group
/Analyses/EventDetection_000/Configuration Group
/Analyses/EventDetection_000/Configuration/abasic_detection Group
/Analyses/EventDetection_000/Configuration/event_detection Group
/Analyses/EventDetection_000/Configuration/hairpin_detection Group
/Analyses/EventDetection_000/Reads Group
/Analyses/EventDetection_000/Reads/Read_151 Group
/Analyses/EventDetection_000/Reads/Read_151/Events Dataset {13837/Inf}
/Sequences               Group
/Sequences/Meta          Group
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group
```

We can compare this to the corresponding base-called file

```sh
h5ls -r Data/read_data/MAP006-1_2100000-2600000_fast5/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch56_file159_strand.fast5
```
```
/                        Group
/Analyses                Group
/Analyses/Basecall_2D_000 Group
/Analyses/Basecall_2D_000/BaseCalled_2D Group
/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment Dataset {8429}
/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq Dataset {SCALAR}
/Analyses/Basecall_2D_000/BaseCalled_complement Group
/Analyses/Basecall_2D_000/BaseCalled_complement/Events Dataset {6677}
/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq Dataset {SCALAR}
/Analyses/Basecall_2D_000/BaseCalled_complement/Model Dataset {4096}
/Analyses/Basecall_2D_000/BaseCalled_template Group
/Analyses/Basecall_2D_000/BaseCalled_template/Events Dataset {7090}
/Analyses/Basecall_2D_000/BaseCalled_template/Fastq Dataset {SCALAR}
/Analyses/Basecall_2D_000/BaseCalled_template/Model Dataset {4096}
/Analyses/Basecall_2D_000/Configuration Group
/Analyses/Basecall_2D_000/Configuration/aggregator Group
/Analyses/Basecall_2D_000/Configuration/basecall_1d Group
/Analyses/Basecall_2D_000/Configuration/basecall_2d Group
/Analyses/Basecall_2D_000/Configuration/calibration_strand Group
/Analyses/Basecall_2D_000/Configuration/general Group
/Analyses/Basecall_2D_000/Configuration/hairpin_align Group
/Analyses/Basecall_2D_000/Configuration/post_processing Group
/Analyses/Basecall_2D_000/Configuration/post_processing.3000Hz Group
/Analyses/Basecall_2D_000/Configuration/post_processing.5000Hz Group
/Analyses/Basecall_2D_000/Configuration/recipes Group
/Analyses/Basecall_2D_000/Configuration/split_hairpin Group
/Analyses/Basecall_2D_000/HairpinAlign Group
/Analyses/Basecall_2D_000/HairpinAlign/Alignment Dataset {6261}
/Analyses/Basecall_2D_000/InputEvents Soft Link {Analyses/EventDetection_000/Reads/Read_151/Events}
/Analyses/Basecall_2D_000/Log Dataset {SCALAR}
/Analyses/Basecall_2D_000/Summary Group
/Analyses/Basecall_2D_000/Summary/basecall_1d_complement Group
/Analyses/Basecall_2D_000/Summary/basecall_1d_template Group
/Analyses/Basecall_2D_000/Summary/basecall_2d Group
/Analyses/Basecall_2D_000/Summary/hairpin_align Group
/Analyses/Basecall_2D_000/Summary/post_process_complement Group
/Analyses/Basecall_2D_000/Summary/post_process_template Group
/Analyses/Basecall_2D_000/Summary/split_hairpin Group
/Analyses/Calibration_Strand_000 Group
/Analyses/Calibration_Strand_000/Configuration Group
/Analyses/Calibration_Strand_000/Configuration/aggregator Group
/Analyses/Calibration_Strand_000/Configuration/basecall_1d Group
/Analyses/Calibration_Strand_000/Configuration/basecall_2d Group
/Analyses/Calibration_Strand_000/Configuration/calibration_strand Group
/Analyses/Calibration_Strand_000/Configuration/general Group
/Analyses/Calibration_Strand_000/Configuration/hairpin_align Group
/Analyses/Calibration_Strand_000/Configuration/post_processing Group
/Analyses/Calibration_Strand_000/Configuration/post_processing.3000Hz Group
/Analyses/Calibration_Strand_000/Configuration/post_processing.5000Hz Group
/Analyses/Calibration_Strand_000/Configuration/recipes Group
/Analyses/Calibration_Strand_000/Configuration/split_hairpin Group
/Analyses/Calibration_Strand_000/InputEvents Soft Link {Analyses/EventDetection_000/Reads/Read_151/Events}
/Analyses/Calibration_Strand_000/Log Dataset {SCALAR}
/Analyses/Calibration_Strand_000/Summary Group
/Analyses/Calibration_Strand_000/Summary/calibration_strand_2d Group
/Analyses/Calibration_Strand_000/Summary/calibration_strand_complement Group
/Analyses/Calibration_Strand_000/Summary/calibration_strand_template Group
/Analyses/EventDetection_000 Group
/Analyses/EventDetection_000/Configuration Group
/Analyses/EventDetection_000/Configuration/abasic_detection Group
/Analyses/EventDetection_000/Configuration/event_detection Group
/Analyses/EventDetection_000/Configuration/hairpin_detection Group
/Analyses/EventDetection_000/Reads Group
/Analyses/EventDetection_000/Reads/Read_151 Group
/Analyses/EventDetection_000/Reads/Read_151/Events Dataset {13837}
/Sequences               Group
/Sequences/Meta          Group
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group
```

We can also use h5ls -d to extract specific datasets from the files (note the /Analyses/Basecall_2D_000/BaseCalled_2D/Fastq appended to the filename):
```sh
 h5ls -d Data/read_data/MAP006-1_2100000-2600000_fast5/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch56_file159_strand.fast5/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq
 ```
 ```
Fastq                    Dataset {SCALAR}
    Data:
        (0) "@1b63c36c-fb28-4cf3-8df1-3eb075eb00b2_Basecall_2D_000_2d LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch56_file159_strand\nCATTTTCTTTCTTACTGATATTAGTTTTTGG
   ...
   <snip>
   ...
   ,*)))+,))*)***++*))**+)))\n"
```
(but beware the odd " '+' repeats 8 times " in the quality socres!


Unsurprisingly h5dump dumps the entire file to STDOUT

```sh
h5dump Data/read_data/MAP006-1_2100000-2600000_fast5/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch56_file159_strand.fast5
```

## Browsing HDF5 files

Any HDF5 file can be opened using hdfview and browsed/edited in a GUI

```sh
# non base-called
hdfview Data/read_data/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch56_file159_strand.fast5 &

# base called
hdfview Data/read_data/MAP006-1_2100000-2600000_fast5/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch56_file159_strand.fast5 &
```

## Basic manipulation in poRe

poRe is a library for R available from [SourceForge](https://sourceforge.net/projects/rpore/) and published in [bioinformatics](http://bioinformatics.oxfordjournals.org/content/31/1/114).  poRe is incredibly simple to install and relies simply on R 3.0 or above and a few additional libraries.

The poRe library is set up to read v1.1 data by default, and offers users parameters to enable reading of v1.0 data.  Let's start it up.

```sh
R
```
```R
library(poRe)
```

#### FASTQ

We'll see how to extract FASTQ from entire directories below, but here are some exemples of single file analysis

```R
f5 <- "Data/read_data/MAP006-1_2100000-2600000_fast5/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch56_file159_strand.fast5"
get_fastq(f5)
```

This returns a list with a single fastq datasets, "2D".  However, where are template and complement?  Unfortunately the defaults are set to v1.1 and this file is v1.0.  We can use the path.t and path.c parameters to tell poRe where the template and complement data are

```R
get_fastq(f5, path.t="/Analyses/Basecall_2D_000/", path.c="/Analyses/Basecall_2D_000/")
fqs <- get_fastq(f5, path.t="/Analyses/Basecall_2D_000/", path.c="/Analyses/Basecall_2D_000/")
names(fqs)
```

If we don't want to extract all 3, we can choose which to extract using the "which" argument

```R
get_fastq(f5, path.t="/Analyses/Basecall_2D_000/", which="template")
```

Working with lists is easy in R, and if you want the FASTQ as a string:

```R
fq <- get_fastq(f5, path.t="/Analyses/Basecall_2D_000/", which="all")
names(fq)
fq$template
```

From here you can see that it's incredibly simple to write a FASTQ extraction script:

```R
# path to fast5 files
f5dir <- "Data/read_data/small_test"
# get vector of all fast5 files
f5files <- dir(f5dir, pattern="\\.fast5$", full.names = TRUE)
# iterate over files
for (f5 in f5files) {
    # extract 2D fastq
    fq <- get_fastq(f5, which="2D")
    # check fq is a list and contains 2D
    if (typeof(fq) == "list" && exists("2D", where=fq)) {
        # cat to "" (STDOUT but could be the name of a file
        # change the cat = "" to a filename to see what happens
        cat(fq[["2D"]], file = "", sep = "\n", fill = FALSE)
    }
}
```

However, we have scripts that do this already, so there is no need to write your own!

#### FASTA data

Simply use get_fasta instead of get_fastq :-)

#### Events data

We extract events data with the get.events function.  As events data are co-located with FASTQ in the FAST5 file, then for Nick's SQK-MAP-006 data we need to give it the paths again.

get.events again returns a list, wth the template and complement events as data frames

```R
f5 <- "Data/read_data/MAP006-1_2100000-2600000_fast5/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch56_file159_strand.fast5"
ev <- get.events(f5, path.t = "/Analyses/Basecall_2D_000/", path.c = "/Analyses/Basecall_2D_000/")
names(ev)
head(ev$template)
head(ev$complement)
```

* start: time in seconds
* mean: the mean signal of the event
* stdev: the standard deviation when sampling
* length: length of the event
* model_state: the predicted kmer
* model_level: the value from the model
* move: how many moves the model has made in that step

#### Extracting the model

Extracting the model itself in poRe is also easy, but relies on a little legacy code that needs to be updated... next version ;-)

```R
mods <- get.models(f5, tmodel = "/Analyses/Basecall_2D_000/BaseCalled_template/Model", 
                        cmodel = "/Analyses/Basecall_2D_000/BaseCalled_complement/Model")
head(mods$template)
head(mods$complement)
```

Here we can see the models that ONT use to turn events into kmers, but this isn't simply a case of looking up the closest match:

```R
head(ev$template)
# first event kmer is GTGTTT, event mean is 50.23671, model mean is 50.36703
     mean    start      stdv      length model_state model_level move
1 50.23671 65351.52 0.6342780 0.049136786      GTGTTT    50.36703    0

mods$template[mods$template$kmer=="GTGTTT",]

# however the model level for this kmer is quite different:
       kmer level_mean level_stdv  sd_mean  sd_stdv   weight
3008 GTGTTT   43.29339    0.55339 0.78968 0.267732 623.0312
```

What gives?  Well, ONT apply three parameters to the model to "scale" the model to fit each read.  These parameters are drift, scale and shift and they are read dependent!

#### Extracting the model parameters

For v1.0:

```R
get.model.params(f5, tsum="/Analyses/Basecall_2D_000/Summary/basecall_1d_template", 
                     csum = "/Analyses/Basecall_2D_000/Summary/basecall_1d_complement")
```

For v1.1:
```R
get.model.params(f5)
```

For this template read, we have:
* drift	0.0003234845
* scale	0.9643541
* shift	8.616875

#### Applying the model parameters

The drift parameter is applied per event and is scaled by the number of seconds since the start of the read.  So "drift * (time - min(time))" can be subtracted from the event mean, or added to the model mean.  As we are dealing with the first event then the drift parameter isn't applied.

Scale and shift are then applied in a classic linear model: (model.mean * scale) + shift.  In our case:

* (43.29339 * 0.9643541) + 8.616875 = 50.36703

Which is the model value that shows up in the events table above. 

#### Plotting Squiggles

Plotting squiggles works directly from the events data.  The parameters minseconds and maxseconds control which events to plot:

```R
plot.squiggle(ev$template)
plot.squiggle(ev$template, minseconds=5)
plot.squiggle(ev$complement)
```

## Run QC in poRe

If you ran pore_rt() during your nanopore run then you will have access to a metadata text file that can be used for run QC etc.  Otherwise we will have to create one (more below).  In the meantime, we have meta data for the pass and fail folders for Nick's SQK-MAP-006 run

```R
# load in pass data
pass <- read.table("Data/run_metadata/pass.meta.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# set standard/expected column names
colnames(pass) <- c("filename","channel_num","read_num","read_start_time",
                    "status","tlen","clen","len2d",
                    "run_id","read_id","barcode","exp_start")

# load in the fail data
fail <- read.table("Data/run_metadata/fail.meta.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# set standard/expected column names
colnames(fail) <- c("filename","channel_num","read_num","read_start_time",
                    "status","tlen","clen","len2d",
                    "run_id","read_id","barcode","exp_start")
head(pass)
head(fail)
```

### Longest high quality read

We can find this from the metadata:

```R
# find the maximum length
max(pass$len2d)

# get the metadata for that read
pass[pass$len2d==max(pass$len2d),]

# We now know the longest read
longest <- "Data/read_data/MAP006-1_2100000-2600000_fast5/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch153_file57_strand.fast5"
lfq <- get_fastq(longest, which="2D")
lfa <- get_fasta(longest, which="2D")

# write fastq and fasta out to files
cat(lfq[["2D"]], file = "longest.fastq", sep = "\n", fill = FALSE)
cat(lfa[["2D"]], file = "longest.fasta", sep = "\n", fill = FALSE)
```


### Yield

Yield over time can be plotted with plot.cumulative.yield

```R
yield.p <- plot.cumulative.yield(pass)
yield.f <- plot.cumulative.yield(fail)
```

The calculated cumulative yields are returned as data.frames

```R
head(yield.p)
head(yield.f)
```

### Read length histogram

We can plot read lengths histograms

```R
plot.length.histogram(pass)
plot.length.histogram(fail)
```

There's quite a long failed template read!

```R
max(fail$tlen)
# [1] 379692
fail [fail$tlen==379692,]
```

Plotting in ggplot2 if you really want to

```R
library(ggplot2)
m <- ggplot(pass, aes(x=len2d))
m + geom_histogram(binwidth=500)
```

### Channel and pore occupancy

The MinION flowcell is arranged into 512 channels in 4 blocks, and we can see the layout using poRe:

```R
show.layout()
```

We first calculate some statistics summarised by channel:

```R
# pass data
pass.s <- summarise.by.channel(pass)
head(pass.s)

fail.s <- summarise.by.channel(fail)
head(fail.s)
```

The rows of the result are the channel numbers, and the columns tell us how many channels appear in our summary data, and either the number (n) or cumulative length (l) of template, complement and 2d reads from each channel.

We can plot these:

```R
# pass

# the number of times the channel appears in any context
plot.channel.summary(pass.s)

# cumulative 2D length
plot.channel.summary(pass.s, report.col="l2d")

# fail

# the number of times the channel appears in any context
plot.channel.summary(fail.s)

# cumulative 2D length
plot.channel.summary(fail.s, report.col="l2d")

```

### Extracting meta-data from fast5
If you haven't run pore_rt(), then you can extract meta-data directly from the fast5 files.  This takes a long time as we have to open each file and extract the attributes.  For the 2666 files we have selected here,  and on this VM, it takes about 4 minutes.

```R
pass <- read.meta.info("Data/read_data/MAP006-1_2100000-2600000_fast5/", 
                        path.t="/Analyses/Basecall_2D_000/", 
                        path.c="/Analyses/Basecall_2D_000/", 
                        pass="P")
                        
plot.cumulative.yield(pass)
```

## Extracting FASTQ from the command-line


Command-line scripts for extracting FASTQ can be pulled from [github](https://github.com/mw55309/poRe_scripts).  There are scripts there that can parse both the old and new format of fast5 file, but we have installed the old format scripts here:

```sh
# 2D
extract2D Data/read_data/MAP006-1_2100000-2600000_fast5/ > MAP006-1.2D.fastq

# template
extractTemplate Data/read_data/MAP006-1_2100000-2600000_fast5/ > MAP006-1.template.fastq

# complement
extractComplement Data/read_data/MAP006-1_2100000-2600000_fast5/ > MAP006-1.complement.fastq
```

FASTQ can be converted to FASTA using a small script

```sh
porefq2fa MAP006-1.2D.fastq > MAP006-1.2D.fasta
```

## poretools

Much of what we did above can also be achieved using poretools.  We have cheated slightly with poRe by extracting most of the metadata we need during the run and metrichor base-calling, which allows us to do things much more quickly.  Below, poretools is going to each fast5 file and ripping out the necessary information which takes a long time.  poRe also takes a long time when having to do this, which is why we chose to pre-extract the necessary meta-data.

For FASTQ extraction, we find that poretools is quicker on pass data, but poRe is quicker on fail data - we're not sure why!

```sh
# extract 2D fastq
poretools fastq --type 2D  Data/read_data/MAP006-1_2100000-2600000_fast5 > MAP006-1.2D.poretools.fastq

# extract 2d fasta
poretools fasta --type 2D Data/read_data/MAP006-1_2100000-2600000_fast5 > MAP006-1.2D.poretools.fasta

# read length histogram
# this doesn't work because we had trouble installing rpy2
poretools hist Data/read_data/MAP006-1_2100000-2600000_fast5

# yield plot
# this doesn't work because we had trouble installing rpy2
poretools yield_plot Data/read_data/MAP006-1_2100000-2600000_fast5

# occupancy plot
# this doesn't work because we had trouble installing rpy2
poretools occupancy Data/read_data/MAP006-1_2100000-2600000_fast5

# find the longest
poretools winner Data/read_data/MAP006-1_2100000-2600000_fast5
```

