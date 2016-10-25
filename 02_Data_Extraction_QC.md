# Edinburgh Genomics - MinION Bioinformatics
## Nanopore Data: FAST5, FASTA/Q, Events and Metadata

## FAST5 Files 

MinION data files are FAST5 files, and FAST5 files are fundamentally [HDF5 files](https://www.hdfgroup.org/HDF5/whatishdf5.html). HDF5 supports storage of a wide variety of data types in the same file, and the format is portable, extensible and [widely supported] (https://www.hdfgroup.org/products/hdf5_tools/SWSummarybyName.htm). HDF5 files are binary files organized in a hierarchical, filesystem-like structure, with two primary object types: *groups* and *datasets*.

* Groups are container structures that can hold datasets and other groups
* Datasets are multidimensional arrays of data
* *Attributes* are user-defined data structures attached to groups and datasets, providing extra information about the HDF5 object to which they are atached.

In FAST5 files major data blocks like FASTQ sequences are stored in datasets arranged in a tree-structure of groups. Metadata are stored in group and dataset attributes. The same file format is used pre- and post-basecalling, with groups, datasets and associated attributes added incrementally.

### A Brief History of FAST5

Oxford Nanopore are very bad at releasing official definitions of file formats and there is consequently a significant amount of guess work and exploration required to work out where all the data is inside FAST5 files.

Most of the early ONT data was released from the SQK-MAP-005 kits - this includes 
[the MARC data](http://f1000research.com/articles/4-1075/v1), [Mick's B fragilis dataset](http://gigadb.org/dataset/100177) and [Nick Loman's first E coli dataset](http://gigadb.org/dataset/100102).  These data were encoded in what can be best described as FAST5 v1.0 (ONT don't actually assign version numbers!)

In early autumn 2015, SQK-MAP-006 came along. This was a major chemistry change that increased throughput, and one important difference was that metrichor switched from a 5mer model to a 6mer model. Early MAP-006 data files were stil in FAST5 v1.0. An example dataset of this type is available here: http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/

In November 2015, ONT released a new file format, which we can call FAST5 v1.1. The major difference was that the template and complement FASTQ and events data were moved to a new group within the FAST5 file, separate from the 2D data.  This actually makes logical sense, but moving things around can complicate data analysis.

#### Major file format differences (v1.0 vs v1.1)

The major difference is where the template and complement data are. In version 1.0 they are all in a group called Basecall_2D_000; however, in v1.1 they have been moved to Basecall_1D_000

FAST5 v1.0
* /Analyses/Basecall_**2D**_000/BaseCalled_2D/
* /Analyses/Basecall_**2D**_000/BaseCalled_template/
* /Analyses/Basecall_**2D**_000/BaseCalled_complement/

FAST5 v1.1
* /Analyses/Basecall_**2D**_000/BaseCalled_2D/
* /Analyses/Basecall_**1D**_000/BaseCalled_template/
* /Analyses/Basecall_**1D**_000/BaseCalled_complement/

####R9

In May 2016, ONT released the R9 pore and SQK-NSK007 kit (replacing the R7 pore and SQK-MAP-006):

![](https://2.bp.blogspot.com/-4Fb6ZFsqPkY/Vt-EBUVLsCI/AAAAAAAAhTs/kT0-tXx1s4A/s1600/2016-03-08-CsgG.PNG "The R9 Pore")

R9 offers greater accuracy and sequences at higher speeds. A major new feature of R9 data files is that they include Raw data in a new top-level group (/Raw).

####Rapid 1D

Due to the noisiness of R7 nanopore reads, virtually all experiments prior to R9 used 2D library prep where reads of both template and complement stands are combined for greater accuracy through consensus. The greater accuracy of R9 means that '1D' reads (template strand only) are sufficiently accurate for many applications, and a transposase-based rapid library prep protocol enables libraries to be prepared in around 5 minutes.

Rapid 1D FAST5 files naturally lack complement and 2D groups. They are otherwise similar to other R9 files.

Today, we will mostly be working with an R9 1D (E coli) dataset released by Nick Loman in July 2016: http://lab.loman.net/2016/07/30/nanopore-r9-data-release/

## File handling

The folder Data/read_data/R9_1D_FAST5 contains a few hundred reads (to save time here), but it is not uncommon to have folders containing 100,000s of FAST5 files, and therefore the ls command can be uninformative/not work.  However, we can use find to reveal directory structure (here using the top-level directory Data)

```sh
find ~/Data -type d
```

## Using HDF5 command line tools

h5ls and h5dump can be quite useful.

h5ls reveals the top-level structure of fast5 files. Here is a 2D R7 file:

```sh
h5ls ~/Data/read_data/r7_2d_ecoli_ch56_file159.fast5
```
```
Analyses                 Group
Sequences                Group
UniqueGlobalKey          Group
```

And here is a 1D rapid R9 example:

```sh
h5ls ~/Data/read_data/r9_1d_ecoli_ch2_read16.fast5
```
```
Analyses                 Group
Raw                      Group
UniqueGlobalKey          Group
```

Adding the -r flag makes the listing recursive. This is a pre-basecalled 2D R9 file (from a zika dataset):
```sh
h5ls -r ~/Data/read_data/r9_2d_zika_ch1_read10_pre.fast5
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
/Analyses/EventDetection_000/Reads/Read_10 Group
/Analyses/EventDetection_000/Reads/Read_10/Events Dataset {2176/Inf}
/Raw                     Group
/Raw/Reads               Group
/Raw/Reads/Read_10       Group
/Raw/Reads/Read_10/Signal Dataset {50722/Inf}
/Sequences               Group
/Sequences/Meta          Group
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group
```

We can compare this to the corresponding (2D R9) base-called file. Note all the basecalling data added into the /Analyses/Basecall_1D_000 and /Analyses/Basecall_2D_000 groups.

```sh
h5ls -r ~/Data/read_data/r9_2d_zika_ch1_read10.fast5
```
```
/                        Group
/Analyses                Group
/Analyses/Basecall_1D_000 Group
/Analyses/Basecall_1D_000/BaseCalled_complement Group
/Analyses/Basecall_1D_000/BaseCalled_complement/Events Dataset {2006}
/Analyses/Basecall_1D_000/BaseCalled_complement/Fastq Dataset {SCALAR}
/Analyses/Basecall_1D_000/BaseCalled_template Group
/Analyses/Basecall_1D_000/BaseCalled_template/Events Dataset {2538}
/Analyses/Basecall_1D_000/BaseCalled_template/Fastq Dataset {SCALAR}
/Analyses/Basecall_1D_000/Configuration Group
/Analyses/Basecall_1D_000/Configuration/aggregator Group
/Analyses/Basecall_1D_000/Configuration/basecall_1d Group
/Analyses/Basecall_1D_000/Configuration/basecall_2d Group
/Analyses/Basecall_1D_000/Configuration/calibration_strand Group
/Analyses/Basecall_1D_000/Configuration/components Group
/Analyses/Basecall_1D_000/Configuration/event_detection Group
/Analyses/Basecall_1D_000/Configuration/general Group
/Analyses/Basecall_1D_000/Configuration/hairpin_align Group
/Analyses/Basecall_1D_000/Configuration/post_processing Group
/Analyses/Basecall_1D_000/Configuration/post_processing.4000Hz Group
/Analyses/Basecall_1D_000/Configuration/split_hairpin Group
/Analyses/Basecall_1D_000/Log Dataset {SCALAR}
/Analyses/Basecall_1D_000/Summary Group
/Analyses/Basecall_1D_000/Summary/basecall_1d_complement Group
/Analyses/Basecall_1D_000/Summary/basecall_1d_template Group
/Analyses/Basecall_2D_000 Group
/Analyses/Basecall_2D_000/BaseCalled_2D Group
/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment Dataset {3111}
/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq Dataset {SCALAR}
/Analyses/Basecall_2D_000/Configuration Group
/Analyses/Basecall_2D_000/Configuration/aggregator Group
/Analyses/Basecall_2D_000/Configuration/basecall_1d Group
/Analyses/Basecall_2D_000/Configuration/basecall_2d Group
/Analyses/Basecall_2D_000/Configuration/calibration_strand Group
/Analyses/Basecall_2D_000/Configuration/components Group
/Analyses/Basecall_2D_000/Configuration/event_detection Group
/Analyses/Basecall_2D_000/Configuration/general Group
/Analyses/Basecall_2D_000/Configuration/hairpin_align Group
/Analyses/Basecall_2D_000/Configuration/post_processing Group
/Analyses/Basecall_2D_000/Configuration/post_processing.4000Hz Group
/Analyses/Basecall_2D_000/Configuration/split_hairpin Group
/Analyses/Basecall_2D_000/HairpinAlign Group
/Analyses/Basecall_2D_000/HairpinAlign/Alignment Dataset {2177}
/Analyses/Basecall_2D_000/Log Dataset {SCALAR}
/Analyses/Basecall_2D_000/Summary Group
/Analyses/Basecall_2D_000/Summary/basecall_2d Group
/Analyses/Basecall_2D_000/Summary/hairpin_align Group
/Analyses/Basecall_2D_000/Summary/post_process_complement Group
/Analyses/Basecall_2D_000/Summary/post_process_template Group
/Analyses/Calibration_Strand_000 Group
/Analyses/Calibration_Strand_000/Configuration Group
/Analyses/Calibration_Strand_000/Configuration/aggregator Group
/Analyses/Calibration_Strand_000/Configuration/basecall_1d Group
/Analyses/Calibration_Strand_000/Configuration/basecall_2d Group
/Analyses/Calibration_Strand_000/Configuration/calibration_strand Group
/Analyses/Calibration_Strand_000/Configuration/components Group
/Analyses/Calibration_Strand_000/Configuration/general Group
/Analyses/Calibration_Strand_000/Configuration/genome_mapping Group
/Analyses/Calibration_Strand_000/Configuration/hairpin_align Group
/Analyses/Calibration_Strand_000/Configuration/post_processing.3000Hz Group
/Analyses/Calibration_Strand_000/Configuration/split_hairpin Group
/Analyses/Calibration_Strand_000/Log Dataset {SCALAR}
/Analyses/Calibration_Strand_000/Summary Group
/Analyses/EventDetection_000 Group
/Analyses/EventDetection_000/Configuration Group
/Analyses/EventDetection_000/Configuration/abasic_detection Group
/Analyses/EventDetection_000/Configuration/event_detection Group
/Analyses/EventDetection_000/Configuration/hairpin_detection Group
/Analyses/EventDetection_000/Reads Group
/Analyses/EventDetection_000/Reads/Read_10 Group
/Analyses/EventDetection_000/Reads/Read_10/Events Dataset {2176}
/Analyses/EventDetection_001 Group
/Analyses/EventDetection_001/Configuration Group
/Analyses/EventDetection_001/Configuration/aggregator Group
/Analyses/EventDetection_001/Configuration/basecall_1d Group
/Analyses/EventDetection_001/Configuration/basecall_2d Group
/Analyses/EventDetection_001/Configuration/calibration_strand Group
/Analyses/EventDetection_001/Configuration/components Group
/Analyses/EventDetection_001/Configuration/event_detection Group
/Analyses/EventDetection_001/Configuration/general Group
/Analyses/EventDetection_001/Configuration/hairpin_align Group
/Analyses/EventDetection_001/Configuration/post_processing Group
/Analyses/EventDetection_001/Configuration/post_processing.4000Hz Group
/Analyses/EventDetection_001/Configuration/split_hairpin Group
/Analyses/EventDetection_001/Log Dataset {SCALAR}
/Analyses/EventDetection_001/Reads Group
/Analyses/EventDetection_001/Reads/Read_10 Group
/Analyses/EventDetection_001/Reads/Read_10/Events Dataset {4563}
/Analyses/EventDetection_001/Summary Group
/Analyses/EventDetection_001/Summary/event_detection Group
/Analyses/Hairpin_Split_000 Group
/Analyses/Hairpin_Split_000/Configuration Group
/Analyses/Hairpin_Split_000/Configuration/aggregator Group
/Analyses/Hairpin_Split_000/Configuration/basecall_1d Group
/Analyses/Hairpin_Split_000/Configuration/basecall_2d Group
/Analyses/Hairpin_Split_000/Configuration/calibration_strand Group
/Analyses/Hairpin_Split_000/Configuration/components Group
/Analyses/Hairpin_Split_000/Configuration/event_detection Group
/Analyses/Hairpin_Split_000/Configuration/general Group
/Analyses/Hairpin_Split_000/Configuration/hairpin_align Group
/Analyses/Hairpin_Split_000/Configuration/post_processing Group
/Analyses/Hairpin_Split_000/Configuration/post_processing.4000Hz Group
/Analyses/Hairpin_Split_000/Configuration/split_hairpin Group
/Analyses/Hairpin_Split_000/Log Dataset {SCALAR}
/Analyses/Hairpin_Split_000/Summary Group
/Analyses/Hairpin_Split_000/Summary/split_hairpin Group
/Raw                     Group
/Raw/Reads               Group
/Raw/Reads/Read_10       Group
/Raw/Reads/Read_10/Signal Dataset {50722/Inf}
/Sequences               Group
/Sequences/Meta          Group
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group
```

Finally, here is the recursive listing for a recent 1D R9 file. Note all the changes:

```sh
h5ls -r ~/Data/read_data/r9_1d_ecoli_ch2_read16.fast5
```
```
/                        Group
/Analyses                Group
/Analyses/Basecall_1D_000 Group
/Analyses/Basecall_1D_000/BaseCalled_template Group
/Analyses/Basecall_1D_000/BaseCalled_template/Events Dataset {17733}
/Analyses/Basecall_1D_000/BaseCalled_template/Fastq Dataset {SCALAR}
/Analyses/Basecall_1D_000/Configuration Group
/Analyses/Basecall_1D_000/Configuration/aggregator Group
/Analyses/Basecall_1D_000/Configuration/basecall_1d Group
/Analyses/Basecall_1D_000/Configuration/calibration_strand Group
/Analyses/Basecall_1D_000/Configuration/components Group
/Analyses/Basecall_1D_000/Configuration/event_detection Group
/Analyses/Basecall_1D_000/Configuration/general Group
/Analyses/Basecall_1D_000/Configuration/split_hairpin Group
/Analyses/Basecall_1D_000/Log Dataset {SCALAR}
/Analyses/Basecall_1D_000/Summary Group
/Analyses/Basecall_1D_000/Summary/basecall_1d_template Group
/Analyses/Calibration_Strand_000 Group
/Analyses/Calibration_Strand_000/Configuration Group
/Analyses/Calibration_Strand_000/Configuration/aggregator Group
/Analyses/Calibration_Strand_000/Configuration/basecall_1d Group
/Analyses/Calibration_Strand_000/Configuration/basecall_2d Group
/Analyses/Calibration_Strand_000/Configuration/calibration_strand Group
/Analyses/Calibration_Strand_000/Configuration/components Group
/Analyses/Calibration_Strand_000/Configuration/general Group
/Analyses/Calibration_Strand_000/Configuration/genome_mapping Group
/Analyses/Calibration_Strand_000/Configuration/hairpin_align Group
/Analyses/Calibration_Strand_000/Configuration/post_processing.3000Hz Group
/Analyses/Calibration_Strand_000/Configuration/split_hairpin Group
/Analyses/Calibration_Strand_000/Log Dataset {SCALAR}
/Analyses/Calibration_Strand_000/Summary Group
/Analyses/EventDetection_000 Group
/Analyses/EventDetection_000/Configuration Group
/Analyses/EventDetection_000/Configuration/aggregator Group
/Analyses/EventDetection_000/Configuration/basecall_1d Group
/Analyses/EventDetection_000/Configuration/basecall_2d Group
/Analyses/EventDetection_000/Configuration/calibration_strand Group
/Analyses/EventDetection_000/Configuration/components Group
/Analyses/EventDetection_000/Configuration/event_detection Group
/Analyses/EventDetection_000/Configuration/general Group
/Analyses/EventDetection_000/Configuration/hairpin_align Group
/Analyses/EventDetection_000/Configuration/post_processing Group
/Analyses/EventDetection_000/Configuration/post_processing.4000Hz Group
/Analyses/EventDetection_000/Configuration/split_hairpin Group
/Analyses/EventDetection_000/Log Dataset {SCALAR}
/Analyses/EventDetection_000/Reads Group
/Analyses/EventDetection_000/Reads/Read_16 Group
/Analyses/EventDetection_000/Reads/Read_16/Events Dataset {18469}
/Analyses/EventDetection_000/Summary Group
/Analyses/EventDetection_000/Summary/event_detection Group
/Analyses/Segment_Linear_000 Group
/Analyses/Segment_Linear_000/Configuration Group
/Analyses/Segment_Linear_000/Configuration/aggregator Group
/Analyses/Segment_Linear_000/Configuration/basecall_1d Group
/Analyses/Segment_Linear_000/Configuration/calibration_strand Group
/Analyses/Segment_Linear_000/Configuration/components Group
/Analyses/Segment_Linear_000/Configuration/general Group
/Analyses/Segment_Linear_000/Configuration/split_hairpin Group
/Analyses/Segment_Linear_000/Log Dataset {SCALAR}
/Analyses/Segment_Linear_000/Summary Group
/Analyses/Segment_Linear_000/Summary/split_hairpin Group
/Raw                     Group
/Raw/Reads               Group
/Raw/Reads/Read_16       Group
/Raw/Reads/Read_16/Signal Dataset {177981/Inf}
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group
```

We can also use h5ls -d to extract specific datasets from the files (note the /Analyses/Basecall_1D_000/BaseCalled_template/Fastq appended to the filename):
```sh
h5ls -d ~/Data/read_data/r9_1d_ecoli_ch2_read16.fast5/Analyses/Basecall_1D_000/BaseCalled_template/Fastq
 ```
 ```
Fastq                    Dataset {SCALAR}
    Data:
        (0) "@a8964d35-bbf8-4af8-afe6-2dfd7b57e7ae_Basecall_Alignment_template nanopore2_20160728_FNFAB24462_MN17024_sequencing_run_E_coli_K12_1D_R9_SpotOn_2_40525_ch2_read16_strand\n
        GTGTTGTACTTCGTTCAGTTACGTATTGCTGTTTTCGCATTTATCGTGAAACGC
   ...
   <snip>
   ...
   32+')&'')(+.+)',$'')*)+*.-1*+*+*2,.-,,/-*-.(-&#%\n"
```
From experience though, there can be strange formating quirks in the output from this command. We will demonstrate more convenient ways to extract the fastq sequences below.

Unsurprisingly h5dump dumps the entire file to STDOUT. There is quite a bit to get through, and even more in a 2D file.

```sh
h5dump ~/Data/read_data/r9_1d_ecoli_ch2_read16.fast5
```

## Browsing HDF5 files

Any HDF5 file can be opened using hdfview and browsed/edited in a Java GUI. The example files can be opened from the terminal as shown below. Alternatively, files can also be opened via the hdfview GUI.

```sh
# non base-called (2D R9)
hdfview ~/Data/read_data/r9_2d_zika_ch1_read10_pre.fast5 &

# base called (2D R9)
hdfview ~/Data/read_data/r9_2d_zika_ch1_read10.fast5 &

# base called (1D R9)
hdfview ~/Data/read_data/r9_1d_ecoli_ch2_read16.fast5 &

# base called (2D R7)
hdfview ~/Data/read_data/r7_2d_ecoli_ch56_file159.fast5 &
```

Have a browse through these files and look for the main differences (hint - check the logs and summary groups).

## Basic manipulation in poRe

poRe is a library for R available from [SourceForge](https://sourceforge.net/projects/rpore/) and published in [bioinformatics](http://bioinformatics.oxfordjournals.org/content/31/1/114).  poRe is incredibly simple to install and relies simply on R 3.0 or above and a few additional libraries. Current documentation is on [Github](https://github.com/mw55309/poRe_docs).

The poRe library is set up to read v1.1 data by default, and offers users parameters to enable reading of v1.0 data.  Let's start it up.

```sh
R
```
```R
library(poRe)
```

#### FASTQ

We'll see how to extract FASTQ from entire directories below, but here are some exemples of single file analysis. We start with the get_fastq() function. For a 1D R9 file, this returns the template fastq sequence. 

```R
f5 <- "~/Data/read_data/r9_1d_ecoli_ch2_read16.fast5"
get_fastq(f5)
```

For a 2D file, it returns a list with template, complement and 2D sequences where available.

```R
f5_2d <- "~/Data/read_data/r9_2d_zika_ch1_read10.fast5"
get_fastq(f5_2d)
```

Additional arguments to this function allow for selection of a specific fastq dataset, and the specification of alternative dataset paths (useful for older v1.0 files).

```R
? get_fastq()
```

Press q to exit the help information.

Working with lists is easy in R, and if you want the FASTQ as a string:

```R
fq <- get_fastq(f5_2d)
names(fq)
fq$template
```

From here it is quite simple to write a FASTQ extraction script. Here is a little example:

```R
# path to fast5 files
f5dir <- "Data/read_data/R9_1D_FAST5"
# get vector of all fast5 files
f5files <- dir(f5dir, pattern="\\.fast5$", full.names = TRUE)
# iterate over files
for (f5 in f5files) {
    # extract template fastq
    fq <- get_fastq(f5)
    # check fq is a list and contains template data
    if (typeof(fq) == "list" && exists("template", where=fq)) {
        # cat to "" (STDOUT but could be the name of a file
        # change the cat = "" to a filename to see what happens
        cat(fq[["template"]], file = "", sep = "\n", fill = FALSE)
    }
}
```

However, we have scripts that do this already, so there is no need to write your own!

#### FASTA data

Simply use get_fasta() instead of get_fastq(). The arguments are the same.

```R
? get_fasta()
```

#### Events data

We extract events data with the get.events() function. For 2D data files, get.events again returns a list, wth the template and complement data frames. For 1D files you just get template events:

```R
f5 <- "~/Data/read_data/r9_1d_ecoli_ch2_read16.fast5"
ev <- get.events(f5)
names(ev)
head(ev$template)

ev2 <- get.events(f5_2d)
names(ev2)
head(ev2$complement)
```

* start: time in seconds
* mean: the mean signal of the event
* stdv: the standard deviation when sampling
* length: length (duration) of the event
* model_state: the predicted kmer
* move: how many moves the model has made in that step

#### Extracting the model

For the Hidden Markov Model (HMM), R7 basecaller, the parameters of the model used to convert events to fastq sequences were stored alongside the event tables and fastq sequences. The get.models() and get.model.params() functions are provided to extract this information for R7 files. However, for R9 files basecalled using the current recurrent neural network (RNN) methods, the model is more of a black box and the parameters are not available for inspection in the FAST5 files. Since R7 is deprecated, we won't go into further detail on this here.

<!---

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

--->

#### Plotting Squiggles

Squiggle plots are staircase step plots of event mean values with event length as the step length. Plotting squiggles works directly from the events data extracted using get.events(). The parameters minseconds and maxseconds control the time interval plotted:

```R
plot.squiggle(ev$template)
plot.squiggle(ev$template,minseconds=3,maxseconds=4)

plot.squiggle(ev2$complement,maxseconds=0.5)
```

## Run QC in poRe

As we saw ealier, the poRe real-time GUI, pore_rt() outputs a metadata text file during a nanopore run and this can be used for run QC and further analysis. If pore_rt() was not run, alternative methods are available for extracting the same metadata post-run. More on this below. 

For demonstration purposes, we have uploaded 

If you ran pore_rt() during your nanopore run then you will have access to a metadata text file that can be used for run QC etc.  Otherwise we will have to create one (more below). In the meantime, we have put some example metadata files in the ~/Data/run_metadata folder. We focus on a 2D R9 run here.

```R
# load in pass data
pass <- read.table("Data/run_metadata/r9_pass_meta.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# set standard/expected column names for legacy poRe code
colnames(pass) <- c("filename","channel_num","read_num","read_start_time",
                    "status","tlen","clen","len2d",
                    "run_id","read_id","barcode","exp_start")

# load in the fail data
fail <- read.table("Data/run_metadata/r9_fail_meta.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)


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
```

<!---
# We now know the longest read
longest <- "Data/read_data/MAP006-1_2100000-2600000_fast5/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch153_file57_strand.fast5"
lfq <- get_fastq(longest, which="2D")
lfa <- get_fasta(longest, which="2D")

# write fastq and fasta out to files
cat(lfq[["2D"]], file = "longest.fastq", sep = "\n", fill = FALSE)
cat(lfa[["2D"]], file = "longest.fasta", sep = "\n", fill = FALSE)
--->

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

As is often the case, there is quite a long failed template read. These are usually nonsense:

```R
max(pass$tlen)
# [1] 37395

max(fail$tlen)
# [1] 130789
fail [fail$tlen==130789,]
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
plot.channel.summary(fail.s, report.col="ltemplate")

```

### Extracting meta-data from fast5
If you haven't run pore_rt(), then you can extract meta-data directly from the fast5 files. This takes a little time as we have to open each file and extract the attributes. For the 878 R9 1D files we have selected here, and on this VM, it takes about 15 seconds. Scaling up to a large dataset, that's around 200k files per hour.

```R
pass <- read.meta.info("Data/read_data/R9_1D_FAST5/")
                        
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

## poRe Parallel GUI

A convenient alternative method of data extraction and basic run QC is offered by a second pore GUI: pore_parellel(). The parallel GUI is designed to process batches of base-called files, rapidly extracting FASTQ and metadata using parallel code.
This GUI also provides plotting output similar to pore rt(). We can launch the parellel GUI from a terminal:

```sh
R
library(poRe)
pore_parallel()
```

There are two operating modes:

Data Extraction
* Select source, target folders
* Select output file type(s) and Dataset(s)
* Run Data Extraction

Plotting
* Choose Metadata File
* Show Plots (for selected metadata)

We'll try extracting FASTQ data, and loading up a metadata file for visualisation.

## poretools

Much of what we did above can also be achieved using poretools.  We have cheated slightly with poRe by extracting most of the metadata we need during the run and metrichor base-calling, which allows us to do things much more quickly.  Below, poretools is going to each fast5 file and ripping out the necessary information which takes a long time.  poRe also takes a long time when having to do this, though the parellel GUI is faster.

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

