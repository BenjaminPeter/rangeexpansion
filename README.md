# Range Expansion analysis
####################################################
this script serves as an example on how to analyze
range expansion data, how to infer the origin and 
to infer the strength of the founder effect. Details
of the methods are found in Peter & Slatkin (2013), Evolution
and Peter & Slatkin, biorXiv

if you are interested in analysing a data set, please
check out the pipeline below. If you just want to 
calculate `psi` from allelic data, check out the next paragraph.


### calculation of `\psi`
 
the script implents a basic pipeline from a genetic data
set to a graphical output. *If you are only interested in
calculating the `\psi` statistic, use the `get_psi` function in
re_functions.r*

the parameters for the function are:
   Parameters:                                        
       - fi : vector of int                           
           - vector of derived allele counts for each 
               snp in population i                    
       - fj : vector of int                           
           - vector of derived allele counts for each 
               snp in population j                    
       - ni : vector of int                           
           - vector of total number of genotypes in   
               population i                           
       - nj : vector of int                           
           - vector of total number of genotypes in   
               population j                           
       -n : the number of samples to downsample to    

all vectors have length equal to the number of snp



### Pipeline

 there are two data sets:
 1. the Arabidopsis thaliana data set, by Bergelson et al.,
 available on their website
 http://bergelson.uchicago.edu/regmap-data/, which I
 analyze in my biorxiv paper.

 2. a dummy data set with the same input format as 
 snapp, which is (hopefully) easyier to generate and
 intended to be used for other analyses.
 
 The script works by loading all data into
 memory so make sure that the computer you are 
 running it on has enough RAM.

I am currently moving comments and documentation from
re_analysis.r to this readme file, so please also check
the initial comments there.


## File format descriptions 
####################################################
 there are two required files and on optional file:
 1. snp_file, containing the genetic data
 2. coords_file, containing location information
 3. outgroup_file (optional) contains outgroup info

 The data formats are the following:

 1a) snp_file (arabidopsis):
 in arabidopsis
 each row is a SNP, each column except the first
 three an individual.
   first row is a header with individual ids, preceeded
       by an X to allow for numerical ids.
   first 3 columns are:
   1: SNP id
   2: chromosome
   3: snp position

 1b) snp_file (snapp, default),
 each row is an individual, each column a SNP, and fields are
 comma separated, a `?`, denotes missing data, 0,1, 2 denote
 0,1 or 2 copies of the allele. The very first column gives
 the name of the individual


 2) coords_file 
 coords_file: each sample is a row, except header
   (first row). 1st column is the id of the individual,
   which should match the snp_file. Columns titled `latitude`
   and `longitude` give sample coordinates, others are ignored


 3) outgroup_file: outgroup data, optional, only for arabidopsis
   as the analysis requires knowledge of the ancestral state
   of a snp, some outgroups might be used to infer that. If
   the data is already encoded in derived allele frequency,
   you won't need this file.

   the first 3 columns are:
   1: SNP id (same as in the snp_file)
   2: chromosome
   3: snp position


Overview of steps
-----------------

 there are three main steps to this program: 
 1. data is loaded from individuals, then 
   population level data is generated from individual 
   level data.
   SNP data for each individual, groups them in pop-
   ulations of arbitrary size (but at least 1 diploid 
   per population is required). 
 2. After that, a directionalit statistc is calculated 
   for all pairs of populations, and a file with the 
   location for each population is generated. 

 3. These two files are then used for the actual inference
 in the find_origin function (step 3). here,
 I only included one sample set of populations in 
 the analyses, but I recommend chaing some parameters here
 to see how differerent sets of individuals can be analyzed.



## Paramters for the analysis (example)
*NOTE: chance these for your data set in re_analysis.r*

the following lines contain the arguments that might be set
alternatively, the snp_file and coords_file are loaded from
the command line as the first two arguments, i.e. running

    Rscript re_analysis.r [snp_file] [coords_file]

the name of the input file, 

    snp_file <- "example_data/example_snp.snapp"

the name of the file specifying location, with extension

    coords_file <- "example_data/example_coordinates.csv" #replaced by cmdline arg 2

whether data set needs to be loaded, if they are already in the R environment, setting
this to `FALSE` will save a lot of time

    load_data <- True

if `FALSE` functions are loaded, but no code is executed

    run_analysis <- True

##### Regions
Regions are sets of populations that can be analyzed independently, i.e. they correspond
to clusters that a priory are thought to have a different origin.
Each list entry corresponds to an analysis, i. e. the following command
will analyze the populations `REGION_1`, `REGION_2` and `REGION_3` individually, but will
also jointly analyze `REGION_1` and `REGION_2`.

    regions_to_analyze <- list("REGION_1", "REGION_2", "REGION_3", 
                        c("REGION_1", "REGION_2"))

if `TRUE`, heterozygostity and FST plots are generated (note that I use BEDASSLE
to calculate pairwise FST, and that BEDASSLE is, as of Dec 3 2014, not yet available
for 3.1.1, but this should work for older versions of R)

    run_additional_analyses <- FALSE

    n_points <- 4   #how many points on the map should be evaluated
                 # higher numbers increase runtime and accuracy

    ploidy <- 2  #set ploidy of individuals. 1=haploid, 2 =diploid

which columns contain outgroup individuals (snapp format)
to be used for polarization of SNPs. If SNP are already polarized
or no outgroups are present, set this to `NULL`:
    # outgroup_columns <- NULL 

    outgroup_columns <- 1:2  

the maximum number of snp to analyze, NULL loads all SNP

    nsnp <- NULL

if you want to run the arabidopsis example instead, set this to True
, file names will be adjusted

    run_arabidopsis_example <- TRUE


downloads arabidopsis data, requires wget. As the arabidopsis data set is around
400MB, I did not include it, run this *once* to download the data

    download_arabidopsis_data <- FALSE


