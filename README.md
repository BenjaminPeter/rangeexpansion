# Range Expansion analysis
####################################################
this package serves as an example on how to analyze
range expansion data, how to infer the origin and 
to infer the strength of the founder effect. Details
of the methods are found in Peter & Slatkin (2013), Evolution
and Peter & Slatkin (2015), Evolution

The major uses of the package are to either test
isolation-by-distance on some sets of populations, or
to infer the origin of an expansion, and possibly the
strength of a founder effect.

The current version uses the snpStats package from
Bioconductor to handle SNP data efficiently. See
http://www.bioconductor.org/packages/release/bioc/html/snpStats.html
on how to install it.


### File Format Descriptions & Modifiers
- the required files are genetic data in the plink
.bed/.bim/.fam file format. 
- In addition, we require a file containing information on the geographical
sampling coordinates, and possibly other information

 The data formats are the following:
###### snp_file (bed format)
see http://pngu.mgh.harvard.edu/~purcell/plink2/formats.html
for reference. This file is perhaps easiest generated using
plink 1.9. For example, from vcf data, use
plink --vcf my.vcf --set-missing-var-ids "@_#" --make-bed --out new --allow-extra-chr
to generate a bed file.
###### coords_file 
a tab-delimited file with a header (first row)
(first row). Headers are (by default specified as follows):
- individual (or 1st column): is the id of the individual
- longitude: the longitude or x-coordinate of sample location
- latitude: the latitude or y-coordinate of the sample
- outgroup (optional): Whether the individual should be treated as an outgroup
to polarize SNP
- region (optional): In the case of multiple origins, this allows for an
factor for which individuals are included in a given run
- population (optional): If custom population assignments are chosen, the 
value here assigns individuals to a population
- country (optional): The country of origin.


### Pipeline Overview:
there are three main steps to this program: 
1. data is loaded from individuals, then population level
 data is generated from individual level data
2. After that, a directionalit statistc is calculated 
 for all pairs of populations, and a file with the 
 location for each population is generated. 
3. These population matrices are used for the actual inference

### Modelling decisions and options

###### Outgroups
Outgroup Individuals are individuals assumed to be ancestral to the population of interest.
Their alleles can be used to determine the most likely ancestral state of an allele.
If more than one outgroup is specified, there might be disagreements between the state
of the outgroup. In this case, there are two options:
    outgroup.disagreement = 'remove' #removes SNP whose ancestral state is ambiguous
    outgroup.disagreement = 'majority' # remove ties, keep SNP with majority for an allele
######    Regions
Regions are sets of populations that can be analyzed independently, i.e. they correspond
to clusters that a priory are thought to have a different origin.
Each list entry corresponds to an analysis, i. e. the following command
will analyze the populations `REGION_1`, `REGION_2` and `REGION_3` individually, but will
also jointly analyze `REGION_1` and `REGION_2`.

    regions_to_analyze <- list("REGION_1", "REGION_2", "REGION_3", 
                        c("REGION_1", "REGION_2"))

######    Bounding box

### calculation of ψ 
 
the script implents a basic pipeline from a genetic data
set to a graphical output. If you are only interested in
calculating the ψ statistic, use the `get_psi` function

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


Overview of steps
-----------------

 there are three main steps to this program: 
 1. data is loaded from individuals, then 
   population level data is generated from individual 
   level data.
   SNP data for each individual, groups them in pop-
   ulations of arbitrary size (but at least 1 diploid 
   per population is required). 
 2. After that, a directionality statistic is calculated 
   for all pairs of populations, and a file with the 
   location for each population is generated. 

 3. These two files are then used for the actual inference
 in the find_origin function (step 3). here,
 I only included one sample set of populations in 
 the analyses, but I recommend chaing some parameters here
 to see how differerent sets of individuals can be analyzed.



## Example analysis
#### Specify File Names
First, we need to specify the files
    snp.file <- "example_data/example_snp.bed"
    coord.file <- "example_data/example_coordinates.csv" 

Note that the bim and fam files will be automatically found by the `snpStats::read.plink` function.

#### Specify Options

    ploidy <- 2 #assume our organism is diploid
    nsnp <- NULL # assume we want to analyze all SNPs in the data set

    # we want two runs, one with all individuals in REGION_1, and one
    # with individuals in either REGION_2 or REGION_3
    regions <- list("REGION_1",
                    c("REGION_2", "REGION_3"))

    n_points <- 20 # at how many points should the function by analyzed?


#### Load Individual Level Data
    snp.data <- load.plink.file(snp.file)
    coord.data <- load.coord.file(coord.file, sep=',')
    raw.data <- check.missing(snp.data, coord.data)
    raw.data <- set.outgroups(raw.data, ploidy)
    pops <- make.pops(raw.data)
#### Generate Population Level Data
    pop_data <- make_pop_data_from_pops( pops, raw.data )
    pop_coords <- make_pop_coords_from_pops( pops, coords)
    pop_ss <- make_pop_ss_from_pops( pops, data, ploidy=ploidy)
    pop_coords <- cbind( pop_coords, hets=get_heterozygosity(pop_data, pop_ss))
#### Run Analyses
    all_psi <- get_all_psi(pop_data, pop_ss ,n=2)
    # write.table(all_psi,psi_name,col.names=F,row.names=F)
    # write.table(pop_coords, coords_name, col.names=T,
    #        row.names=F)

    run_region( region=region, xlen=n_points, ylen=n_points )
    write.table(res.tbl, sprintf("table_%s.txt", out_file_id), col.names=T,row.names=F)



