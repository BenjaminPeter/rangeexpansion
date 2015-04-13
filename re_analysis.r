####################################################
# this script serves as an example on how to analyze
# range expansion data, how to infer the origin and 
# to infer the strength of the founder effect. 
# 
# there are two data sets:
# 1. the Arabidopsis thaliana data set, by Bergelson et al.,
# available on their website
# http://bergelson.uchicago.edu/regmap-data/, which I
# analyze in my biorxiv paper.
#
# 2. a dummy data set with the same input format as 
# snapp, which is (hopefully) easyier to generate and
# intended to be used for other analyses.
# 
# The script works by loading all data into
# memory so make sure that the computer you are 
# running it on has enough RAM.
####################################################



####################################################
### File format descriptions 
####################################################
# there are two required files and on optional file:
# 1. snp_file, containing the genetic data
# 2. coords_file, containing location information
# 3. outgroup_file (optional) contains outgroup info

# The data formats are the following:
#------------------------------------------------------------
# 1a) snp_file (arabidopsis):
# in arabidopsis
# each row is a SNP, each column except the first
# three an individual.
#   first row is a header with individual ids, preceeded
#       by an X to allow for numerical ids.
#   first 3 columns are:
#   1: SNP id
#   2: chromosome
#   3: snp position
#
#------------------------------------------------------------
# 1b) snp_file (snapp, default),
# each row is an individual, each column a SNP, and fields are
# comma separated, a `?`, denotes missing data, 0,1, 2 denote
# 0,1 or 2 copies of the allele. The very first column gives
# the name of the individual
#
#
#------------------------------------------------------------
# 2) coords_file 
# coords_file: each sample is a row, except header
#   (first row). 1st column is the id of the individual,
#   which should match the snp_file. Columns titled `latitude`
#   and `longitude` give sample coordinates, others are ignored
#
#
#------------------------------------------------------------
# 3) outgroup_file: outgroup data, optional, only for arabidopsis
#   as the analysis requires knowledge of the ancestral state
#   of a snp, some outgroups might be used to infer that. If
#   the data is already encoded in derived allele frequency,
#   you won't need this file.
#
#   the first 3 columns are:
#   1: SNP id (same as in the snp_file)
#   2: chromosome
#   3: snp position
####################################################

#
#
# there are three main steps to this program: 
# 1., data is loaded from individuals, then 
#   population level data is generated from individual 
#   level data.
#   SNP data for each individual, groups them in pop-
#   ulations of arbitrary size (but at least 1 diploid 
#   per population is required). 
# 2. After that, a directionalit statistc is calculated 
#   for all pairs of populations, and a file with the 
#   location for each population is generated. 
#
# These two files are then used for the actual inference
# in the find_origin function (step 3). here,
# I only included one sample set of populations in 
# the analyses, but I recommend chaing some parameters here
# to see how differerent sets of individuals can be analyzed.
####################################################


####################################################
### Paramters for the analysis (example)
####################################################
# NOTE: chance these for your data set

# the following lines contain the arguments that might be set
# alternatively, the snp_file and coords_file are loaded from
# the command line as the first two arguments, i.e. running
#
# Rscript re_analysis.r [snp_file] [coords_file]

#the name of the input file, 
snp_file <- "example_data/example_snp.snapp"

#the name of the file specifying location, with extension
coords_file <- "example_data/example_coordinates.csv" #replaced by cmdline arg 2

#whether data set needs to be loaded
load_data <- T  

#if false, only functions are loaded
run_analysis <- T  

# each list entry is a set of populations to analyze independently, i.e this
# will analyze the populations REGION_1, REGION_2 and REGION_3 individually, but will
# also jointly analyze REGION_1 and REGION_2.
regions_to_analyze <- list("REGION_1", "REGION_2", "REGION_3", 
                        c("REGION_1", "REGION_2"))

# if you instead want to infer populations from location, set this to
#pop_to_analyze <- "infer from location"

#if TRUE, heterozygostity and FST plots are generated
run_additional_analyses <- FALSE

n_points <- 4   #how many points on the map should be evaluated
                 # higher numbers increase runtime and accuracy

ploidy <- 2  #set ploidy of individuals. 1=haploid, 2 =diploid

#which columns contain outgroup individuals (snapp format)
# to be used for polarization of SNPs. If SNP are already polarized
# or no outgroups are present, set this to `NULL`:
# outgroup_columns <- NULL 
outgroup_columns <- 1:2  

# the maximum number of snp to analyze, NULL loads all SNP
nsnp <- NULL

# if you want to run the arabidopsis example instead, set this to True
#, file names will be adjusted
run_arabidopsis_example <- TRUE


# downloads arabidopsis data, requires wget
download_arabidopsis_data <- FALSE





####################################################
# minimal changes should be required from here on

if( run_arabidopsis_example ){
    ploidy <- 1
    snp_file <- "arabidopsis/athal_snps_031110_agdp"
    coords_file <- "arabidopsis/accession_coordinates.csv"
    outgroup_file <- "arabidopsis/huTsnps.txt"
    ploidy <- 2
    regions_to_analyze <- "infer from location"

    nsnp <- 1000
}


out_file_id = sprintf( "out_%s", basename( snp_file ) )
#------------------------------------------------------------
# parsing command line: 
#------------------------------------------------------------
if(length(commandArgs(T)) >0 ){
    file_id = commandArgs(T)[1]
}
if(length(commandArgs(T)) >1 ){
    location_file = commandArgs(T)[2]
}

if(length(commandArgs(T)) >2 ){
    nsnp = as.numeric(commandArgs(T)[3])
    out_file_id <- sprintf("%s_nsnp%d", out_file_id, nsnp)
} else {
    nsnp = NULL
}

psi_name <- sprintf("psi%s.txt", out_file_id)
coords_name <- sprintf("pop_coords_psi%s.txt", out_file_id)



#--------------------------------------------------
# step 0, installing libraries
#--------------------------------------------------
# libraries, if they are not already installed
#--------------------------------------------------
source("re_functions.r")
packages(sp)
packages(rworldmap)
packages(rworldxtra)
packages(geosphere)
packages(plyr)
packages(maps)
packages(fossil) 

#free some names
if(!is.data.frame(data))data <- c()
if(!is.data.frame(bbox))bbox <- c()

##################################################
# part 1: loading data
##################################################
# the goal of this step is it to load all data,
# and calculate the psi statistic for all pairs of
# individuals. Those can then be used for subseq.
# analysis (in part 2)
##################################################

#this only needs to be run once, and takes a lot of time
#after the script ran through part 1 once, set load_data
#to false to save time
if( load_data){


#--------------------------------------------------
# step 0.2, preparation: go to the following website: 
#--------------------------------------------------
# http://bergelson.uchicago.edu/regmap-data/
# Download the samples with high-quality geographic
# data, and the latitude and longitude for those accessions
# on linux systems with wget installed, the function can be
# used, users of other os will have to downlaod and
# extract the files manually. Also, the file
# accession_coordinates.xls needs to be transformed into
# a csv file, e.g. using libreoffice or excel.
# you should end up with the following files
#   - athal_snps_031110_agdp
#   - accession_coordinates.csv
#   - huTsnps.txt
#
#    
# (uncomment if needed)
#--------------------------------------------------
if( download_arabidopsis_data ){
    download_data()
}

#--------------------------------------------------
# step 1: load data into memory
# (this may take a while, and requires ~6GB of RAM)
#--------------------------------------------------

if( run_arabidopsis_example ){
    read_data_arab(snp_file, coords_file, outgroup_file,
                   nsnp)
} else {
    read_data_snapp(snp_file=snp_file,
                    coords_file=coords_file, nsnp=nsnp,
                    outgroup_columns=outgroup_columns
                   )
}

print(c("file id is", out_file_id))


#--------------------------------------------------
# step 2: define populations
#--------------------------------------------------
# now that all data is loaded, we need to assign
# individuals to populations, which will most
# likely be different for each data set. For the
# A. thaliana data set, I create populations based
# on the sample locations. For many applications, it
# might be easier to make the required data structure
# by hand. It's format is a python list, where each
# entry is a population. For example, the list
#
# [[1]]
# 1, 6, 7, 8
#
# [[2]]
# 2, 3, 4, 5
#
# would assign individuals 2-5 to population 2, 
# and individuals 1, 6, 7 and 8 to population 1.
# 
# if each individual is its own population, use the
# make_pops_auto function
# the mape_pops function assigns individuals into
# populations based on sampling locations
#--------------------------------------------------

if( run_arabidopsis_example ){
    pops <- make_pops( coords)
} else {
    pops <- make_pops_auto( coords )
}

print( "finished step 2")

#--------------------------------------------------
# step 3: make population data set
#--------------------------------------------------
# in this step, the goal is to creat two new data
# structures, representing the population level
# data and location from the individual level data
# read from the files, and the population structure
# data generated in step 2.
#
#   the following data sets are generated
#       -pop_data: a n x m matrix, where each row
#           is a SNP, and each column is a pop.
#           Entries are the absolute allele freq
#           at a SNP in each pop
#       -pop_coords
#           data frame that gives
#               latitude, longitude, country, 
#               population and sample size for each
#--------------------------------------------------
pop_data <- make_pop_data_from_pops( pops, data )
pop_coords <- make_pop_coords_from_pops( pops, coords)
pop_ss <- make_pop_ss_from_pops( pops, data, ploidy=ploidy)
pop_coords <- cbind( pop_coords, hets=get_heterozygosity(pop_data, pop_ss))

print( "finished step 3")
#--------------------------------------------------
# step 4: calculate psi statistic
#--------------------------------------------------
# the next step is to calculate the psi statistic
# for all pairs of population. This function will 
# return a n x n matrix with the pairwise psi values.
#--------------------------------------------------
all_psi <- get_all_psi(pop_data, pop_ss ,n=2)


#--------------------------------------------------
# step 5: save the data structures
#--------------------------------------------------
#after this part, we have all data condensed in a 
# n x n matrix of psi values, and a n x 5 dataframe
# of population coordinates, which can be used for
# various analyses. they are now saved so they do
# not have to be recomputed
#--------------------------------------------------
write.table(all_psi,psi_name,col.names=F,row.names=F)
write.table(pop_coords, coords_name, col.names=T,
        row.names=F)
}

##################################################
# part 3: analysis
##################################################
# in this part, we will run the various analyses, 
# to infer the origin, founder effect strength,
# and genrate various plots.
##################################################

#--------------------------------------------------
# step 2.0: reload data
#--------------------------------------------------
# first we load the data made in the first step
# from the hard drive
#--------------------------------------------------
all_psi <- read.table(psi_name, header=F)
pop_coords <- read.table(coords_name, header=T )
#--------------------------------------------------
# step 2.1: infer origin
#--------------------------------------------------
# then, we find the origin for various sets of pops
# using the find_origin function, 
# arguments:
#    pop_coords, all_psi: these are the two data 
#     files generated in the first step. 
#
#    region: selects the regions being plotted. This
#         is compared against the  pop_coords$pop
#        argument, all populations that are in
#        the region are  included in the analysis.
#
#    xlen, ylen: the number of points analysed in 
#        x and y direction, the total number
#        of points analyzed is xlen x ylen,
#        the more points are analyzed, the 
#        longer this function will run, but
#        the more accurate the output woll be.
#
#    doPlot: should a plot be produced?
#    
#    doPdf: should a pdf be generated? if yes, the
#         argument is the name of the file,
#        set to F to generate local plot.
#
#
#Returns:
#    a vector whose first two elements are the latitde and
#    longitude of the inferred origin, elements 3-7 are
#    various expressions of the founder effect strength,
#    elelemtn 8 is the r squared value and element 9 is
#    the bonferroni corrected p value.
#--------------------------------------------------





if( run_analysis ){

res.tbl <- data.frame()

for( region in regions_to_analyze ){
    print(region)
    run_region( region=region, xlen=n_points, ylen=n_points )
}




#------------------------------------------------------------
# after the script finished running, there should be a table
# with the founder effect strength and the estimated origin
# for all populations in res.tbl
#------------------------------------------------------------
write.table(res.tbl, sprintf("table_%s.txt", out_file_id), col.names=T,row.names=F)

}


if( run_additional_analyses ){


    # analysis 1: print pairwise fst
    pdf(sprintf("plots/fst_pw_%s.pdf", out_file_id))
    fst.mat <- get_all_pairwise_fst( pop_data, pop_ss)
    colnames(fst.mat) <- paste(pop_coords$pop, coords[,1])
    packages(lattice)
    o <- order(pop_coords$pop, colMeans(fst.mat))
    l <- levelplot( fst.mat[o,o], at=c(-1,seq(0,1,0.01)), 
              col.regions=c( "black", heat.colors(100)) )
    print( l )
    dev.off()
    rownames(fst.mat) <- paste(pop_coords$pop, coords[,1])
    write.table( fst.mat, sprintf("fst_table_%s.pdf", out_file_id), quote=F)

    # analysis 2: plot population Heterozygosity
    pdf(sprintf("plots/het_%s.pdf", out_file_id))
    print( xyplot(pop_coords$hets ~ pop_coords$pop) )
    dev.off()


}
