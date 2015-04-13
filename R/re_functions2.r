require(snpStats)
require(geosphere)

#' @import snpStats
#' @import geosphere
#' @import rworldmap
NULL

#' Reads VCF File
#'
read.data.vcf <- function(){
    stop("Function Not Yet Implemented")
    require(GGtools)
}

#' Reads data from SNAPP file
#' 
#' Deprecated, use plink input format instead
#' @example examples/example_1.r
load.snapp.file <- function(snp.file="../example_data/example_snp.snapp",
                            n.snp=-1){
    data <- read.table(snp.file, sep =",", header=F, strings=F,
                       na.strings="?", nrow=n.snp, row.names=1)
    raw.data <- list()
    raw.data$genotypes <- as(as.matrix(data), "SnpMatrix")
    raw.data$fam <- rownames(data)
    return(raw.data)
}

#' Wrapper around snpStats::read.plink to read data in bed/bim/fam format
#'
#' Data can be converted to bed/bim/fam format e.g.
#' using plink 1.9 or bcftools. A useful command line
#' to use is
#' \code{plink --vcf my.vcf --set-missing-var-ids "@@.#" --make-bed --out new --allow-extra-chr}
#' @export
#' @inheritParams snpStats::read.plink
#' load.plink.file()
#' @examples
#' \dontrun{
#' f = system.file("examples/example.data.bed", package="rangeExpansion")
#' load.plink.file(f, sep=",")
#' }
load.plink.file <- function(..., n.snp=NULL){
    if(!is.null(n.snp)){
        f <- read.plink(..., select.snps=1:n.snp)
    } else {
        f <- read.plink(...) 
    }
    return(f)
}



#' loads a data set from a plink and coords file
#'
#'
#' This is the main function to load a data set. It requires a 
#' plink file in bed format and a coords file
#' @param plink.file the name of a plink bed file, passed to read.plink
#' @param coords.file a file containing coordinates and other info 
#'      about the sample
#' @param n.snp the maximum number of snp to be read. if NULL, all
#'      snps are read
#' @param ploidy the ploidy of the organism, 1 for haploids, 2 for diploid
#' @param ... further arguments passed to load.coord.file
#' @return an object of type origin.data with entries `genotypes`,
#'      containing the genetic data and entry `coord` containing location
#'      data
#' @export
load.data <- function(plink.file, coords.file, n.snp=NULL, 
                      ploidy=2, ...){
    f <- load.plink.file(plink.file, n.snp=n.snp)
    coords <- load.coord.file(coords.file, ...)
    raw.data <- check.missing(f, coords)
    raw.data <- set.outgroups(raw.data)
    class(raw.data) <- 'origin.data'
    return(raw.data)
}

#' loads a data set from a snapp and coords file
#'
#'
#' This is a function that reads the old data format. Presented for backwards
#' compatibility
#' @param snapp.file a file in snapp data format
#' @param coords.file a file containing coordinates and other info 
#'      about the sample
#' @param n.snp the maximum number of snp to be read. if -1, all
#'      snps are read
#' @param ploidy the ploidy of the organism, 1 for haploids, 2 for diploid
#' @param ... further arguments passed to load.coord.file
#' @return an object of type origin.data with entries `genotypes`,
#'      containing the genetic data and entry `coord` containing location
#'      data
#' @example examples/example_1.r
#' @export
load.data.snapp <- function(snapp.file, coords.file, n.snp=-1, 
                      ploidy=2, ...){
    f <- load.snapp.file(snapp.file, n.snp=n.snp)
    coords <- load.coord.file(coords.file, ...)
    raw.data <- check.missing(f, coords)
    raw.data <- set.outgroups(raw.data, ploidy)
    class(raw.data) <- 'origin.data'
    return(raw.data)
}


#' Reads MS File
#'
read.data.ms <- function(){
    stop("Function Not Yet Implemented")
}

#' Loads a coordinate data set
#'
#' Loads a file that specifies location data and other sample specific information
#' A
#'
#' @param file: the name of the file which the data are to be read from.
#' @param ... Additional arguments passed to read.table
#' @return A data.frame with columns sample, longitude, latitude, region, 
#'outgroup and country 
#' @example examples/example_1.r
#' @export
load.coord.file <- function(file, ...){
    coords <- read.table(file, header=T, strings=F, ...)
    names(coords)[1] <- 'id'
    return(coords)
}

#' Cleans genotype and coordinate files for consistency
#'
#' 
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param snp.data output from load.plink.file
#' @param coord.data output from load.coord.file
#' @return An object with entries 'genotypes' and 'coords' containing the
#'    input data sets
#'
#' @example examples/example_1.r
#' @export
check.missing <- function(snp.data, coord.data){
    ids <- rownames(snp.data$genotypes)
    to.remove <- !unlist(lapply(ids,
            function(x)x%in%coord.data$id))
    res.data <- list()
    res.data$coords <- coord.data
    res.data$genotypes <- snp.data$genotypes[!to.remove,]
    ids <- ids[!to.remove]

    data.ordering <- res.data$coords$id
   
    res.data$genotypes <- res.data$genotypes[data.ordering,]
    
    return(res.data)
}

#' Sets columns of outgroup individuals
#'
#' As psi requires the derived state of an allele to be known, this function
#' provides an easy way to polarize SNP with outgroup individuals. From the
#' outgroup columns in the coord file, outgroup individuals are removed
#' and SNPs swithced if necessary
#'
#' @param data output from check missing
#' @param ploidy the number of copies per individual
#' @return An object with entries 'genotypes' and 'coords' containing the
#'    input data sets
#'
#' @example examples/example_1.r
#' @export
set.outgroups <- function(data, ploidy=2
                          ){
    if(is.null(data$coords$outgroup)){ 
        return(data)
    }

    outgroup.columns <- which(as.logical(data$coords$outgroup))
    n.outgroups <- length(outgroup.columns)
    if( n.outgroups > 0 ){
        data$outgroup <- data$genotypes[outgroup.columns,]
        data$outgroup <- as(data$outgroup, 'numeric')
        data$genotypes <- data$genotypes[-1* outgroup.columns,]

        n.alleles <- colSums(!is.na(data$outgroup))
	derived <- colSums(data$outgroup==ploidy, na.rm=T) == n.alleles
	anc <- colSums(data$outgroup==0, na.rm=T)== n.alleles
	no.outgrp <- is.na(derived+anc) | !(anc | derived)  

        data$genotypes <- switch.alleles(data$genotypes,derived)

	data$genotypes <- data$genotypes[,!no.outgrp ]
        data$coords <- data$coords[-outgroup.columns,]
    }

    return(data)
}


#' Puts individuals into populations for further analysis
#'
#' There are two main modes of grouping individuals into populations: 
#' The first one is based on the location, all individuals with the 
#' same sample location are assigned to the same population.
#' the other one is based on a column `population` in 
#'
#'
#' @param raw.data Output from set.outgroups or check.data
#' @param mode One of 'coord' or 'custom'. If 'custom', populations are
#'   generated from the population column in the coords file. If 'coord',
#'   Individuals are grouped according to their coordinates
#' @return A list with each entry being a population
make.pops <- function(raw.data, mode='coord'){
    if(mode == 'coord'){
        s <- c('latitude', 'longitude')
        rowSums <- function(x)200*x[,1]+x[,2]
        raw.data$coords$pop <-  match(rowSums(raw.data$coords[,s]), 
                                      rowSums(raw.data$coords[,s]))
    } 
    return(raw.data)
}



#' Sets up the popultion structure data sets
#'
#' Groups up raw data into a data structure that 
#' 
#' @param raw.data Output from set.outgroups or check.data
#' @param ploidy the ploidy of organisms
#' @return A list with entries
#' - data : an n x p matrix with derived allele counts
#' - coords : an p x 3 matrix of coords
#' - n : number of populations
#' @example examples/example_1.r
#' @export
make.pop <- function(raw.data, ploidy=2){
    if(is.null(raw.data$coord$pop)){
        raw.data<- make.pops(raw.data)
    }
    pop <- list()
    pop$data <- make.pop.data(raw.data)[,-1]
    pop$ss <- make.pop.ss(raw.data, ploidy)[,-1]
    pop$coords <- make.pop.coords(raw.data)
    pop$coords$hets <- rowMeans(hets(pop$data/pop$ss), na.rm=T)
    pop$n <- nrow(pop$coords)
    class(pop) <- 'population'
    return(pop)
}
make.pop.data <- function(raw.data){
    data <- as(raw.data$genotypes, 'numeric')
    aggregate(data, by=list(raw.data$coords$pop), FUN=sum, na.rm=T)
}
make.pop.ss <- function(raw.data, ploidy=2){
    data <- as(raw.data$genotypes, 'numeric')
    aggregate(data, by=list(raw.data$coords$pop), 
              FUN=function(x)ploidy*sum(!is.na(x)))
}
make.pop.coords <- function(raw.data){
    lat <- aggregate(raw.data$coords$latitude, 
                     by=list(raw.data$coords$pop), 
                     FUN=mean)
    names(lat) <- c("pop", "latitude")
    long <- aggregate(raw.data$coords$longitude, 
                      by=list(raw.data$coords$pop), 
                      FUN=mean)
    names(long) <- c("pop", "longitude")
    pop.coords <- merge(lat,long)
    region.id <- match(pop.coords$pop, raw.data$coords$pop)
    pop.coords$region <- raw.data$coords$region[region.id]
    return(pop.coords)
}


#' Heterozygosity function
#' @param x the allele frequency to calculate heterozygosity for
#' @return Heterozygosity = x(1-x)
#' @examples
#' freqs <- 0:10/10
#' h <- hets(freqs)
hets <- function(x)x * (1-x)


#' calculates the psi matrix for a pop object
#' 
#' @param pop population data object from make.pop
#' @param n the sample size which we downsample to
#' @param resampling mode of resampling. Currently, only
#'    hyper is supported
#' @return A n x n matrix of psi values
#' @example examples/example_1.r
#' @export
get.all.psi <- function(pop, n=2,
                    subset=NULL,
                    resampling="hyper"){
#this function calculates psi for columns i,j, both resampled down to
# n samples
    
    if(is.null(subset))
        subset <- 1:pop$n
    if( is.logical( subset) )
        subset <- which(subset)
    n.pops <- length(subset)
    mat = matrix( 0, nrow=n.pops, ncol=n.pops )
    for(i in 1:(n.pops-1)){
        for(j in (i+1):n.pops){
	    ii <- subset[i]
	    jj <- subset[j]
            ni <- unlist(pop$ss[ii,])
            nj <- unlist(pop$ss[jj,])
            fi <- unlist(pop$data[ii,])
            fj <- unlist(pop$data[jj,])
            mat[j,i] <- get.psi( ni, nj, fi, fj, 
                                resampling=resampling, n=n )
            mat[i,j] <- -mat[j,i]
	    #print( c(ii, jj))
        }
    }
    
    return(mat)
}

#' the psi statistic calculation
#' the actual calculation of the psi statistic for
#' a single pair of populations. The function requires
#' 4 vectors, all of length equal to the number of snps
#' to be analyzed.
#
#'
#' @param ni the number of sampled haplotypes in population i
#' @param nj the number of sampled haplotypes in population j
#' @param fi the number of derived alleles in population i
#' @param fj the number of derived alleles in population j
#' @param n the number of samples to downsample to
#' @example examples/example_1.r
#' @return psi a matrix of pairwise psi values
get.psi <- function (ni, nj, fi, fj,
                     n=2, resampling="hyper"){

        fn <- cbind( fi, ni, fj, nj)

        tbl <- table(as.data.frame(fn))
        tbl <- as.data.frame( tbl )
        tbl <- tbl[tbl$Freq > 0,]

        tbl$fi <- as.integer(as.character( tbl$fi ))
        tbl$fj <- as.integer(as.character( tbl$fj ))
        tbl$ni <- as.integer(as.character( tbl$ni ))
        tbl$nj <- as.integer(as.character( tbl$nj ))
        
        to.exclude <- tbl$fi == 0 | tbl$fj == 0 | 
            tbl$ni < n | tbl$nj < n

        tbl <- tbl[! to.exclude, ]

        if(nrow(tbl)==0){ return(NaN)}


        poly.mat <- matrix(0,nrow=n+1, ncol=n+1)
        poly.mat[2:(n+1),2:(n+1)] <- 1
        poly.mat[n+1,n+1] <- 0

        #psi.mat is the contribution to psi for each entry
        psi.mat <- outer(0:n,0:n,FUN=function(x,y)(y-x))
        psi.mat[1,] <- 0
        psi.mat[,1] <- 0

        f.contribution <- function(row, b=2){
            a <- 0:b
            f1 <- row[1]
            n1 <- row[2]
            f2 <- row[3]
            n2 <- row[4]
            cnt <- row[5]
            q1 <- choose(b, a) * choose(n1-b, f1-a)/choose(n1,f1)
            q2 <- choose(b, a) * choose(n2-b, f2-a)/choose(n2,f2)
            return( cnt * outer(q2, q1) )
        }


        resampled.mat <- matrix(rowSums(apply(tbl, 1,
                                                  f.contribution)),nrow=n+1)
        
        return( sum(resampled.mat * psi.mat) / sum(resampled.mat * poly.mat) )
}

#' runs TDOA analysis for multiple regions
#' 
#' runs the analysis for all individuals in a region. See
#' run.single.region for further details on argumetns
#' @param a list of vectors, where each vector contains elements
#'   that are compared to pop@@coords@@region
#'   to determine whether they should be included
#'   if null, all are run
#' @param pop population data object from make.pop
#' @param psi a matrix of psi values
#' @param output.file.name name of output file
#' @param xlen number of points to evaluate function in x direction
#' @param ylen number of points to evaluate function in y direction
#' @param ... further arguments passed to plot
#' @example examples/example_1.r
#' @export
run.regions <- function(region, ...){
    res <- list()
    tbl <- lapply(region, run.single.region, ...)
    res$tbl <- tbl
    res$regions <- region
    class(res) <- 'origin.result.list'
    return(res)
}

#' runs analysis for a single region
#' 
#' runs the analysis for all individuals in a region
#' @param pop population data object from make.pop
#' @param psi a matrix of psi values
#' @param region a vector of entries that are compared to pop@@coords@@region
#'   to determine whether they should be included
#' @param loc.file.id output stuff
#' @param xlen number of points to evaluate function in x direction
#' @param ylen number of points to evaluate function in y direction
#' @param bbox type of the bounding box to limit region
#' @param ... further arguments passed to plot
#' @export
run.single.region <- function(region, pop, psi,
                              output.file.name='default', 
                              xlen=100, ylen=100,bbox=NULL,
                              ...){
    reg.str <- paste(region, collapse="+")
    print(c("fid", output.file.name))
    dir.create('plots', showWarnings = F)
    pdf.name<- sprintf("plots/orig.%s.%s.%d.pdf", reg.str, output.file.name,
                       xlen)
    print(pdf.name)
    if(is.null(region)){
        subsample <- T
    } else {
        subsample <- pop$coords$region %in% region
    }
    subset.coords <- pop$coords[subsample,]
    subset.psi <- psi[subsample, subsample]
    origin1 <- find.origin(subset.coords, subset.psi,
                  region=region,
                   xlen=xlen,ylen=ylen,
                   doPlot=T, doPdf=pdf.name)
    return(origin1)
}


#' finds the origin for a set of Populations
find.origin <- function(coords, psi, 
                        xlen=50,  ylen=50, doPlot=F, doPdf=F, 
             ...){
    tdoa.data <- prep.tdoa.data(coords, psi)
    bbox <- get.sample.bbox(coords)

    res <- single.origin(tdoa.data,
                           bbox = bbox,
                           coords,
                           xlen=xlen,
                           ylen=ylen, ...)


    return( res )
}

#' preps data for tdoa analysis
#'
#' makes a data structure of format xi, yi, xj, yj, psi
#' 
#' @param pop.coords coordination data with latitude and longitude
#' @param all.psi matrix with psi values
#' @param region region or set of regions to run analysis on
#' @param countries countries to restrict analysis to
#' @param xlen number of x points
#' @param ylen number of y points
#'
prep.tdoa.data <- function(coords, psi){
    locs <- coords[,c('longitude', 'latitude')]
    n.locs <- nrow(locs)

    tdoa.data <- c()
    for(i in 1:n.locs){
        for(j in (i+1):n.locs){
            if( i>=j | j >n.locs) break
            tdoa.data <- rbind( tdoa.data, c(locs[i,], locs[j,], psi[i,j]))
        }
    }

    tdoa.data <- matrix(unlist(tdoa.data), ncol=5)

    return(tdoa.data)
}


#' Runs the entire analysis.
#'
#' See individual functions for details
#' @param snp.file the bed file with genetic data
#' @param coord.file file with coordinate info
#' @param pop.mode How populations are created. Default is that individuals
#' with same location are merged
#' @param regions Independent regions for which analysis is supposed to be ran
#' @param n.points How many points should be generated for plotting
#' @param n.snp Maximum number of snp to be used
run.analysis <- function(snp.file, coord.file, sep=',',
                         ploidy=2,
                         pop.mode='coord',
                         regions=NULL,
                         n.points=20,
                         n.snp=NULL
                         ){
    
    raw.data <- load.data(snp.file, coord.file, n.snp=n.snp, ploidy=ploidy,
                          sep=sep)
    pop <- make.pop(raw.data, ploidy)
    psi <- get.all.psi(pop ,n=2)

    res <- run.regions(region=region, pop=pop, psi=psi, xlen=10,ylen=20)

}

#' coords2country finds country for a given coordinate
#' @param points coordinates to get coordinate from
coords2country = function(points){  
    countriesSP <- getMap(resolution='high')
    pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
    
     
    # use 'over' to get indices of the Polygons object containing each point 
    indices = over(pointsSP, countriesSP)
     
    indices$ADMIN  
}


#' Finds origin for a region
#' @param tdoa.data Object of type tdoa.data
#' @param bbox bounding box describing the location where inference should
#'   be done in
#' @param pop.coords coordinations
#' @param pct threshold parameter in model.1d
#' @param xlen, ylen parameters describing the number of points to use
#' @param exclude.ocoean boolean, whether points not on land should be
#'    excluded
#' @param exclude.land boolean, whether points on land should be
#'    excluded
single.origin <- function(tdoa.data, bbox,  pop.coords,
                 pct=0.01,
            xlen=100, ylen=100, 
            exclude.ocean=T,
            exclude.land=F,
            ...){
    
    
    #define locs for estimate
    s1<-seq(bbox[1,1],bbox[1,2], length.out=xlen)
    s2<-seq(bbox[2,1],bbox[2,2], length.out=ylen)
    coords <- expand.grid(s1,s2)
    ij <- expand.grid(1:length(s1), 1:length(s2))

    if(exclude.ocean){
        cc <- coords2country(coords)
        to.keep <- !is.na(cc)
    } else {
        to.keep <- rep(T, nrow(coords))
    }
    if(exclude.land){
        cc <- coords2country(coords)
        to.keep <- is.na(cc)
    }
    # init output
    d0 <- matrix(NA, ncol=ylen, nrow=xlen)
    rsq <- matrix(NA, ncol=ylen, nrow=xlen)
    mdlq <- list()
    for(i in 1:xlen){
        mdlq[[i]] <- list()
    }
    
    
    for(r in 1:nrow(coords)){
        i <- ij[r,1]
        j <- ij[r,2]
        x <- coords[r,1]
        y <- coords[r,2]

        if(to.keep[r]){
            mdl <- model.1d( xy=c(x,y), data=tdoa.data, pct=pct)
            d0[i, j] <- mdl$f.e
            rsq[i, j] <- mdl$rsq
            mdlq[[i]][[j]] <- mdl
            print(i)
        }
    }
    res <- list( d0=d0, rsq=rsq, mdlq=mdlq, bbox=bbox, xlen=xlen, 
                ylen=ylen, coords=pop.coords)
    class(res) <- 'origin.results'
    return(res)
}



#' gets a bounding box around the populations in samples
get.sample.bbox <- function(samples){
    s <- c('longitude', 'latitude')
        mins <- apply(samples[,s],2,min)
        maxs <- apply(samples[,s],2,max)
        return(cbind(mins, maxs))
}

#' calculates distance from a given point
#' @param xy coordinate of point to evalute function at
#' @param data a 5 col data frame with columns xi, yi, xj, yj, psi, with 
#'   coords from the two sample location and their psi statistic.
#'   best generated unsg prep.tdoa.data
#' @param f.dist the distance function to use. the default 'haversine'
#'   uses the haversine distance. Alternatively, 'euclidean' uses 
#'   Euclidean distance
#' @return an object of type lm describing fit
model.1d <- function(xy, data, pct=0.01, f.dist="haversine"){
    if (f.dist=="euclidean"){
    f.dist <- function(i, j){
        sqrt((ix -jx)^2 + (iy-jy)^2 )
    }
    }else{ if(f.dist=="haversine"){
        f.dist <- distHaversine
    }}

    y = xy[2] 
    x = xy[1]
    ixy = data[,1:2]
    jxy = data[,3:4]
    psi = data[,5]


    d = f.dist(ixy, c(x,y)) - f.dist(jxy, c(x,y))
    l = lm( psi ~ d ) 
    l$f.e = .5 * pct/ l$coefficients[2]
    l$rsq = summary(l)$r.squared

    return (l)
}

#' given output from single.origin, finds the origin
#' @param res an object of type origin.results to be summarized
#' @return a data frame with some summary statistics and the
#'    most likely origin
#' @export
summary.origin.results <- function(res){
    bbox <- res$bbox
    orig <- which(t(res[[2]])==max((res[[1]]>0)*res[[2]],
                                   na.rm=T), arr.ind=T)        
    s1<-seq(bbox[1,1],bbox[1,2],length.out=res$xlen) 
    s2<-seq(bbox[2,1],bbox[2,2],length.out=res$ylen) 

    coords_orig <- c(s1[orig[1,2]], s2[orig[1,1]])
    max.reg <- res[[3]][[orig[1,2]]][[orig[1,1]]]

    smry <- summary(max.reg)
    q <- smry$coefficients[2,1]
    r1 = 1/(1+2*1000*q)
    r2 = 1/(1+2*10000*q)
    r3 = 1/(1+2*100000*q)
    r <- 0.99
    d1pc <- (1-r)/(2*q*r)

    res.data <- ( c( coords_orig, q*1000, r1, r2, r3, d1pc/1000, smry$adj.r.squared,
      smry$coefficients[2,4] * 10^4) )
    res.data <- as.data.frame(t(res.data))
    names(res.data) <- c("longitude", "latitude", "q", 
                         "r1","r10","r100",
                         "d1","rsq","pval")
    return(res.data)
}


#' summarizes origin.result.list object
#' given output from origin inference, returns a table with statistics
#' @param res an object of type origin.results to be summarized
#' @return a data frame with some summary statistics and the
#'    most likely origin
#' @export
summary.origin.result.list <- function(res){
    a <- sapply(res$tbl, summary)
    a <- data.frame(a)
    na <- sapply(res$regions, paste, collapse="+")
    na[na==""] <- "ALL"
    names(a) <- na
    return(a)
}

#' plots a set of origin-inference results
#'
#' simply loops over all different regions and calls the plot function
#' for all of them
#' @param res an object of type origin.result.list
#' @param ... further objects passed. to plot.origin.result
#' @export
plot.origin.result.list <- function(res, ...){
    lab <- sapply(res$regions, paste, collapse="+")
    lab[lab==""] <- "ALL"
    n.plots <- length(lab)
    if(n.plots > 1 && n.plots <6){
        par(mfrow=c(n.plots %/% 2, 2))
    } else if(n.plots >1){
        par(mfrow=c(3,2))
    }
    
    for(i in 1:n.plots){
        plot(res$tbl[[i]], main=lab[i], ...)

    }
}



#' plots the output of find.origin as a heatmap using .filled.contour
#' @param x an object of type origin.results, as obtained by 
#' @param n.levels the number of color levels
#' @param color.function a function that takes an integer argument and
#'    returns that many colors
#' @param color.negative a single color to be used for negative values
#' @param add.map boolean whether a map should be added
#' @param add.likely.origin boolean, whether origin should be marked with an
#' @param asp aspect ratio, set to 1 to keep aspect ratio with plot
#'   X
#' @export
plot.origin.results <- function(x, n.levels=100, color.function=heat.colors,
                                color.negative='grey',
                                add.map=T,
                                add.samples=T,
                                add.sample.het=T,
                                add.likely.origin=T,
                                asp=1,
                                ...){
    plot.default(NA, xlim=x$bbox[1,], ylim=x$bbox[2,], 
         xlab="", ylab="",
         xaxt="n", yaxt="n",
         xaxs='i', yaxs='i',
         asp=asp, ...)

    s1<-seq(x$bbox[1,1],x$bbox[1,2],length.out=x$xlen)
    s2<-seq(x$bbox[2,1],x$bbox[2,2],length.out=x$ylen)

    rel <- (x[[1]]>0) * (x[[2]]-min(x[[2]],na.rm=T)) /
                   (max(x[[2]],na.rm=T)-min(x[[2]],na.rm=T))+0.001

    levels <- c(0,quantile(rel[rel>0.001], 0:n.levels/n.levels, na.rm=T) + 
                1e-6 * 0:n.levels/n.levels)
    cols <- c(color.negative, color.function(n.levels-1))


    .filled.contour(s1, s2, rel, levels, cols)


    # rect(x$bbox[1,1], x$bbox[2,1], x$bbox[1,2], x$bbox[2,2], border='black',
    #      lwd=2, col=NULL)

    if(add.likely.origin){
        points(summary(x)[,1:2], col='black', pch='x', cex=2)
    }
    if(add.map){
        require(rworldmap)
        m <- getMap("high")
        plot(m, add=T, lwd=1.3)
    }

    if(add.sample.het){
        samples <- x$coords
        hets <- (samples$hets - min(samples$hets) )/(
            max(samples$hets) - min(samples$hets))
        points( samples$longitude, samples$latitude,
            pch=16, cex=3, col=grey(hets) )
        points( samples$longitude, samples$latitude,
            pch=1, cex=3, col="black",lwd=2 )
    }
    else if(add.samples){
        points( samples$longitude, samples$latitude,
            pch=16, cex=1, col="black",lwd=1 )
    }
}



############################################################
# unfinished
############################################################
permute.data <- function(data, pops, n.permutations=1000){
    permutation.matrix <- sapply(1:n.permutations, 
                                 function(x)sample(1:ncol(data)))
    permutations <- array(0, dim=c(dim(all.psi), n.permutations))

    for( i in 1:n.permutations){
        per.data <- make.pop.data.from.pops(pops, data[permutation.matrix[,i]])
        per.ss <- make.pop.ss.from.pops(pops, data[permutation.matrix[,i]])
        permuted.psi <- get.all.psi(per.data, per.ss)
        permutations[,,i] = permuted.psi
        print(i)
    }

    permutations
}


#' tests for isolation by distance by doing permutations on data
#' @export
test.ibd <- function(data, all.psi, pops, n.permutations=1000){
    permuted.data <- permute.data(data, pops, n.permutations)
    n.samples <- dim(all.psi)[1]
    res <- c()
    for( i in 1:(n.samples-1)){
        for (j in (i+1):n.samples){
            l <- sum(all.psi[i,j] < permuted.data[i,j,])
            u <- sum(all.psi[i,j] > permuted.data[i,j,])
            res <- rbind(res, c(i,j, min(l, u)))
        }
    }
    q <- list()
    q[['indiv']] = res
    q[['p']] = sum(res[,3]) / (nrow(res) * n.permutations)
    return(q)
}
