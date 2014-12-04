#--------------------------------------------------
# function to automatically load packages
#--------------------------------------------------
packages<-function(x){
    x<-as.character(match.call()[[2]])
    if (!require(x,character.only=TRUE)){
        install.packages(pkgs=x,repos="http://cran.r-project.org")
        require(x,character.only=TRUE)
    }
}

#--------------------------------------------------
#download and load data using wget
#--------------------------------------------------
download_data <- function(){
    if( Sys.info()['sysname']=="Linux" ){
        system("wget http://bergelson.uchicago.edu/Members/mhorton/resources/snps/accession_coordinates.xls")
        system("wget http://bergelson.uchicago.edu/Members/mhorton/resources/snps/snps.tgz")
        system("tar xzvf snps.tgz")
    }
}

#--------------------------------------------------
# function to read all the data and coords
#   this function loads the three required files:
#       - snp_file, the file with the snps
#       - coords_file, the file with sample coords
#       - outgroup_file, file to polarize snps
#
#   will create 2 data frames in the global name space:
#       data: which contains the snp data
#       coords: which contains info on the pop
#           locations
#--------------------------------------------------
read_data_arab <- function(snp_file="athal_snps_031110_agdp",
          coords_file="accession_coordinates.csv",
          outgroup_file=NULL,  nsnp=NULL, ...){

    if(is.null(nsnp)){
	data <- read.table(snp_file, header=T, strings=F,  ...)
    }else{
	data <- read.table(snp_file, header=T, strings=F, nrows=nsnp, ...)
    }

    print("finished reading main data file")

    coords <- read.csv(coords_file, header=T, strings=F)
    print("finished reading other data file")

    data <- polarize_data( data, outgroup_file )

    #remove individuals without coords
    ids <- as.numeric(substring(names(data),2))
    to.remove <- which(!unlist(lapply(ids,function(x)x%in%coords$ecotype_id)))
    data<-data[-to.remove]
    ids <- as.numeric(substring(names(data),2))

    data.ordering<-unlist(lapply(coords$ecotype_id,function(x)which(x==ids)))
    data <<- data[,data.ordering]
    coords$ecotype_id <- paste0("X",coords$ecotype_id)
    coords$id <- coords$ecotype_id

    coords <<- coords
}

get_snp_state <- function(snp_csv_file, snp_annotation_file){
    #file_id <- "camax_sd5484_miss5_1pl_ga"                                
    #a <- read.csv(sprintf("%s.csv",file_id),header=T, strings=F)          
    #snp <- read.table("snp_camax_out.txt")                                
    snp <- read.table(snp_annotation_file)                                
    a <- read.csv(snp_csv_file,header=F, strings=F, nrow=2)          
    aa <- t(a[,-1])
    print("done annotation merging")
    names(snp) <- c("chr","pos","state")                                  
    pos <- as.numeric(aa[,2])                                  
    chr <- aa[,1]
    b <- data.frame(chr, pos, id=1:length(chr))                           
    bb <- merge(b, snp)                                                   
                                                                          
    o <- order(bb[,3])                                                    
                                                                          
    bb <- bb[o,-3]                                                        
                                                                          
    return(bb)
    print("done annotation merging")
}

# reads data from a snapp file. this assumes that the first two
# data columns contain outgroup individuals that can be used
# to polarize the SNP
read_data_snapp <- function(snp_file="athal_snps_031110_agdp",
          coords_file="accession_coordinates.csv", nsnp=NULL, 
          annotation=NULL, outgroup_columns=NULL,...){

    data <- read.table(snp_file, header=F, strings=F,
               sep=",", na.strings="?",...)
    print(c("finished reading main data file", dim(data) ))
    n_outgroups <- length(outgroup_columns)

    if( n_outgroups > 0 ){
	    outgrp <- data[outgroup_columns,-1]
	    data <- data[-1* outgroup_columns,]

	    if ( !is.null(nsnp ) ){
	    outgrp <- outgrp[,1:nsnp]
	    }
    }
    hdr <- data[,1]
    data <- t(data[,-1])
    data <- as.data.frame(data)
    names( data ) <- hdr

    if ( !is.null(nsnp ) ){
    print(sprintf("keeping only %s snp", nsnp))
    data <- data[1:nsnp,]
    }



    if( n_outgroups  > 0 ){
	derived <- colSums(outgrp==2, na.rm=T)== n_outgroups
	anc <- colSums(outgrp==0, na.rm=T)== n_outgroups
	no.outgrp <- is.na(derived + anc ) | derived+anc==0

	data[derived,] <- 2 - data[derived,]
	data <- data[!no.outgrp, ]
    }


    coords <<- read.csv(coords_file, header=T, strings=F)
    print(c("finished reading other data file", dim(data)))

    #remove individuals without coords
    #ids <- as.numeric(substring(names(data),2))
    ids <- names(data)
    to.remove <- !unlist(lapply(ids,
            function(x)x%in%coords$id))
    data<-data[,!to.remove]
    ids <- names(data)

    data.ordering<-unlist(lapply(
        coords$id,function(x)which(x==ids)))
    data <- data[,data.ordering]
    if(!is.null(annotation)){
        annotation <<- annotation[!no.outgrp]
    }
    data <<- data
}

#--------------------------------------------------
# polarizes data using the outgroup data set
#--------------------------------------------------
polarize_data <- function(data, outgroup_file){
    if( is.null( outgroup_file ) ) return( data )

    anc_data <- read.table( outgroup_file, header=T, 
                            strings=F )
    anc_data <<- anc_data

    full_data <- merge(anc_data,data)
    print(c("done merging",ncol(data), ncol(full_data)))

    #sort the data after merging
    chr_name <- names(data)[2]
    pos_name <- names(data)[3]
    o <- order(  get( chr_name, full_data  ), 
         get( pos_name, full_data  ))
    full_data <- full_data[o,] 

    #polarize snps
    first_snp_pos <- ncol( anc_data ) + 1
    snp_pos <- first_snp_pos:ncol(full_data)
    to_switch <- full_data$minor_allele == full_data$lyrata_consensus
    full_data[ to_switch, snp_pos ] <- 1- 
        full_data[ to_switch, snp_pos ]
    full_data <<- full_data

    loc_data <<- data[,1:first_snp_pos]
    data <- full_data[,snp_pos]

    return( data )
}
#--------------------------------------------------
#the two following functions create the data frame
# from a population set. pops is a list, where each
# entry is a vector with the individuals in a given
# population. For example, the list
#
# [[1]]
# 1, 6, 7, 8
#
# [[2]]
# 2, 3, 4, 5
#
# would assign individuals 2-5 to population 2, 
# all other individuals are in population 1
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
make_pop_data_from_pops <- function( pops, data ){
    pop_data  <- sapply(pops,function(x)rowSums(data[x],
                                                na.rm=T))
}

make_pop_ss_from_pops <- function( pops, data, ploidy=2 ){
    count_nas <- function(x){
    if( length(x) > 1 ){
        return (ploidy * rowSums((!is.na(data[,x])))) 
    } else {
        return( ploidy * (!is.na(data[,x]))) 
    }
    }
    pop_ss  <- sapply(pops, count_nas)

}
make_pop_coords_from_pops <- function( pops, coords ){
    pop_coords <- t(sapply(pops,get_coords_for_pop,coords))
    pop_coords <- data.frame(pop_coords,n=sapply(pops,length))  
    pop_coords <- apply(pop_coords,2,unlist)
    pop_coords<-as.data.frame(pop_coords)
    pop_coords[,c('latitude','longitude')] <- 
	apply(pop_coords[,c('latitude','longitude')],2,as.numeric)
    regions <- unlist(unique(pop_coords[,3]))
    pop_coords <- cbind(pop_coords, id=sapply(pops, function(x) coords[x[1],]$id) )
    pop_coords <<- pop_coords
    return(pop_coords)
}

get_coords_for_pop<- function(pop,coords=coords){
    return(coords[pop[1],c('latitude','longitude','pop')])
}

#--------------------------------------------------
#calculates psi for all pairs of populations
#       -pop_data: a m x n matrix, where each row
#           is a SNP, and each column is a pop.
#           Entries are the absolute allele freq
#           at a SNP in each pop
#       -pop_ss: a m x n matrix of number of genotype
#           calls made for a population at a given SNP
#       -pop_coords
#           data frame that gives
#               latitude, longitude, country, 
#               population and sample size for each
#               population
#       -n : the number of samples to downsample to
#       -subset: subset of pops to calculate psi over
#--------------------------------------------------
get_all_psi <- function(pop_data, pop_ss, n=2,
                    subset=1:ncol(pop_data),
                    resampling="hyper"){
#this function calculates psi for columns i,j, both resampled down to
# n samples

    if( is.logical( subset) )
        subset <- which(subset)
    n_pops <- length(subset)
    mat = matrix( 0, nrow=n_pops, ncol=n_pops )
    for(i in 1:(n_pops-1)){
        for(j in (i+1):n_pops){
	    ii <- subset[i]
	    jj <- subset[j]
            ni <- pop_ss[,ii]
            nj <- pop_ss[,jj]
            fi <- pop_data[,ii]
            fj <- pop_data[,jj]
            mat[i,j] <- get_psi( ni, nj, fi, fj, 
                                resampling=resampling, n=n )
            mat[j,i] <- -mat[i,j]
	    print( c(ii, jj))
        }
    }
    
    return(mat)
}

#--------------------------------------------------
#the actual calculation of the psi statistic for
#a single pair of populations. The function requires
#4 vectors, all of length equal to the number of snps
#to be analyzed.
#
#   Parameters:
#       - fi : vector of int
#           - vector of derived allele counts for each 
#               snp in population i
#       - fj : vector of int
#           - vector of derived allele counts for each 
#               snp in population j
#       - ni : vector of int
#           - vector of total number of genotypes in 
#               population i
#       - nj : vector of int
#           - vector of total number of genotypes in 
#               population j
#       -n : the number of samples to downsample to
#--------------------------------------------------
get_psi <- function (ni, nj, fi, fj,
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


        poly_mat <- matrix(0,nrow=n+1, ncol=n+1)
        poly_mat[2:(n+1),2:(n+1)] <- 1
        poly_mat[n+1,n+1] <- 0

        #psi_mat is the contribution to psi for each entry
        psi_mat <- outer(0:n,0:n,FUN=function(x,y)(x-y))
        psi_mat[1,] <- 0
        psi_mat[,1] <- 0

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


        resampled_mat <- matrix(rowSums(apply(tbl, 1,
                                                  f.contribution)),nrow=n+1)
        
        return( sum(resampled_mat * psi_mat) / sum(resampled_mat * poly_mat) )
    }

#this function will group individuals into populations,
# which are then subsequently used for further analysis
make_pops <- function( coords, n=2 ){
    get_pops_w_multiple_samples <- function(coords,cutoff=3){
        #this function returns all locations (in lat, long) that have more than
        #cutoff samples. I'll use these as locations for analysis
        #the return is the id (i.e. col in data, row in coords) of the populations
        
        pop_list = list()
        pop1 <- unlist(lapply(coords[,3],function(x)which(x==unique(coords[,3]))))
        pop2 <- unlist(lapply(coords[,4],function(x)which(x==unique(coords[,4]))))
        pop <- 10000*pop1 + pop2

        ids_to_keep <- as.numeric(names(table(pop))[table(pop)>cutoff])
        for(i in 1:length(ids_to_keep)){
            pop_list[[i]] <- which(pop==ids_to_keep[i])
        }

        return(pop_list)
    }
    get_pops_w_multiple_samples <<- get_pops_w_multiple_samples


    #get all populations from sample data
    pops <- get_pops_w_multiple_samples(coords,n)
}

#makes a pop list with one entry per pop, use this for individual
#based analyses
make_pops_auto <- function( coords ){
    pops <- list()
    for(i in 1:nrow(coords) ){
    pops[[i]] <- c( i )
    }
    pops
}

#------------------------------------------------------------
# the function that infers the estimate of the expansion
# origin for a single group of populations, habitat borders
# are taken from the worldmap
# arguments:
#   region, countries: which region or countries are supposed
#            to be plotted, compared to pop_coords
#   pop_coords, all_psi: data objects
#   xlen, ylen: number of origins to consider in each direction
#           increasing this will reduce speed
#
# this is done using three different functions: the first,
# prep_tdoa_data, is reformatting the data and selecting the
# samples from the requested populations
# - the second, model_countries, goes through the points and
# calls the third function, model_1d, if a point should be considered for the
# origin. model_1d calculates the goodness of fit for a single
# candidate point
#
#------------------------------------------------------------
find_origin <- function(pop_coords, all_psi, 
                        region=NULL, countries=NULL,
                        xlen=50,  ylen=50, doPlot=F, doPdf=F, 
             ...){
    data_countries <- prep_tdoa_data( pop_coords, all_psi, region,
                          countries, xlen, ylen, ... )

    samples <- pop_coords[pop_coords$pop %in% region,]
    bbox <- get_sample_bbox( samples )
    #bbox <- get_country_bbox(data_countries[[2]], pop_coords, ...)

    print(bbox)

    data_countries <<- data_countries
    res <- model_countries(data = data_countries[[1]],
               bbox = bbox,
                           pop_coords,
                           countries= data_countries[[2]],
                           xlen=xlen,
                           ylen=ylen, ...)


    orig <- find_orig_table(bbox, res, xlen, ylen)

    if(doPlot !=F){
        plot_dataset(bbox=bbox, res=res, samples=samples,
                     xlen=xlen, ylen=ylen,
                     orig=orig, doPdf=doPdf )
    }

    return( orig )
}
prep_tdoa_data <- function(pop_coords, all_psi, 
                        region=NULL, countries=NULL,
                        xlen=50,  ylen=50, ...){

    if(is.null(countries)){
	countries <- coords2country( pop_coords[, c('longitude','latitude')])
	pop_coords$country <- as.character(countries)
	countries <- countries[!is.na(countries)]
        countries <- as.character(unique(countries))
        ip <- pop_coords$pop %in% region
    } else {
        ip <- pop_coords$country %in% countries
    }
    print( sum(ip))
    locs <- pop_coords[ip,1:2]
    psi <- all_psi[ip,ip]
    n_locs <- nrow(locs)

    tdoa_data <- c()
    for(i in 1:n_locs){
        for(j in (i+1):n_locs){
            if( i>=j | j >n_locs) break
            tdoa_data <- rbind( tdoa_data, c(locs[i,], locs[j,], psi[i,j]))
            }}
    write.table(tdoa_data, "tmp2", col.names=F, row.names=F)

    tdoa_data <<- tdoa_data

    tdoa_data <- matrix(unlist(tdoa_data), ncol=5)
    print(countries)

    res <- list(tdoa_data, countries)
    return(res)
}

model_countries <- function(data, bbox,  pop_coords,
                countries, pct=0.01,
            xlen=100, ylen=100, ...){


    #define locs for estimate
    s1<-seq(bbox[1,1],bbox[1,2], length.out=xlen)
    s2<-seq(bbox[2,1],bbox[2,2], length.out=ylen)

    # init output
    d0 <- matrix(NA, ncol=ylen, nrow=xlen)
    rsq <- matrix(NA, ncol=ylen, nrow=xlen)
    mdlq <- list()
    for(i in 1:xlen){
        mdlq[[i]] <- list()
    }
    
    #stupid progress bar
    #plot(NA,xlim=c(s1[1],s1[xlen]), ylim=c(s2[1], s2[ylen]))
    
    for( i in 1:length(s1)){
        for( j in 1:length(s2) ){
            x <- s1[i]
            y <- s2[j]
            #print(c(x,y))
            #stupid progress bar
            #points(x,y)

            if( ! is.na(coords2country(data.frame(x,y)))){
                mdl <- model_1d( xy=c(x,y), data=data, pct=pct)
                d0[i, j] <- mdl$f_e
                rsq[i, j] <- mdl$rsq
                mdlq[[i]][[j]] <- mdl
            }
        }
    print(i)
    }
    return( list( d0, rsq, mdlq, bbox))
}


#------------------------------------------------------------
# find the origin from the table generated previously
#------------------------------------------------------------
find_orig_table <- function(bbox, res, xlen, ylen){
    orig <- which(t(res[[2]])==max((res[[1]]>0)*res[[2]],
                                   na.rm=T), arr.ind=T)        
    s1<-seq(bbox[1,1],bbox[1,2],length.out=xlen) 
    s2<-seq(bbox[2,1],bbox[2,2],length.out=ylen) 

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

#------------------------------------------------------------
# make the plot for a data set
#------------------------------------------------------------
plot_dataset<- function(bbox, res, samples, 
            xlen, ylen, orig, doPdf="plot.doPdf"){
    if(doPdf != F) pdf(file=doPdf, width=8, height=6)

    plt <- palette()
    palette(c("#999999",heat.colors(100)))


    plot_map( bbox, res, xlen, ylen)
    points(orig[1],orig[2], col="black", pch="x", cex=2)


    hets <- (samples$hets - min(samples$hets) )/(
        max(samples$hets) - min(samples$hets))
    points( samples$longitude, samples$latitude,
        pch=16, cex=3, col=grey(hets) )
    points( samples$longitude, samples$latitude,
        pch=1, cex=3, col="black",lwd=2 )

    palette(plt)
    if(doPdf != F )dev.off()
}


plot_map <- function( bbox, res, xlen=100, ylen=100){
    m <- getMap("high")
    #m <- map(xlim=bbox[1,], ylim=bbox[2,])
    plot(NA, xlim=bbox[1,], ylim=bbox[2,], xlab="", ylab="", xaxt="n", yaxt="n")

    s1<-seq(bbox[1,1],bbox[1,2],length.out=xlen)
    s2<-seq(bbox[2,1],bbox[2,2],length.out=ylen)
    xc<-c(); yc<-c(); for(i in s1){ for(j in s2){xc <- c(xc, i); yc <-c(yc, j)}}
    rel<- (res[[1]]>0)*(res[[2]]-min(res[[2]],na.rm=T))/(max(res[[2]],na.rm=T)-min(res[[2]],na.rm=T))
    rel <<- rel
    e <- cbind(xc, yc, c(t(rel+0.01)))
    points(e[,1], e[,2], col=e[,3]*100, pch=".", cex=14)

    plot(m, xlim=bbox[1,], ylim=bbox[2,], add =T)
    #rel[is.na(rel)] <- N
    e <<- e
    #points(e[,1], e[,2], col=e[,3]>0, pch=15, cex=1)
    #points(pop_coords[,2], pop_coords[,1], pch=16,cex=1)
}
#------------------------------------------------------------
# for a given coordinate, gives the country the point is in
#------------------------------------------------------------
coords2country = function(points){  
    countriesSP <- getMap(resolution='low')
    pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
    
     
    # use 'over' to get indices of the Polygons object containing each point 
    indices = over(pointsSP, countriesSP)
     
    indices$ADMIN  
}
#------------------------------------------------------------
# function that calculates the regression for a single origin.
# xy are the coordinates of the candidate point, 
# data is the tdoa_data table,
# pct is the percent size for the founder effects
# f.dist is the distance measure function, with options "haversine" and
# "euclidean"
#
# returns an lm object, with the additional attributes 
#   f_e, for the founder effect at the given pct level, and
#   rsq, for the r squared of the stuff
#------------------------------------------------------------
model_1d <- function(xy, data, pct=0.01, f.dist="haversine"){
    if (f.dist=="euclidean"){
    f.dist <- function(i, j){
        sqrt((ix -jx)^2 + (iy-jy)^2 )
    }
    }else{ if(f.dist=="haversine"){
        f.dist <- distHaversine
    }}

    y = xy[2] 
    x = xy[1]
    ixy = data[,2:1]
    jxy = data[,4:3]
    psi = data[,5]


    d = f.dist(ixy, c(x,y)) - f.dist(jxy, c(x,y))
    l = lm( psi ~ d ) 
    l$f_e = .5 * pct/ l$coefficients[2]
    l$rsq = summary(l)$r.squared

    return (l)
}



#------------------------------------------------------------
# bounding box functions; is used to find the proper plotting
#  region.  if bbox_country=T, the whole country is plotted,
# otherwise it is restricted to the region spanned by the 
# samples.
#------------------------------------------------------------
get_country_bbox <- function(country, pop_coords, bbox_country=F){
    get_single_country_bbox_ <- function(country){
        m <- getMap("high")
        q <- m@polygons[m$ADMIN %in% country]
        q2 <- try(q[[1]]@Polygons, silent=T)
        if( class(q2) == "try-error")
        q2 <- try(q[[2]]@Polygons, silent=T)
        crds <- q2[[1]]@coords
        mins <- apply(crds, 2, min)
        maxs <- apply(crds, 2, max)
        return(cbind(mins,maxs))
    }

    get_single_country_bbox <- function(country, pop_coords){
        samples <- pop_coords[pop_coords$country %in% country,]
        mins <- apply(samples[,2:1],2,min)
        maxs <- apply(samples[,2:1],2,max)
        return(cbind(mins, maxs))
    }



    bbox <- cbind(c(1000,1000),c(-1000,-1000))
    for( c in country){
    if(bbox_country){
        new_bbox <- get_single_country_bbox_( c )
    }else{
        new_bbox <- get_single_country_bbox( c, pop_coords )
    }

    bbox[1,1] <- min(bbox[1,1], new_bbox[1,1])
    bbox[1,2] <- max(bbox[1,2], new_bbox[1,2])
    bbox[2,1] <- min(bbox[2,1], new_bbox[2,1])
    bbox[2,2] <- max(bbox[2,2], new_bbox[2,2])
    }
    bbox
    bbox <<- bbox
}

get_sample_bbox <- function(samples){
        mins <- apply(samples[,2:1],2,min)
        maxs <- apply(samples[,2:1],2,max)
        return(cbind(mins, maxs))
}

#------------------------------------------------------------
# function to analyze a given region. This function will run
# the analysis for all samples in a given region, and append
# to the table res.tbl the founder effect sizes. If doPlot
# is set to true, the function will also produce a plot with
# the most likely origin of the expansion.

# arguments:
#   pop_coords, all_psi are the data set to analyze
#   region gives the region to pick samples from
#   xlen, ylen are the number of points to evaluate
#   doPlot: should a plot be generated?
#   doPdf: should the plot be saved to a pdf?


#   returns: none; however, the table res.tbl is appended a
#   line with the output for the region
#------------------------------------------------------------
run_region <- function( region, loc_file_id=out_file_id, 
                       xlen=100, ylen=100,...){
    reg.str <- paste(region, collapse="+")
    print(c("fid", loc_file_id))
    dir.create('plots', showWarnings = F)
    pdf.name<- sprintf("plots/orig_%s_%s_%d.pdf", reg.str, loc_file_id, 
                       xlen)
    print(pdf.name)
    origin1 <- find_origin(pop_coords, all_psi,
                  region=region,
                   xlen=xlen,ylen=ylen,
                   doPlot=T, doPdf=pdf.name)
    res.line<- cbind(region=reg.str, origin1)
    res.tbl <<- rbind( res.tbl, res.line )
}



do_pca <- function(pop_data, pop_coords){
    plt <- palette()
    palette(c("black", "pink", "blue", "brown", "black", "cyan", "lightblue",
              "red", "orange", "darkgreen","darkgreen","darkgreen"))

    pop_data_n <- apply(pop_data,1,function(x)x/pop_coords$n) 
    p <- prcomp(pop_data_n[,1:100000]) 
    plot(p$x, col=pop_coords$pop, pch=15)
    palette(plt)
}


fst2 <- function(a, b){
    h_1 <- mean(a==1)
    h_2 <- mean(b==1)
    h_b <-mean( a/2 * (1-b/2) + b/2 * (1-a/2) )
    h_w <- mean(h_1 + h_2 ) /2

    fst <-  1- h_w / h_b
    return( c(fst, h_b, h_w,h_1, h_2) )
}

get_all_pairwise_fst2 <- function( pop_data ){
    n <- ncol( pop_data )
    fst.mat <- sapply(1:n, function(x) sapply(1:n, function(y) 
                fst(pop_data[,x], pop_data[,y])[1]))
}

fst<-function(x,y)calculate.pairwise.Fst(t(pop_data[,c(x,y)]),t(pop_ss[,c(x,y)]))

get_all_pairwise_fst<- function( pop_data, pop_ss ){
    packages( BEDASSLE )
    n <- ncol( pop_data )
    sapply(1:n, function(x) sapply(1:n, function(y) fst(x,y)) )
}

get_heterozygosity <- function( pop_data ){
    hets <- colMeans(pop_data==1, na.rm=T)
    hets
}
