## test compiled functions..
source("../../R_intron_alignments_2/functions.R")

bl.data.3 <- readRDS('../../R_intron_alignments_summary_2/bl_data_3.rds')
bl.files.3 <- readRDS('../../R_intron_alignments_summary_2/bl_files_3.rds')
aligns.short.2.top <- readRDS('../../R_intron_alignments_summary_2/aligns_short_2_top.rds')
aligns.top <- readRDS('../../R_intron_alignments_summary_2/aligns_top.rds')

bl.cov.3 <- readRDS('../../R_intron_alignments_summary_2/bl_cov_3.rds')

ql.tel <- sapply( aligns.short.2.top, function(x){
    nchar( degap.seq( x$teleostei$seq[[1]][1] )) })

i <- which(bl.files.3$cl == 'teleostei' & bl.files.3$sp == 'danio.rerio')
## 27..
## which is also the one with the largest numbers of rows; by far...

dyn.load( 'bl_cov.so' )

bl <- bl.data.3[[i]]

tmp <- .Call( 'query_covs', ql.tel,
              bl$qseqid, bl$qlen, bl$qstart, bl$qend, bl$sseqid, "ALT" )

tmp.t <- sapply( 1:length(tmp), function(i){
    all(tmp[[i]] == bl.cov.3$teleostei$danio.rerio[[i]])
})

all(tmp.t) ## TRUE. So that is rather good.

system.time(
    tmp <- .Call( 'query_covs', ql.tel,
                 bl$qseqid, bl$qlen, bl$qstart, bl$qend, bl$sseqid, "ALT" )
)

##  user  system elapsed 
## 6.285   0.100   6.386 

## we can also do that for all of the data, but we need to derive the ql and so on

system.time( 
    tmp.2 <- lapply( c('teleostei', 'mammalia'), function(cl){
        ql <- as.integer(sapply( aligns.short.2.top, function(x){
            nchar( degap.seq( x[[cl]]$seq[[1]][1] )) }) )
        i <- which( bl.files.3$cl == cl )
        lapply(i, function(j){
            bl <- bl.data.3[[j]]
            cat("calling query_covs for ", j, " : ", bl.files.3$sp[j], "\n")
            .Call( 'query_covs', ql,
                  bl$qseqid, bl$qlen, bl$qstart, bl$qend, as.character(bl$sseqid), "ALT" )
        })
    })
)
##   user  system elapsed 
## 11.964   0.136  12.101 

## in R this would have taken several hours. Here it takes 12 seconds.
## I do like C, even if does take a while.. 
