## Extend the sampled set of intron alignments in order to increase
## the coverage at longer mn reaches.
##
## Otherwise follow the procedure established for:
## ../R_intron_alignments_2/extracted_seqs/blast/blast_alignments.R
##
## Specifically in order to produce a better:
##
## alignment_scores_filtered.pdf
## summary figure

## run alignments of long introns in parallel

library(parallel)

dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")
source("../R_intron_alignments_2/functions.R")

max.l <- 1e8
## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2


## this makes use of a number 
## of global variables; use with care.. 
## max.l is the maximum length of the score tables
## i is the row int fam.id / orth.i / from which
## we extract b.seq
## if j == i then we align b.seq to the orthologous
## intron is al.sp; otherwise as a control we align to
## a control intron chosen to have a simlar size.
align.introns <- function(i, j=i, al.sp='danio.rerio',
                          gap=c(-10,-1), min.width=15, min.score=40,
                          sm=sub.matrix, max.l=max.l, r.c=FALSE){
    fam.id <- orth.fam[i]
    tr <- ortho$tr[i,]
    int.i <- ortho$i[i,]
    
    int.seq <- read.exons( intron.f[ fam.id ] )
    int.meta <- sapply( int.seq, function(x){ c('id'=x$id, 'tr'=x$tr, 'sp'=x$sp, 'n'=length(x$l) )})
    ## select the sequences that we will use
    b <- int.meta['tr',] %in% tr
    int.seq <- int.seq[b]
    int.meta <- int.meta[,b]
    ## so we can order by species
    names(int.seq) <- sub(" ", ".", int.meta['sp',])
    colnames(int.meta) <- sub(" ", ".", int.meta['sp',])
    
    int.i <- int.i[ names(int.i) %in% names(int.seq) ]
    al.i <- which( names(int.i) == al.sp )
    if(i == j){
        a.seq <- int.seq[[al.sp]]$e[ as.numeric(int.i[al.i]) ]
    }else{
        a.int.seq <- read.exons( intron.f[ orth.fam[j] ] )
        a.int.meta <- sapply( a.int.seq, function(x){ c('id'=x$id, 'tr'=x$tr,
                                                        'sp'=x$sp, 'n'=length(x$l) )})
        a.seq.b <- a.int.meta['tr',] == ortho$tr[j,al.sp]
        a.seq <- a.int.seq[[ which(a.seq.b) ]]$e[ ortho$i[j,al.sp] ]
    }
    ##        
    b.seq <- sapply( names(int.i[-al.i]), function(sp){
        int.seq[[sp]]$e[ as.numeric(int.i[sp]) ] })
    b.seq <- b.seq[ !is.na(b.seq) ]
    b.sp <- names(b.seq)
    ## reverse complement if r.c
    if(r.c){
        b.seq <- rev.comp(b.seq)
        names(b.seq) <- b.sp
    }
    ##
    aligns <- vector(mode='list', length=length(b.seq))
    names(aligns) <- names(b.seq)
    for(sp in names(b.seq)){
        if(as.double(nchar(a.seq)) * as.double(nchar(b.seq[sp])) <= max.l )
            aligns[[sp]] <- local.aligns( a.seq, b.seq[sp], sm$offset,
                                        sm$size, sm$sm,
                                        gap=gap, min.width=min.width, min.score=min.score )
    }
    list(sp1=al.sp, sp=names(b.seq), al=aligns, l1=nchar(a.seq), l2=unname(sapply(b.seq, nchar)),
         i=i, j=j, al.i=al.i, seq.meta=int.meta )
}


align.introns.ctl <- function( i, max.l, gap, ... ){
    ## find a suitable j to make the alignment from
    rank <- ga.int.l.r[ i ]
    rank.r <- -1:1 + rank
    rr.i <- ga.int.l.or[ rank.r ]
    rr.i <- rr.i[!is.na(rr.i) & rr.i != i]
    align.introns(i=i, j=rr.i[1], max.l=max.l, gap=gap, ...)
}


conv.sp.names <- function(x){
    x <- sub(".", "_", x, fixed=TRUE)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}


## variance in teleosts and 
tel.var <-  read.table("../R_172_genomes/dr_teleost_var.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
colnames(tel.var) <- sub("\\.$", "", colnames(tel.var)) 


## the files containing intron sequences
intron.f <-  read.table("../R_172_genomes/intron_files.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
tmp <- rownames(intron.f)
intron.f <- intron.f[,1]
names(intron.f) <- tmp
rm(tmp)

## We will probably need these later on. At the very minimum we will need the family
## one to obtain ortholog sequences

ortho.f <- list.files("../R_172_genomes/", pattern="dr.+?_ortho.+", full.names=TRUE)
names(ortho.f) <- c('i', 'id', 'l', 'tr')
ortho <- lapply(ortho.f, read.table, header=TRUE, sep="\t", stringsAsFactors=FALSE)

orth.fam <- read.table("../R_172_genomes/dr_intron_fam.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)[[1]]

## for controls I need the following
ga.int.l <- ortho$l[,'danio.rerio']
ga.int.l[ is.na(ga.int.l) ] <- 0
ga.int.l.o <- order( ga.int.l )
ga.int.l.r <- rank( ga.int.l, ties.method='random' )
ga.int.l.or <- (1:length(ga.int.l))[ga.int.l.o]


## let us read in the previous sampled alignment in order to get the row indices
## from that. That way we can subsample the others..

align.f <- list.files("../R_intron_alignments/", pattern="rds", full.names=TRUE)
names(align.f) <- c('sampled', 'sampled.ctl', 'long', 'long.ctl')

align.d <- lapply(align.f, readRDS)

align.rows <- lapply( align.d, function(x){
    sapply(x, function(y){ y$i }) })

## Define a set of control (sampled) rows:
ctl.b <- tel.var[,'X50'] < 8 & ortho$l[,'danio.rerio'] > 1024 & tel.var[,'n'] > 10
ctl.i <- which(ctl.b)
length(ctl.i)
## 3130

sum( align.rows[['sampled']] %in% ctl.i )
## 500. That is good.
ctl.i <- setdiff( ctl.i, align.rows[['sampled']] )
length(ctl.i) ## 2630

## Since I am going to be spending a few days writing, I might
## as well run alignments for all of these
## unfortunately, odin is under rather heavy load these days
## so things are likely to take a bit longer than otherwise.

sub.matrix <- make.sub.matrix()
gap <- c(-10, -2)

tmp <- align.introns(ctl.i[1], al.sp='danio.rerio', gap=gap, max.l=max.l )
## seems OK
rm(tmp)

ctl.aligns <- mclapply( ctl.i, align.introns, al.sp='danio.rerio', gap=gap, max.l=max.l, mc.cores=10 )

### since this doesn't actually take that much time let us also expand the selected set
### let us call this long.2 (bad, but... )
long.b <- tel.var[,'min'] > 10 & !is.na(tel.var[,'sd']) & tel.var[,'n'] > 10
sum(long.b) ## 1255
long.i <- which(long.b)
long.o <- order(tel.var[ long.b, 'sd' ])

## we have already done the first 500, so do the rest of the alignments here:
long.2.aligns <- mclapply( long.i[ long.o[501:length(long.o)] ], align.introns, al.sp='danio.rerio',
                          gap=gap, max.l=max.l, mc.cores=10 )

saveRDS( long.2.aligns, "long_2_aligns.rds" )

### We can also run sets for shorter sequences:
### we call this long.b.3 since I used long.2 for the second set of long
### alignments. Not good, 
long.b.3 <- tel.var[,'min'] > 8 & !is.na(tel.var[,'sd']) & tel.var[,'n'] > 10
sum(long.b.3)  ## 6323
long.i.3 <- which(long.b.3)
sum( long.i.3 %in% long.i ) ## 1255
long.i.3 <- setdiff( long.i.3, long.i )
length(long.i.3)
## 5068.. This should have min lengths between 256 and 1024.
## do not order these, but simply run alignments for the first 2000
## of these.

long.3.aligns <- mclapply( long.i.3[1:2000], align.introns, al.sp='danio.rerio',
                           gap=gap, max.l=max.l, mc.cores=10 )

saveRDS( long.3.aligns, "long_3_aligns.rds" )
save.image()

## have a look at these:
sp.col.3 <- readRDS("../R_trees_distances/sp_col_3.rds")
names(sp.col.3) <- sub("_", ".", names(sp.col.3))

for(i in 1:length(long.2.aligns)){
    plot.alignments(long.2.aligns[[i]])
    input <- readline( paste(i, ":"))
    if(input == 'q')
        break
}

for(i in 1:length(long.3.aligns)){
    plot.alignments(long.3.aligns[[i]])
    input <- readline( paste(i, ":"))
    if(input == 'q')
        break
}


## extract scores:
## a simple species classification
sp.class <- get.sp.taxonomy()

ctl.aligns.top.cl <- lapply( colnames(sp.class), function(cl){
    sp.b <- sp.class[,cl]
    names(sp.b) <- rownames(sp.class)
    extract.top.al.l( ctl.aligns, max.l=max.l, sp.b=sp.b )
})
names(ctl.aligns.top.cl) <- colnames(sp.class)

saveRDS(ctl.aligns, 'ctl_aligns.rds')
saveRDS(ctl.aligns.top.cl, 'ctl_aligns_top_cl.rds')


long.2.aligns.top.cl <- lapply( colnames(sp.class), function(cl){
    sp.b <- sp.class[,cl]
    names(sp.b) <- rownames( sp.class )
    extract.top.al.l( long.2.aligns, max.l=max.l, sp.b=sp.b )
})
names(long.2.aligns.top.cl) <- colnames(sp.class)
## remember that extract.top.al.l has a bug in it that results
## in some returned elements not having the same structure. That
## gives some weird situations, where we have empt sequences.

saveRDS( long.2.aligns.top.cl, 'long_2_aligns_top_cl' )

## to extract the best aligned sequences we do as before:
for( cl in names(ctl.aligns.top.cl) ){
    fname <- paste("extracted_seqs/", cl, "_ctl", ".fasta", sep="")
    lines <- do.call(c,
                     lapply(ctl.aligns.top.cl[[cl]], function(x){
                             if(!length(x)) return(NULL)
                             c(paste(c(">", ortho$id[ x[[1]]$i, x[[1]]$sp1],
                                       "_", ortho$i[ x[[1]]$i, x[[1]]$sp1],
                                       " ", ortho$tr[ x[[1]]$i, x[[1]]$sp1],
                                       " ", ortho$id[ x[[1]]$i, x[[1]]$sp2],
                                       "_", ortho$i[ x[[1]]$i, x[[1]]$sp2],
                                       " ", ortho$tr[ x[[1]]$i, x[[1]]$sp2],
                                       " ", x[[1]]$sp1, " ", x[[1]]$sp2, " ",
                                       paste(x[[1]]$coords, collapse=" ")), collapse=""),
                               degap.seq( x[[1]]$seq[1] ))
                     }))
    writeLines(lines, fname)
}

## and the same for the long.2 set of alignments:
for( cl in names(long.2.aligns.top.cl) ){
    fname <- paste("extracted_seqs/", cl, "_long_2", ".fasta", sep="")
    lines <- do.call(c,
                     lapply(long.2.aligns.top.cl[[cl]], function(x){
                             if(!length(x)) return(NULL)
                             c(paste(c(">", ortho$id[ x[[1]]$i, x[[1]]$sp1],
                                       "_", ortho$i[ x[[1]]$i, x[[1]]$sp1],
                                       " ", ortho$tr[ x[[1]]$i, x[[1]]$sp1],
                                       " ", ortho$id[ x[[1]]$i, x[[1]]$sp2],
                                       "_", ortho$i[ x[[1]]$i, x[[1]]$sp2],
                                       " ", ortho$tr[ x[[1]]$i, x[[1]]$sp2],
                                       " ", x[[1]]$sp1, " ", x[[1]]$sp2, " ",
                                       paste(x[[1]]$coords, collapse=" ")), collapse=""),
                               degap.seq( x[[1]]$seq[1] ))
                     }))
    writeLines(lines, fname)
}


## blast is now running for these; at this point it makes more sense to import the
## resulting data sets into
## ../R_intron_alignments_2/extracted_seqs/blast/blast_alignments.R
## and modify the figure appropriately.

