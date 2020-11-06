## Try to obtain some quality check on the intron orthology in order to
## make sure that the observations from the sequence alignments is OK.

## in case they come in useful;
source('../R_common_functions/functions.R') ## for ensembl functions
source("../R_intron_alignments_2/functions.R")
require("RMySQL")
require("parallel")

dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")


## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2

## the max.nm allowed
max.l <- 1e8
seed.sp <- c("danio rerio", "gasterosteus aculeatus", "takifugu rubripes", "mus musculus")

## the matrix used to define the alignments
sub.matrix <- make.sub.matrix()
sub.matrix.2 <- make.sub.matrix( letters=c('A', 'C', 'G', 'T'), match=c(4,4,4,4) )

## the exon alignments from which the intron orthology was defined:
ex.align.id <- readRDS("../R_172_genomes/ex_align_id.rds")

## for further analyses I will also need:
gene.stats <- read.table("../family_members/orthologue_transcripts/exon_intron_stats.csv", sep="\t", stringsAsFactors=FALSE)
colnames(gene.stats) <- c('sp', 'db', 'family', 'gene', 'transcript', 'chr', 'strand', 'pos', 'ex.s', 'in.s')

dim(gene.stats)
## [1] 1063202      10
## about a million rows

exon.f <- list.files("../family_members/orthologue_transcripts/", pattern="exon",  full.names=TRUE)
exon.f <- exon.f[!grepl("stats\\.csv", exon.f)]  ## contains the stats file
names(exon.f) <- sub(".+?exons_(.+)$", "\\1", exon.f)
## 6114 files


## 6114 entries
##   4  (no names unfortualy, but sp1 is, 'danio.rerio', 'g aculeatus', 'takifugu rubripes' 'mus'
##      8: sp.i, sp.1, sp.2, al, id1, id2, tr1, tr2
##        al: stats, seq, a.pos, b.pos, al.i

## we need to make a function that takes two aligned sequences and calculates
## a sliding window alignment score and then we can look at how alignments at
## positions

###################################
##################################
## test the local score function.
plot( local.score( ex.align.id[[1]][[1]]$al[[1]]$seq, 10L, -8L, sub.matrix ), type='l' )
plot( local.score( ex.align.id[[1]][[1]]$al[[2]]$seq, 10L, -8L, sub.matrix ), type='l' )

## that seems to do the correct thing.
h <- 8
for(i in 1:length( ex.align.id[[h]][[1]]$al)){
    plot( local.score( ex.align.id[[h]][[1]]$al[[i]]$seq, 10L, -8L, sub.matrix.2 ), type='l' )
    abline(h=0, lty=2)
    with( ex.align.id[[h]][[1]]$al[[i]], segments( a.pos, 0, a.pos, 4, col='red', lwd=2))
    with( ex.align.id[[h]][[1]]$al[[i]], if(length(b.pos)) segments( b.pos, -8, b.pos, 0, col='blue', lwd=2))
    input <- readline(paste(i, ":"))
}



####################################
## Seems OK, with some worryingly bad alignments..
## but at least there is some contrast..

## then we try to make a similar orthology tables that we created before, but
## with the alignment scores at the intron sites.
dbs <- strsplit(readLines( "../family_members/vertebrate_family_members_1_130.txt", n=1 ), "\t")[[1]]
names(dbs) <- sub("_core_98_\\d+", "", dbs )
sp.names <- sub("_", " ", names(dbs))


## this is very slow.. 
intron.orth <- lapply(1:length(seed.sp), function(i){
    tmp <- lapply(ex.align.id, function(x){
        ## here I will only consider introns which could be aligned to the seed species
        ## one could try to make this for the full set of seed species. Then one could
        ## also consider where it differs. But for now we will leave it as it is..
        m.i <- matrix(ncol=length(sp.names), nrow=length( x[[i]]$al[[1]]$a.pos ) ) ## the indice
        colnames(m.i) <- sp.names
        m.id <- m.i  ## the genes
        t.id <- m.i ## the transcripts
##        m.l <- m.i   ## the lengths
        m.scores <- m.i  ## localised intron scores
        scores <- rep(0, length(sp.names))
        names(scores) <- sp.names
        sp1 <- x[[i]]$sp.1
        id1 <- x[[i]]$id1
        if(is.null(x[[i]]$al))
            return( list(m.i, m.id, t.id, m.scores) )
        for(j in 1:length(x[[i]]$al)){
            sp2 <- x[[i]]$sp.2[j]
            id2 <- x[[i]]$id2[j]
            tr2 <- x[[i]]$tr2[j]
            y <- x[[i]]$al[[j]]
            if(scores[ sp2 ] < y$stats['score'] ){
                local.scores <- local.score( y$seq, 10L, -8L, sub.matrix.2 )
                m.i[ , sp2 ] <- NA  ## overwrite
                m.i[ y$al.i[,1] + 1, sp2 ] <- y$al.i[,2] + 1
                m.id[ , sp2 ] <- id2
                t.id[ , sp2 ] <- tr2
                m.scores[ , sp2 ] <- local.scores[ y$a.pos + 1 ]
                scores[ sp2] <- y$stats['score']
            }
        }
        list(m.i, m.id, t.id, m.scores)
    })
    names(tmp) <- names(exon.f)
    tmp
})
names(intron.orth) <- seed.sp

## and create full tables from these:
## create a table of orthologous introns from these
intron.orth.t <- vector(mode='list', length=length(seed.sp))
names(intron.orth.t) <- seed.sp
for(i in 1:length(intron.orth.t)){
    intron.orth.t[[i]] <- vector(mode='list', length=4)
    names(intron.orth.t[[i]]) <- c('i', 'id', 'tr', 'score')
    intron.orth.t[[i]][[1]] <- do.call( rbind, lapply( intron.orth[[i]], function(x){ x[[1]] }) )
    intron.orth.t[[i]][[2]] <- do.call( rbind, lapply( intron.orth[[i]], function(x){ x[[2]] }) )
    intron.orth.t[[i]][[3]] <- do.call( rbind, lapply( intron.orth[[i]], function(x){ x[[3]] }) )
    intron.orth.t[[i]][[4]] <- do.call( rbind, lapply( intron.orth[[i]], function(x){ x[[4]] }) )
}


## we can have a look at the distributions of the intron scores for different species
for(j in 1:ncol( intron.orth.t[[1]]$i )){
    with( intron.orth.t[[1]], hist( score[ !is.na( i[,j] ), j ], main=colnames(intron.orth.t[[1]]$i)[j]) )
    input <- readline(paste(j, ":"))
}

saveRDS( intron.orth.t, 'intron_orth_lscore.rds' )
