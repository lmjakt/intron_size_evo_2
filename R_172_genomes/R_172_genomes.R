## This extends analysis done on my workstation to included alignments of
## sequences from 172 genomes..
## Oriiginal analysis in ../R_Mosjoen/ on my workstation


source("../R/functions.R")
source("~/R/general_functions.R")

gene.stats <- read.table("../family_members/orthologue_transcripts/exon_intron_stats.csv", sep="\t", stringsAsFactors=FALSE)
colnames(gene.stats) <- c('sp', 'db', 'family', 'gene', 'transcript', 'chr', 'strand', 'pos', 'ex.s', 'in.s')

dim(gene.stats)

## [1] 1063202      10
## about a million rows


t.lengths <- sapply( gene.stats[,'ex.s'], function(x){
    sum( as.numeric( strsplit( x, ',' )[[1]] ) ) })

sum(t.lengths > 10000) / length(t.lengths)
## [1] 0.006477603

sum(t.lengths > 15000) / length(t.lengths)
## [1] 0.0007223463

## how many families have a maximum length
## could probably have done this using the t.lengths structure.. 
fam.tlengths <- tapply( 1:nrow(gene.stats), gene.stats[,'family'], function(i){
    sapply( gene.stats[i,'ex.s'], function(x){
        sum( as.numeric( strsplit( x, ',' )[[1]] )) })
})

fam.max.tlengths <- sapply( fam.tlengths, max )
sum( fam.max.tlengths > 10000 ) ## 506
sum( fam.max.tlengths > 15000 ) ## 73
sum( fam.max.tlengths > 15000 ) / length( fam.max.tlengths) ## [1] 0.01193981

## for the alignments we may want to reduce the number of families anyway,
## but lets look at what we want to do here.

dbs <- strsplit(readLines( "../family_members/vertebrate_family_members_1_130.txt", n=1 ), "\t")[[1]]
names(dbs) <- sub("_core_98_\\d+", "", dbs )
sp.n <- names(dbs)

length(dbs)
## 172

## I don't have the dbs installed on odin, and I have not set them up for access, so I
## will not actually use them here.

fam.sizes <- tapply( 1:nrow(gene.stats), gene.stats$family, length )
sp.gene.n <- tapply( 1:nrow(gene.stats), gene.stats$sp, length )
hist(sp.gene.n)

par('mar'=c(8.1, 4.1, 4.1, 2.1))
barplot( sp.gene.n, names.arg=sub("_[^_]+$", "", names(sp.gene.n)), las=2)

exon.s <- lapply( strsplit( gene.stats$ex.s, "," ), as.numeric )
intron.s <- lapply( strsplit( gene.stats$in.s, ","), as.numeric )

names(exon.s) <- gene.stats$gene
names(intron.s) <- gene.stats$gene

dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")
### for general alignments we want to run a set of seed species against everything else
seed.sp <- c("danio rerio", "gasterosteus aculeatus", "takifugu rubripes", "mus musculus")
sub.matrix <- make.sub.matrix()

## Try to get orthologous introns..
## First read in the exon sequence data.
exon.f <- list.files("../family_members/orthologue_transcripts/", pattern="exon",  full.names=TRUE)
exon.f <- exon.f[!grepl("stats\\.csv", exon.f)]  ## contains the stats file
names(exon.f) <- sub(".+?exons_(.+)$", "\\1", exon.f)
## 6114 files

intron.f <- list.files("../family_members/orthologue_transcripts", pattern="intron",  full.names=TRUE)
intron.f <- intron.f[!grepl("stats\\.csv", intron.f)]  ## contains the stats file
names(intron.f) <- sub(".+?introns_(.+)$", "\\1", intron.f)


## guard against stupid memory usage:
max.slen = 1.5e4  ## this loses about 0.7 percent of families.
## and we should not need more than about 2GB of memory for doing the alignment
##
### that sort of looks ok. but one intron looks a bit wrong..
### well,lets do all of them...

## This will run around, 4 * 172 * 6114 alignments
## i.e. ~ 420,000
## It will return the sequences for each of those alignments
## in addition to other things. That is a minimum of
## sum(t.lengths) * 2 * 4 / 1e9
## 21.2 Gbytes. In reality we can probably double that
## for terminal gaps...
## So we are talking about 40 Gbytes to hold these alignments. I can do this
## on odin, but I really should rewrite the function to return cigar strings and
## add a function to keep the dat
ex.align <- lapply( exon.f, function(f){
    tmp <- read.exons(f);
    tmp.sp <- sapply( tmp, function(x){ x$sp })
    print(paste("processing", f))
    seed.al <- vector(mode='list', length=length(seed.sp))
    if( max( sapply(tmp, function(x){ nchar(x$s) })) > max.slen )
        return(seed.al)

    for(k in 1:length(seed.sp)){
        seed <- seed.sp[k]
        i <- which( tmp.sp == seed )[1] ## take the first one only. This
        if(!is.na(i)){
            sp.index <- cbind( rep(i, length(tmp)), 1:length(tmp))
            b.seq <- sapply( tmp, function(x){ x$s })
            al <- align.seqs.mt( tmp[[i]]$s, b.seq, al.offset=sub.matrix$offset, al.size=sub.matrix$size,
                                 sub.matrix=sub.matrix$sm, gap=as.integer(c(-10, -1)), tgaps.free=TRUE, sp.char='I', thread.n=20 )
            seed.al[[k]] <- list('sp.i'=sp.index, 'sp.1'=seed.sp[k], 'sp.2'=tmp.sp, 'al'=al)
        }
    }

    seed.al
})
##
##
## 20:
##    user  system elapsed 
## 849.470  43.801  52.730
## 849.470/52.730  ## 16x speedup
## 10:
##    user  system elapsed 
## 812.533  32.056  90.642 
## 812.533 / 90.642 ## 8.9x speedup
## Note that the first time I ran this though I did not seem to get the same CPU utilisation.
## I do not know why.. 
save.image()

### Unfortunately I forgot to include the sequence ids which I need for further steps. If I go through the data I should
### be able to enter them into a similar structure.
## this was
## ex.align.id <- ... before, but that isn't so useful.
ex.align.g.id <- lapply( exon.f, function(f){
    tmp <- read.exons(f);
    tmp.sp <- sapply( tmp, function(x){ x$sp })
    tmp.id <- sapply( tmp, function(x){ x$id })
    tmp.tr <- sapply( tmp, function(x){ x$tr })
    print(paste("processing", f))
    seed.al <- vector(mode='list', length=length(seed.sp))
    if( max( sapply(tmp, function(x){ nchar(x$s) })) > max.slen )
        return(seed.al)
    for(k in 1:length(seed.sp)){
        seed <- seed.sp[k]
        i <- which( tmp.sp == seed )[1] ## take the first one only. This
        if(!is.na(i)){
            sp.index <- cbind( rep(i, length(tmp)), 1:length(tmp))
            ## b.seq <- sapply( tmp, function(x){ x$s })
            ## al <- align.seqs.mt( tmp[[i]]$s, b.seq, al.offset=sub.matrix$offset, al.size=sub.matrix$size,
            ##                      sub.matrix=sub.matrix$sm, gap=as.integer(c(-10, -1)), tgaps.free=TRUE, sp.char='I', thread.n=20 )
##            seed.al[[k]] <- list('sp.i'=sp.index, 'sp.1'=seed.sp[k], 'sp.2'=tmp.sp, 'al'=al)
            seed.al[[k]] <- list('sp.i'=sp.index,'id1'=tmp[[i]]$id, 'id2'=tmp.id,
                                 'tr1'=tmp[[i]]$tr, 'tr2'=tmp.tr)
        }
    }
    seed.al
})

for(i in 1:length(ex.align)){
    for(j in 1:length(ex.align[[i]])){
        if(is.null(ex.align[[i]][[j]]))
            next
        if(all(ex.align[[i]][[j]]$sp.i == ex.align.g.id[[i]][[j]]$sp.i)){
            ex.align.id[[i]][[j]] <- c(ex.align[[i]][[j]], ex.align.g.id[[i]][[j]][
                                                               c('id1', 'id2', 'tr1', 'tr2') ])
        }else{
            warning(paste("sp.i does not match for", i, "and", j))
            break
        }
    }
}

saveRDS( ex.align.id, "ex_align_id.rds" )

## create a datastructure holding information about orthologous introns
sp.names.lc <- sub("_", " ", sp.n)
sp.names <- sp.names.lc
substr(sp.names, 1, 1) <- toupper( substr(sp.names, 1, 1))

intron.orth <- lapply(1:length(seed.sp), function(i){
    tmp <- lapply(ex.align.id, function(x){
        ## here I will only consider introns which could be aligned to the seed species
        ## one could try to make this for the full set of seed species. Then one could
        ## also consider where it differs. But for now we will leave it as it is..
        m.i <- matrix(ncol=length(sp.names), nrow=length( x[[i]]$al[[1]]$a.pos ) ) ## the indice
        colnames(m.i) <- sp.names.lc
        m.id <- m.i  ## the genes
        t.id <- m.i ## the genes
        m.l <- m.i   ## the lengts
        scores <- rep(0, length(sp.names))
        names(scores) <- sp.names.lc
        sp1 <- x[[i]]$sp.1
        id1 <- x[[i]]$id1
        if(is.null(x[[i]]$al))
            return( list(m.i, m.id, m.l) )
        for(j in 1:length(x[[i]]$al)){
            sp2 <- x[[i]]$sp.2[j]
            id2 <- x[[i]]$id2[j]
            tr2 <- x[[i]]$tr2[j]
            y <- x[[i]]$al[[j]]
            if(scores[ sp2 ] < y$stats['score'] ){
                m.i[ , sp2 ] <- NA  ## overwrite
                m.id[ , sp2 ] <- NA  ## overwrite
                t.id[ , sp2 ] <- NA 
                m.l[ , sp2 ] <- NA  ## overwrite
                m.i[ y$al.i[,1] + 1, sp2 ] <- y$al.i[,2] + 1
                m.id[ , sp2 ] <- id2
                t.id[ , sp2 ] <- tr2
                m.l[ y$al.i[,1] + 1, sp2 ] <- intron.s[[ id2 ]][ y$al.i[,2] + 1 ]
                scores[ sp2] <- y$stats['score']
            }
        }
        list(m.i, m.id, m.l, t.id)
    })
    names(tmp) <- names(ex.align)
    tmp
})


## create a table of orthologous introns from these
intron.orth.t <- vector(mode='list', length=length(seed.sp))
names(intron.orth.t) <- seed.sp
for(i in 1:length(intron.orth.t)){
    intron.orth.t[[i]] <- vector(mode='list', length=4)
    names(intron.orth.t[[i]]) <- c('i', 'id', 'l', 'tr')
    intron.orth.t[[i]][[1]] <- do.call( rbind, lapply( intron.orth[[i]], function(x){ x[[1]] }) )
    intron.orth.t[[i]][[2]] <- do.call( rbind, lapply( intron.orth[[i]], function(x){ x[[2]] }) )
    intron.orth.t[[i]][[3]] <- do.call( rbind, lapply( intron.orth[[i]], function(x){ x[[3]] }) )
    intron.orth.t[[i]][[4]] <- do.call( rbind, lapply( intron.orth[[i]], function(x){
        if(length(x) == 4)
            return(x[[4]])
        x[[2]] }) )
}


## This is rather waseteful, but in order to directly find which family a particular entry in the
## intron.orth.t table correspond to we need something like this.
intron.orth.fam <- lapply( intron.orth, function(x){
    unlist( lapply( 1:length(x), function(i){
        rep( names(x[i]), nrow(x[[i]][[1]]))
    }))
})
names(intron.orth.fam) <- seed.sp

tmp <- read.table('genome_sizes.txt' )
genome.sz <- tmp[,1]
names(genome.sz) <- rownames(tmp)
genome.sz <- sort(genome.sz)
rm(tmp)

## order the columns by genome size, so that evertying follows afterwards..
o <- match(names(genome.sz), sp.n)
for(i in 1:length(intron.orth.t)){
    for(j in 1:length(intron.orth.t[[i]])){
        intron.orth.t[[i]][[j]] <- intron.orth.t[[i]][[j]][,o]
    }
}


### I can now now calculate mutal entropy for all of these
require('entropy')

l <- length(sp.n)
int.sp.mi.dr <- vector(mode='list', length=(l^2 - l)/2)

k <- 0
for(i in 2:length(sp.n) - 1){
    for(j in (i+1):length(sp.n)){
        k <- k + 1
        mi <- mutual.info(intron.orth.t[[1]]$l[,i],
                          intron.orth.t[[1]]$l[,j],
                          numBins=20)
        sp <- colnames(intron.orth.t[[1]]$l)[c(i,j)]
        int.sp.mi.dr[[k]] <- c('i'=i, 'j'=j, 'sp1'=sp[1],
                               'sp2'=sp[2], mi)
    }
}


int.sp.mi.dr.m <- matrix(0, nrow=length(sp.n), ncol=length(sp.n))
int.sp.mi.dr.n <- matrix(0, nrow=length(sp.n), ncol=length(sp.n))
for(i in 1:length(int.sp.mi.dr)){
    x <-  int.sp.mi.dr[[i]]
    int.sp.mi.dr.m[ x$i, x$j ] = x$mi
    int.sp.mi.dr.m[ x$j, x$i ] = x$mi
    int.sp.mi.dr.n[ x$i, x$j ] = x$n
    int.sp.mi.dr.n[ x$j, x$i ] = x$n
}

rownames(int.sp.mi.dr.m) <- colnames( intron.orth.t[[1]][[1]] )
colnames(int.sp.mi.dr.m) <- colnames( intron.orth.t[[1]][[1]] )

saveRDS( int.sp.mi.dr.m, "int_sp_mi_dr_m.rds" )

y <- matrix( 1:length(sp.n), nrow=length(sp.n), ncol=length(sp.n), byrow=TRUE )
x <- matrix( 1:length(sp.n), nrow=length(sp.n), ncol=length(sp.n), byrow=FALSE )

par(mfrow=c(1,1))
cm <- 0.3  ## something to moderate the color
par(mar=c(1.1,1.1,1.1,1.1))
plot.new()
plot.window(xlim=c(-10,length(sp.n)), ylim=c(1, 10+length(sp.n)))
rect( x, y, x+1, y+1, col=hsvScale(int.sp.mi.dr.m, val=(cm + int.sp.mi.dr.m)/(cm + max(int.sp.mi.dr.m))  ), border='grey', lwd=0.25 )
text( 0.95, y[1,] + 0.5, sub("_[^_]+$", "", rownames(int.sp.mi.dr.m) ), adj=c(1,0.5), cex=0.8)
text( x[,1]+0.5, length(sp.n) + 1.05, sub("_[^_]+$", "", colnames(int.sp.mi.dr.m)), adj=c(0,0.5), srt=90, cex=1)

## there are some weird numbers here. In particular hucho hucho doesn't look much like a fish..
## hucho hucho vs esox lucius (pike)
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot( log2(intron.orth.t[[1]]$l[,'hucho hucho']), log2(intron.orth.t[[1]]$l[,'esox lucius']) )

## that is nicely correlated. No bimodality though
plot( log2(intron.orth.t[[1]]$l[,'hucho hucho']), log2(intron.orth.t[[1]]$l[,'danio rerio']) )
## not well correlated
plot( log2(intron.orth.t[[1]]$l[,'hucho hucho']), log2(intron.orth.t[[1]]$l[,'parambassis ranga']) )

## some other examples:
plot( log2(intron.orth.t[[1]]$l[,'gasterosteus aculeatus']), log2(intron.orth.t[[1]]$l[,'cottoperca gobio']) )

h.sp <- names( sort( int.sp.mi.dr.m['homo sapiens',], decreasing=TRUE))
for(i in 1:length(h.sp)){
    plot( log2(intron.orth.t[[1]]$l[,'homo sapiens']), log2(intron.orth.t[[1]]$l[,h.sp[i]]), xlab='homo sapiens', ylab=h.sp[i] )
    inpt <- readline("next: ")
    if(inpt == 'q')
        break
}

h.sp <- names( sort( int.sp.mi.dr.m['gasterosteus aculeatus',], decreasing=TRUE))
for(i in 1:length(h.sp)){
    plot( log2(intron.orth.t[[1]]$l[,'gasterosteus aculeatus']), log2(intron.orth.t[[1]]$l[,h.sp[i]]), xlab='stickleback', ylab=h.sp[i] )
    inpt <- readline(h.sp[i])
    if(inpt == 'q')
        break
}

tf.sp <- names( sort( int.sp.mi.dr.m['takifugu rubripes',], decreasing=TRUE))
for(i in 1:length(h.sp)){
    plot( log2(intron.orth.t[[1]]$l[,'takifugu rubripes']), log2(intron.orth.t[[1]]$l[,tf.sp[i]]), xlab='takifugue rubripes', ylab=tf.sp[i],
         cex=0.5, xlim=c(4, max(log2(intron.orth.t[[1]]$l[,'takifugu rubripes']), na.rm=TRUE)), ylim=c(4, max(log2(intron.orth.t[[1]]$l[,tf.sp[i]]), na.rm=TRUE)) )
    inpt <- readline(tf.sp[i])
    if(inpt == 'q')
        break
}


fam.db.count <- tapply( 1:nrow(gene.stats), gene.stats$family, function(i){
    tbl <- table( c(dbs, gene.stats$db[i]) ) - 1
    sum( tbl > 0 & tbl < 3 )
})

## we do 1 or 2 entries as using one gives use very little. This could
## be due to the presence of a salmonid that has about twice as many genes
## in its genome as the other species.

## This is for a single family member from each species
## sort( table( fam.db.count ) )
## fam.db.count
##  78 112 116 117 122 126 172 123 124 127 125 171 128 129 170 169 130 168 132 141 
##   1   1   1   1   1   1   1   2   3   3   4   5   6   9  16  26  44  66  71  78 
## 133 139 131 135 134 136 167 137 144 142 146 140 145 143 147 138 150 148 149 166 
##  81  84  85  86  87  92  92  93  95 102 102 110 115 119 133 136 148 152 155 160 
## 151 165 152 154 153 156 155 163 159 164 161 157 158 160 162 
## 164 169 173 197 198 234 240 243 248 252 266 271 287 288 317 

## there is exactly one family that is represented exactly once across the
## full data set.. 

sort(table(fam.db.count))
## fam.db.count
## 131 132 133 134 135 136 137 138 140 144 139 172 142 146 141 145 148 143 147 151 
##   1   5  10  18  21  29  36  40  42  42  43  47  50  51  56  61  64  65  66  69 
## 153 149 150 152 154 155 157 156 158 171 159 160 161 162 170 163 164 169 165 166 
##  77  78  81  87 110 115 129 131 138 151 191 204 207 268 274 332 376 428 439 468 
## 167 168 
## 507 507 


sum(fam.db.count > 150 ) ## 5255
sum(fam.db.count > 164 ) ## 2821
sum(fam.db.count > 169 ) ## 472
sum(fam.db.count > 170 ) ## 198

(472 * (172^2 - 172) /2) / 1e6
## 6.9 million alignments... Is this really feasible?
## seems that I actually have 4 million alignments in ex.align, so
## this is actually quite doable.

((172^2 - 172) /2 * 472 * 10000) / 1e9
## We may need something like 70GB of memory..
## but much less if we do not keep the alignments themselves but only the relevant alignment statistics
## but we have to consider how we handle duplicates. Easiest is to simply use the first one rather than
## to to include a whole set of other information. That can be done by using match on the sequences.
exon.f2 <- paste("../family_members/orthologue_transcripts/exons_", names(fam.db.count)[ fam.db.count >= 170 ], sep="")
all(file.exists(exon.f2))  ## TRUE!
names(exon.f2) <- names(fam.db.count)[ fam.db.count >= 170 ]

## all of these alignments should consist of a list of lists of alignments statistics returned by the
## align.seqs.mt.

al.stats.tmp <- lapply( 2:length(sp.n) - 1, function(i){
    v <- vector(mode='list', length=(length(sp.n) - i))
    names(v) <- sp.n[ (i+1):length(sp.n) ]
    v
})
names(al.stats.tmp) <- sp.n[ 2:length(sp.n)-1 ]

ex.align.2 <- lapply(exon.f2, function(f){
    al.stats <- al.stats.tmp
    tmp <- read.exons(f)
    tmp.sp <- sub(" ", "_", sapply( tmp, function(x){ x$sp }))
    o <- match(sp.n, tmp.sp)
    print(f)
    for(i in 1:(length(o)-1) ){
        if(is.na(o[i]))
            next
        ## i1 and i2 are used to extract the sequences to be aligned
        ## i and (i+1):length(o) - i give the corresponding indices for the
        ## al.stats vector which we will return.. We will encode these in j1 and j2
        ## and then remove those which we do not have an index for
        i1 <- o[i]
        i2 <- o[(i+1):length(o)]
        j1 <- i
        j2 <- (i+1):length(o) - i
        ## remove entries where we do not have sequences
        j2 <- j2[ !is.na(i2) ]
        i2 <- i2[ !is.na(i2) ]
        if(length(i2) < 1)
            next
        seq1 <- tmp[[i1]]$s
        seq2 <- lapply( tmp[i2], function(x){ x$s })
        seq2 <- unlist(seq2)
        id1 <- tmp[[i1]]$id
        id2 <- sapply( tmp[i2], function(x){ x$id })
        ## There should now not be any NA sequences. We can call error if this is the case
        if( sum(is.na( seq2 )) > 0 ){
            warning(paste("NA sequences found for ", f))
            next
        }
        cat("  ", length(seq2))
        seq.al <- align.seqs.mt( seq1, seq2, al.offset=sub.matrix$offset, al.size=sub.matrix$size,
                                sub.matrix=sub.matrix$sm, gap=as.integer(c(-10, -1)), tgaps.free=TRUE,
                                sp.char='I', thread.n=20)
        for(k in 1:length(seq.al)){
            al.stats[[j1]][[ j2[k] ]] <- c(seq.al[[k]][c('stats', 'a.pos', 'b.pos', 'al.i')], 'id1'=id1, 'id2'=id2[k])
        }
    }
    cat("\n")
    al.stats
})

save.image()

saveRDS( ex.align.2, "ex_align_2.rds" )

## I didn't record the species which I used there. But never mind.
## We want to define a structure containing
## ex.align.2.stats[[sp1]][[sp2]][[ family ]] = stats..

ex.align.2.stats <- lapply( sp.n, function(a){
    l <- lapply(sp.n, function(b){
        v <- vector(mode='list', length=length(ex.align.2))
        names(v) <- names(ex.align.2)
        v
    })
    names(l) <- sp.n
    l
})
names(ex.align.2.stats) <- sp.n

## then we fill in this data set in a somewhat wasteful fashion..
for(i in 1:length(ex.align.2)){
    for(sp1 in names(ex.align.2[[i]])){
        for(sp2 in names(ex.align.2[[i]][[sp1]])){
            if(is.null(ex.align.2[[i]][[sp1]][[sp2]]))
                next
            stats <- align.mt.to.stats( ex.align.2[[i]][[sp1]][[sp2]])
            ex.align.2.stats[[ sp1]][[ sp2]][[i]] <- stats
            ex.align.2.stats[[ sp2]][[ sp1]][[i]] <- stats
        }
    }
}


ex.align.2.jc <- matrix(0, nrow=length(sp.n), ncol=length(sp.n))
dimnames(ex.align.2.jc) <- list(names(genome.sz), names(genome.sz))
ex.align.2.jcg <- ex.align.2.jc
ex.align.2.k2 <- ex.align.2.jc

for(sp1 in names(ex.align.2.stats)){
    for(sp2 in names(ex.align.2.stats[[sp1]])){
        if(i == j)
            next
        mgd <- merge.stats( ex.align.2.stats[[sp1]][[sp2]] )
        ex.align.2.jc[ sp1, sp2 ] <- jukes.cantor(mgd)
        ex.align.2.jcg[ sp1, sp2 ] <- jukes.cantor.indel(mgd)
        ex.align.2.k2[ sp1, sp2 ] <- kimura.two(mgd)
    }
}

saveRDS( ex.align.2.jc, "ex_align_2_jc.rds" )
saveRDS( ex.align.2.jcg, "ex_align_2_jcg.rds" )
saveRDS( ex.align.2.k2, "ex_align_2_k2.rds" )

y <- matrix( 1:length(sp.n), nrow=length(sp.n), ncol=length(sp.n), byrow=TRUE )
x <- matrix( 1:length(sp.n), nrow=length(sp.n), ncol=length(sp.n), byrow=FALSE )

par(mfrow=c(1,1))
cm <- 0.3  ## something to moderate the color
par(mar=c(1.1,1.1,1.1,1.1))
plot.new()
plot.window(xlim=c(-10,length(sp.n)), ylim=c(1, 10+length(sp.n)))
rect( x, y, x+1, y+1, col=hsvScale(ex.align.2.jc), border='grey', lwd=0.25 )
##rect( x, y, x+1, y+1, col=hsvScale(ex.align.2.jc, val=(cm + ex.align.2.jc)/(cm + max(ex.align.2.jc))  ), border='grey', lwd=0.25 )
text( 0.95, y[1,] + 0.5, sub("_[^_]+$", "", rownames(ex.align.2.jc) ), adj=c(1,0.5), cex=0.8)
text( x[,1]+0.5, length(sp.n) + 1.05, sub("_[^_]+$", "", colnames(ex.align.2.jc)), adj=c(0,0.5), srt=90, cex=1)

par(mar=c(1.1,1.1,1.1,1.1))
plot.new()
plot.window(xlim=c(-10,length(sp.n)), ylim=c(1, 10+length(sp.n)))
rect( x, y, x+1, y+1, col=hsvScale(ex.align.2.jcg), border='grey', lwd=0.25 )
##rect( x, y, x+1, y+1, col=hsvScale(ex.align.2.jc, val=(cm + ex.align.2.jc)/(cm + max(ex.align.2.jc))  ), border='grey', lwd=0.25 )
text( 0.95, y[1,] + 0.5, sub("_[^_]+$", "", rownames(ex.align.2.jc) ), adj=c(1,0.5), cex=0.8)
text( x[,1]+0.5, length(sp.n) + 1.05, sub("_[^_]+$", "", colnames(ex.align.2.jc)), adj=c(0,0.5), srt=90, cex=1)

par(mar=c(1.1,1.1,1.1,1.1))
plot.new()
plot.window(xlim=c(-10,length(sp.n)), ylim=c(1, 10+length(sp.n)))
rect( x, y, x+1, y+1, col=hsvScale(ex.align.2.k2), border='grey', lwd=0.25 )
##rect( x, y, x+1, y+1, col=hsvScale(ex.align.2.jc, val=(cm + ex.align.2.jc)/(cm + max(ex.align.2.jc))  ), border='grey', lwd=0.25 )
text( 0.95, y[1,] + 0.5, sub("_[^_]+$", "", rownames(ex.align.2.jc) ), adj=c(1,0.5), cex=0.8)
text( x[,1]+0.5, length(sp.n) + 1.05, sub("_[^_]+$", "", colnames(ex.align.2.jc)), adj=c(0,0.5), srt=90, cex=1)


## plot genetic distance vs mutual information for the full set.. 
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot( ex.align.2.k2, int.sp.mi.dr.m, xlab='Kimura 2 factor', ylab='intron size mutual information' )

## to make sense of that we need to have the species lineages..
sp.lineages <- readRDS("species_lineages.rds")
## we may need to modify some of the lineage names here

substr( names(sp.lineages), 1, 1 ) <- tolower(substr( names(sp.lineages), 1, 1))

## unfortunately, I didn't manage to be consistent, so
## colnames in ex.align.2.k2, etc, are genus_species
## colnames in int.sp.mi.dr.m are 'genes species'
## names in sp.lineages are 'genes species'
## but since the order is the same, this should be ok

all(colnames( ex.align.2.k2 ) == sub(" ", "_", colnames(int.sp.mi.dr.m) ))  ## TRUE

which(!colnames(int.sp.mi.dr.m) %in% names(sp.lineages))

colnames(int.sp.mi.dr.m)[ !colnames(int.sp.mi.dr.m) %in% names(sp.lineages) ]
## [1] "canis familiaris" "mus pahari"       "mus caroli"       "mus spretus"     
## [5] "cebus capucinus"  "gorilla gorilla" 
## 'canis familiaris' -> canis lupus
## 'mus pahari' -> mus pahari strain PAHARI_EIJ
## 'mus caroli' -> mus caroli strain CAROLI_EIJ
## 'mus spretus' -> mus spretus strain SPRET/EiJ
## 'cebus capucinus' -> cebus capucinus imitator
## 'gorilla gorilla' -> 'gorilla gorilla gorilla'
##
## The easiest thing to do here is to change the names in the sp.lineages

names(sp.lineages)[ !(names(sp.lineages) %in% colnames(int.sp.mi.dr.m)) ]
## [1] "canis lupus"                  "mus pahari strain PAHARI_EIJ"
## [3] "mus caroli strain CAROLI_EIJ" "mus spretus strain SPRET/EiJ"
## [5] "cebus capucinus imitator"     "gorilla gorilla gorilla"     

names(sp.lineages)[ !(names(sp.lineages) %in% colnames(int.sp.mi.dr.m)) ] <- c(
    'canis familiaris', 'mus pahari', 'mus caroli', 'mus spretus', 'cebus capucinus', 'gorilla gorilla')

colnames(int.sp.mi.dr.m)[ !colnames(int.sp.mi.dr.m) %in% names(sp.lineages) ] ## empty vector

## this means that we can now define vectors of different types..
mammal.b <- sapply( colnames(int.sp.mi.dr.m), function(sp){ 'Mammalia' %in% sp.lineages[[sp]][,'node_name'] }) ## 80
aves.b <- sapply( colnames(int.sp.mi.dr.m), function(sp){ 'Aves' %in% sp.lineages[[sp]][,'node_name'] }) ## 22
teleost.b <- sapply( colnames(int.sp.mi.dr.m), function(sp){ 'Osteoglossocephalai' %in% sp.lineages[[sp]][,'node_name'] }) ## 22
sauria.b <-  sapply( colnames(int.sp.mi.dr.m), function(sp){ 'Sauria' %in% sp.lineages[[sp]][,'node_name'] }) ## 22

## set up colours so that:
## red => mammal to mammal
## green => sauria to sauria
## blue => teleost to teleost

square.b <- function(x){
    sapply(1:length(x), function(i){ x[i] & x })
}

image(square.b(mammal.b))
image(square.b(aves.b))
image(square.b(teleost.b))

## seems ok..
red <- as.numeric(square.b( mammal.b ))
green <- as.numeric(square.b( sauria.b ))
blue <- as.numeric(square.b( teleost.b ))


## plot genetic distance vs mutual information for the full set.. 
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot( ex.align.2.k2, int.sp.mi.dr.m, xlab='Kimura 2 factor', ylab='intron size mutual information',
     col=rgb(red, green * 0.7, blue, 0.5), pch=19)

genus <- sub("(\\S+)\\s+\\S+", "\\1",  colnames(int.sp.mi.dr.m))
labels <- paste( matrix(genus, nrow=length(genus), ncol=length(genus)),
                 matrix(genus, nrow=length(genus), ncol=length(genus), byrow=TRUE),
                 sep="\n" )

identify( ex.align.2.k2, int.sp.mi.dr.m, labels=labels )

## it looks like there is a difference between eurasian and other mammal types..
## but I need a better lineage data to do anything systematic about that. Hence
## I will export the taxonomy trees that I have used before..

## This file is from my workstation, at:
## /home/lmj/taxonomy/NCBI/perl
tax.tree <- read.table('bold_eukaryote_coverage_2', header=FALSE, sep="\t",
                       stringsAsFactors=FALSE, quote="", comment.char="" )
colnames(tax.tree) <- c('id', 'p.id', 'bold.n', 'bold.seq.n', 'rank', 'rank.l', 'level', 'name', 'x', 'y', 'tax.n', 'd.n')
rownames(tax.tree) <- tax.tree[,'id']  ## these are numbers, so need to be used carefully.

## we need a function to extract a lineage from this tree..
## this is incredibly slow. The rownames lookup is crazily slow
## seems like it does a linear scan across it.
## Much faster to do which(... ) than to use the rownames attribute
## That is weird as hell
extract.lineage <- function(tree, taxon, capitalise=TRUE){
    if( capitalise )
        substr(taxon, 1, 1) <- toupper( substr(taxon, 1, 1))
    i <- which( tree$name == taxon )
    if(length(i) == 0)
        return(i)
    if(length(i) > 1){
        warning("more than one entry for: ", taxon )
        i <- i[1]
    }
    ii <- i
    ## we might be able to use the level to set the size of the matrix so that
    ## we don't need to use rbind here.. 
    while(TRUE){
        if(tree[i, 'level'] == 0)
            break
        i <- which( tree$id == tree[i, 'p.id'] )
        ii <- c(ii, i)
    }
    tree[ii, ]
}

sum(sp.names %in% tax.tree$name ) ## 171
## so we are missing one. .Which one is it I wonder
sp.names[ !(sp.names %in% tax.tree$name) ]
## Ahh, the old 'Canis familiaris' which should be Canis lupus familiaris..
## The simplest, but not good option here is to change the tax.tree name

tax.tree[ tax.tree$name == 'Canis lupus familiaris', 'name'] <- 'Canis familiaris'
sum(sp.names %in% tax.tree$name ) ## 172

## then we can do this for the 
ncbi.lineages <- lapply( colnames(int.sp.mi.dr.m), function(sp){
    print(paste("getting lineage for", sp))
    extract.lineage(tax.tree, sp)
})
names(ncbi.lineages) <- colnames(int.sp.mi.dr.m)
## that now seems a bit more resonable.. Can I get some information for some of the
## the genuss...
saveRDS(ncbi.lineages, "ncbi_lineages.rds")

metatheria.b <- sapply( ncbi.lineages, function(x){ 'Metatheria' %in% x$name } ) ## 5
eutheria.b <- sapply( ncbi.lineages, function(x){ 'Eutheria' %in% x$name } ) ## 74
xenarthra.b <- sapply( ncbi.lineages, function(x){ 'Xenarthra' %in% x$name } ) ## only 2,
boreoutheria.b <- sapply( ncbi.lineages, function(x){ 'Boreoeutheria' %in% x$name } ) ## 69

image(square.b(metatheria.b))
image(square.b(boreoutheria.b))
image(square.b(xenarthra.b))

red <- as.numeric(square.b( eutheria.b ))
red <- as.numeric(square.b( boreoutheria.b ))
green <- as.numeric(square.b( metatheria.b ))
blue <- as.numeric(square.b( xenarthra.b ))

par(mar=c(5.1, 4.1, 4.1, 2.1))
plot( ex.align.2.k2, int.sp.mi.dr.m, xlab='Kimura 2 factor', ylab='intron size mutual information',
     col=rgb(red, green * 0.7, blue, 0.5), pch=19, cex=1.25)

identify( ex.align.2.k2, int.sp.mi.dr.m, labels=labels )

###### Can we do a PCA with all of the intron sizes:
sum( apply(intron.orth.t[[1]]$l, 1, function(x){ all(!is.na(x)) })) ## 45

## there are only 45 introns that I can find in all species. That's obviously silly...
par(mar=c(14.1, 4.1, 4.1, 2.1))
barplot( apply( intron.orth.t[[1]]$l, 2, function(x){ sum(!is.na(x)) }), las=2 )

sort( apply( intron.orth.t[[1]]$l, 2, function(x){ sum(!is.na(x)) }), decreasing=TRUE )
## we can probably restrict ourselves to species with at least 40,000 introns

b <- apply( intron.orth.t[[1]]$l, 2, function(x){ sum(!is.na(x)) }) > 40000
## that is 167 / 172 species
sum( apply(intron.orth.t[[1]]$l[,b], 1, function(x){ all(!is.na(x)) })) ## 158

## not much of an improvement
b <- apply( intron.orth.t[[1]]$l, 2, function(x){ sum(!is.na(x)) }) > 45000
## that is 160 / 172 species
sum( apply(intron.orth.t[[1]]$l[,b], 1, function(x){ all(!is.na(x)) })) ## 268

b <- apply( intron.orth.t[[1]]$l, 2, function(x){ sum(!is.na(x)) }) > 50000
## that is 139 / 172 species
sum( apply(intron.orth.t[[1]]$l[,b], 1, function(x){ all(!is.na(x)) })) ## 667

## still not so good as that probably does not represent a random sampling..
## Let's set the mean to something else..


## replacing the missing values like this is not a good idea.. 
int.l <- log2( intron.orth.t[[1]]$l )

for(i in 1:nrow(int.l)){
    b <- !(is.na(int.l[i,])) & is.finite(int.l[i,])
    m <- mean( int.l[i,b] )
    if(is.na(m))
        break
    bb <- is.na(int.l[i, ]) & !(is.finite(int.l[i,]))
    int.l[i, bb] <- m
}
int.l.sd <- apply( int.l, 1, sd )


int.l.pca <- prcomp( t( int.l[int.l.sd > 0 & !is.na(int.l.sd), ] ), scale=TRUE )
plot(int.l.pca) ## almost all of the variance in the first dimension. This will be genome size I guess

plot(int.l.pca$x[1,], int.l.pca$x[,2], pch=19, col=hsvScale( log2(genome.sz) ), cex=2, xlab='dim 1', ylab='dim 2')
## ok, that's good we have a small numer of outliers. Removing those is kind of a good idea..
identify(int.l.pca$x[1,], int.l.pca$x[,2], labels=rownames(int.l.pca$x))

## ok, that is actually a bit weird. Variance in distance from the mean is a function of genome size.
## For a large genome size, everything lines up along the 0 on the x-axis (i.e. the center). I suspect that this
## might be related to the scaling. I have to think about the relevant equations.

## we should try without scaling..
int.l.pca.2 <- prcomp( t( int.l[int.l.sd > 0 & !is.na(int.l.sd), ] ), scale=FALSE )
plot(int.l.pca.2) ## again, almost all variance in a single dimension
plot(int.l.pca.2$x[1,], int.l.pca.2$x[,2], pch=19, col=hsvScale( log2(genome.sz) ), cex=2, xlab='dim 1', ylab='dim 2')
identify(int.l.pca.2$x[1,], int.l.pca.2$x[,2], labels=rownames(int.l.pca$x))

## that looks the same; the weird thing is that we see tetraodon and takifuge very far apart in dimenion 1,
## but we expect that they should be well correlated.

plot( int.l[,'takifugu rubripes'], int.l[,'tetraodon nigroviridis'], xlab='fugu', ylab='tetraodon' )
## ok, int.l is screwed up by something..

plot( log2( intron.orth.t[[1]]$l[,'takifugu rubripes']), log2(intron.orth.t[[1]]$l[,'tetraodon nigroviridis']),
     xlim=c(4, max(int.l[,'takifugu rubripes'], na.rm=TRUE)), ylim=c(4, max(int.l[,'tetraodon nigroviridis'], na.rm=TRUE)),
     cex=0.5)

plot( int.l[,'takifugu rubripes'], int.l[,'tetraodon nigroviridis'], xlab='takifugue rubripes', ylab='tetraodon nigroviridis',
     cex=0.5, xlim=c(4, max(int.l[,'takifugu rubripes'], na.rm=TRUE)), ylim=c(4, max(int.l[,'tetraodon nigroviridis'], na.rm=TRUE)) )


plot( log2( intron.orth.t[[1]]$l[,'danio rerio']), log2(intron.orth.t[[1]]$l[,'tetraodon nigroviridis']),
     xlim=c(4, max(int.l[,'danio rerio'], na.rm=TRUE)), ylim=c(4, max(int.l[,'tetraodon nigroviridis'], na.rm=TRUE)),
     cex=0.5)

plot( int.l[,'danio rerio'], int.l[,'tetraodon nigroviridis'], xlab='danio rerio', ylab='tetraodon nigroviridis',
     cex=0.5, xlim=c(4, max(int.l[,'danio rerio'], na.rm=TRUE)), ylim=c(4, max(int.l[,'tetraodon nigroviridis'], na.rm=TRUE)) )

## we can use a PCA function that tries to impute missing and weird data. Not sure how well it will work here, but we can try
## to see what happens.

install.packages("cellWise")

require(cellWise)

inttron.orth.mpca.1 <- MacroPCA( t(log(intron.orth.t[[1]]$l)), k=0 )
inttron.orth.mpca.2 <- MacroPCA( t(log(intron.orth.t[[1]]$l)), k=10 )

with( inttron.orth.mpca.2, {
    plot( scores[,1], scores[,2], col=hsvScale( log(genome.sz[ sub(" ", "_", rownames(scores)) ]) ), pch=19, cex=1.5 )
    identify( scores[,1], scores[,2], labels=rownames(scores) )
})


mpca.pars <- list(scale=FALSE)
intron.orth.mpca.3 <- MacroPCA( t(log(intron.orth.t[[1]]$l)), k=10, MacroPCApars=mpca.pars )

with( intron.orth.mpca.3, {
    plot( scores[,1], scores[,2], col=hsvScale( log(genome.sz[ sub(" ", "_", rownames(scores)) ]) ), pch=19, cex=1.5 )
    identify( scores[,1], scores[,2], labels=rownames(scores) )})
} )

## These plots do not give us very much extra. There is a general correlation with genome size, but also a correlation
## with taxonomy with closely related species clustered together. There are some curious outgroups like rodents and
## primates. Presumably these are outliers not because of biology, but because of the extent of annotation
## but I am not sure how to assess that question.


### We can look at the variance of individual introns; both for the full set and within teleost groups

var.par <- function(x){
    c('n'=sum(!is.na(x)), 'mean'=mean(x, na.rm=TRUE), 'sd'=sd(x, na.rm=TRUE),
      'min'=min(x, na.rm=TRUE), 'max'=max(x, na.rm=TRUE), quantile(x, probs=seq(0,1,0.1), na.rm=TRUE))
}


teleost.var <- t(apply( log2(intron.orth.t[[1]]$l[,teleost.b]), 1, var.par ))
mammal.var <-  t(apply( log2(intron.orth.t[[1]]$l[,mammal.b]), 1, var.par ))
sauria.var <-  t(apply( log2(intron.orth.t[[1]]$l[,sauria.b]), 1, var.par ))

## we are then interested in the distributions of the varios parameters..
hist(teleost.var[,'n'])
range(teleost.var[,'n'])  ## 1 -> 54
t.min.n <- 49
b1 <- teleost.var[,'n'] >= t.min.n
sum(b1) ## 38,815, not too bad.


range(teleost.var[b1,'sd'], na.rm=TRUE)
## [1]  0.06701767 3.07486888
## we have some NA values..
sum( is.na(teleost.var[,'sd']))
## 420
head(which( is.na(teleost.var[,'sd'])))
## [1]  194  338  339  746  815 1530

sum(is.na(teleost.var[,'sd']) & b1 )
## 0

hist( teleost.var[b1,'sd'] ) ## not too bad..
## peak is at about 1.25, meaning a standard deviation of something pretty small

plot( teleost.var[b1, c('mean','sd')] )
## The correlation is restricted to introns of size 2^7.5 and smaller
## longer means decrease in variance, which is nice to see

o1 <- order( teleost.var[,'sd'] )
o2 <- rev(o1)
o1.b <- teleost.var[ o1, 'n'] >= t.min.n
o2.b <- teleost.var[ o2, 'n'] >= t.min.n

head(teleost.var[ o1[o1.b], ])
head(teleost.var[ o2[o2.b], ])
## I want to look at the lengths of these

x <- sapply( sp.lineages[ colnames( intron.orth.t[[1]]$l ) ], function(x){ x[1,'left_index'] })
sp <- colnames(intron.orth.t[[1]]$l)

## Here I try to do the alignments of introns from selected species. This can fail when we have more
## than one transcript from a single gene since I only made use of gene ids; not transcript ids.
## When we have more than one transcript we can look at the expected intron sizes and select
## the proper one; that should work most of the time, but it would be better to fix the code
## Unfortunately that will take a couple of days to run due to it being incredibly slow to
## iterate through all of the alignments in ex.align. (I don't have the information there, but
## I do have the index of the sequences that were aligned. Not as good, but it can be used.
for(i in o1[o1.b]){
    plot(x, log2(intron.orth.t[[1]]$l[ i, ]), col=ifelse(teleost.b, 'red', 'black'), pch=19, cex=1.5 )
    inpt <- readline( "next: ")
    if(grepl('i', inpt))
        identify( x, log2(intron.orth.t[[1]]$l[ i, ]), labels=colnames(intron.orth.t[[1]]$l ), las=2 )
    if(grepl('a', inpt)){
        ## identify is broken; it order the return value by position, not by selection order
        seed.sp.i <- identify( x, log2(intron.orth.t[[1]]$l[ i, ]), labels=sp, col='red' )
        sp.i <- identify( x, log2(intron.orth.t[[1]]$l[ i, ]), labels=sp )
        tmp.s <- read.exons( intron.f[ intron.orth.fam[[1]][i] ] )
        ## then we want to align the correct sequences
        seed.id <- intron.orth.t[[1]]$id[i, seed.sp.i[1]]
        seed.int.i <- intron.orth.t[[1]]$i[i, seed.sp.i[1]]
        id <- intron.orth.t[[1]]$id[i, sp.i]
        int.i <- intron.orth.t[[1]]$i[i, sp.i]
        tmp.s.id <- sapply( tmp.s, function(x){ x$id })
        i1 <- which( tmp.s.id == seed.id )
##        i2 <- which( tmp.s.id %in% id )
        seed.seq <- tmp.s[[i1]]$e[ seed.int.i ]
        q.seq <- sapply( 1:length(id), function(j){
            i2 <-  which( tmp.s.id == id[j] )
            tmp.s[[i2]]$e[ int.i[j] ]
        })
        al <- align.seqs.mt( seed.seq, q.seq, al.offset=sub.matrix$offset, al.size=sub.matrix$size,
                            sub.matrix=sub.matrix$sm, gap=as.integer(c(-10, -1)), tgaps.free=TRUE, sp.char='I', thread.n=20 )
        for(j in 1:length(al)){
            print(paste(sp[seed.sp.i[1]], "vs", sp[sp.i[j]], " :", seed.id, "vs", id[j]))
            align.print( al[[j]]$seq, 80 )
        }
        inpt <- readline("continue or (q)uit? ")
    }
    if(grepl('q', inpt))
        break
}
## introns that are long and have low variance within the teleosts are very much enriched for
## developmental regulators. The question is if that is also true simply for long introns,
## from a given species. Might give the same sort of data...
## I would also have to consider the fact that the intron number affects the likelihood of
## inclusion in any selection based on individual introns... 

### Some of the plots below give a very bimodal pattern. That suggests that we may have an issue with misidentifcation of
### of orthologous introns. That means that we should try to examine the evidence from the original alignments
### We can do this similarly for a set of species that we choose. But then we will want to plot them on a different
### device. So We need to have another plotting device... 
x <- sapply( sp.lineages[ colnames( intron.orth.t[[1]]$l ) ], function(x){ x[1,'left_index'] })
sp <- colnames(intron.orth.t[[1]]$l)
al.cols <- c( rep('grey', 4), 'yellow', 'white', 'red')
names(al.cols) <- c('A', 'C', 'G', 'T', 'N', '-', 'I')

plt.dev <- 3
al.dev <- 4
for(i in o2[o2.b]){
    plot(x, log2(intron.orth.t[[1]]$l[ i, ]), col=ifelse(teleost.b, 'red', 'black'), pch=19, cex=1.5 )
    inpt <- readline( "next: ")
    if(grepl('i', inpt))
        identify( x, log2(intron.orth.t[[1]]$l[ i, ]), labels=colnames(intron.orth.t[[1]]$l ), las=2 )
    if(grepl('a', inpt)){
        ## identify is broken; it order the return value by position, not by selection order
        seed.sp.i <- identify( x, log2(intron.orth.t[[1]]$l[ i, ]), labels=sp, col='red' )
        sp.i <- identify( x, log2(intron.orth.t[[1]]$l[ i, ]), labels=sp )
        tmp.s <- read.exons( intron.f[ intron.orth.fam[[1]][i] ] )
        ## then we want to align the correct sequences
        seed.id <- intron.orth.t[[1]]$id[i, seed.sp.i[1]]
        seed.int.i <- intron.orth.t[[1]]$i[i, seed.sp.i[1]]
        id <- intron.orth.t[[1]]$id[i, sp.i]
        int.i <- intron.orth.t[[1]]$i[i, sp.i]
        tmp.s.id <- sapply( tmp.s, function(x){ x$id })
        i1 <- which( tmp.s.id == seed.id )
##        i2 <- which( tmp.s.id %in% id )
        seed.seq <- tmp.s[[i1]]$e[ seed.int.i ]
        q.seq <- sapply( 1:length(id), function(j){
            i2 <-  which( tmp.s.id == id[j] )
            tmp.s[[i2]]$e[ int.i[j] ]
        })
        al <- align.seqs.mt( seed.seq, q.seq, al.offset=sub.matrix$offset, al.size=sub.matrix$size,
                            sub.matrix=sub.matrix$sm, gap=as.integer(c(-10, -1)), tgaps.free=TRUE, sp.char='I', thread.n=20 )
        for(j in 1:length(al)){
            print(paste(sp[seed.sp.i[1]], "vs", sp[sp.i[j]], " :", seed.id, "vs", id[j]))
            align.print( al[[j]]$seq, 80 )
        }
        inpt <- readline("continue or (q)uit? ")
    }
    if(grepl('x', inpt)){
        ## then plot the original exon alignments for the selected species. Stick to the first seed species
        sp.i <- identify( x, log2(intron.orth.t[[1]]$l[ i, ]), labels=sp )
        sp.id <- sp[sp.i]
        fam.id <- intron.orth.fam[[1]][i]
        ex <- ex.align[[fam.id]][[1]]  ## 
        sp.id <- union(sp.id, ex$sp.1)
        ## we want to draw length(sp.i) alignments..
        ex.i <- which( ex.align[[fam.id]][[1]]$sp.2 %in% sp.id )
        al.lengths <- sapply( ex$al[ex.i], function(x){ nchar( x$seq[1] )} )
        dev.set( al.dev )
        plot.new()
        h <- 5 * length(ex.i)
        plot.window( xlim=c(-0.1*max(al.lengths), 1.1*max(al.lengths)), ylim=c(-1, h + 1) )
        y <- h
        for(j in 1:length(sp.id)){
            for(k in which( ex$sp.2 == sp.id[j])){
                draw.aligns( ex$al[[ k ]], y, y-2, 1, cols=al.cols, sp.a=ex$sp.1, sp.b=sp.id[j] )
                text(0, y-0.5, ex$al[[k]]$stats['score'], adj=c(1.2,0.5), cex=1.5 )
                y <- y - 5
            }
        }
        dev.set(plt.dev)
        inpt <- readline("continue or (q)uit? ")
    }
    if(grepl('q', inpt))
        break
}

## There almost looks like there is a relationship between teleost variance and non teleost variance
## let us chack this..
non.teleost.var <- t(apply( log2(intron.orth.t[[1]]$l[,!teleost.b]), 1, var.par ))
hist(non.teleost.var[,'n'])
range(non.teleost.var[,'n'])  ## 9 -> 118

all.var <- t(apply( log2(intron.orth.t[[1]]$l), 1, var.par ))

nt.min.n <- 103
b2 <- non.teleost.var[,'n'] >= nt.min.n
sum(b2) ## 40,180

sum(b1 & b2)  ## 29,967

## is the variance related?
plot( teleost.var[b1 & b2, 'sd'], non.teleost.var[b1 & b2, 'sd'] ) ## no relationships!
abline(0, 1, col='red')

par(mfrow=c(1,2))
hist( teleost.var[b1 & b2, 'sd'] )
hist( non.teleost.var[b1 & b2, 'sd'] )
## higher variance in teleosts; but lots of closely related mammals might explain
## this (though inclusion of birds and reptiles ought to make this moot).

par(mfrow=c(1,1))
plot( teleost.var[b1 & b2, 'mean'], non.teleost.var[b1 & b2, 'mean'], col=hsvScale( teleost.var[b1 & b2, 'sd'], alpha=0.5), pch=19 )
## and that has the typical plot that we see between teleosts and mammal
## species; A weak relationship for longer exons.
abline(0, 1, col='red', lwd=2)

plot( teleost.var[, c('mean', 'sd')] )
abline(h=2.4, col='red')
hist( teleost.var[ teleost.var[,'sd'] < 1, 'mean'] )
hist( teleost.var[ teleost.var[,'sd'] >= 1, 'mean'] )

hist( non.teleost.var[ teleost.var[,'sd'] < 0.65, 'mean'] )
hist( non.teleost.var[ teleost.var[,'sd'] >= 1.5, 'mean'] )


    
### I can infer the ancestral intron sizes from the set of sizes that I have. But for this I want a
### phylogenetic tree. I can use ape for this

require('ape')

## ex.align.2.k2 contains distances based on exon alignments

ex.align.2.k2.nj <- nj( ex.align.2.k2 )
plot(ex.align.2.k2.nj, cex=0.8)
plot( root( ex.align.2.k2.nj, outgroup='eptatretus_burgeri'), cex=0.8)


## this is true..
all(colnames(intron.orth.t[[1]][[1]]) == sub("_", " ", colnames( ex.align.2.k2 )))

## so we can then make a table of log2 intron sizes that we can encode using functions in
dyn.load("~/R/R_max_parsimony/src/max_parsimony.so")
source("~/R/R_max_parsimony/R/functions.R")
## this might have some functions with the same name as in the exon_aligner functions..

int.s.l <- log2( 1 + intron.orth.t[[1]]$l )
range(int.s.l, na.rm=TRUE)
## [1]  0.00000 20.73986

saveRDS(int.s.l, "int_s_l.rds")

int.s.l.all.h <- hist( int.s.l, breaks=50 )

## lets plot all of these
int.s.l.h <- lapply( colnames(int.s.l), function(sn){
    hist(int.s.l[,sn], breaks=int.s.l.all.h$breaks, plot=FALSE)
})

names(int.s.l.h) <- colnames(int.s.l)


## We want the tree to be rooted on epatretus burgeri or petromyozon marinus
## which are both jawless vertebrates.

ex.align.2.k2.nj <- root( ex.align.2.k2.nj, outgroup='eptatretus_burgeri')

## save for posterity, so I can use in other sessions
saveRDS(ex.align.2.k2.nj, "ex_align_2_k2_nj.rds")

dev.set(2)
plot(ex.align.2.k2.nj, cex=0.8)
usr <- par("usr")


tree.lines <- edge.lines( ex.align.2.k2.nj )

dev.set(3)
plot.new()
plot.window(xlim=usr[1:2], ylim=usr[3:4], xaxs='i', yaxs='i')
with(tree.lines, segments(x[,1], y, x[,2], y))
with(tree.lines, segments(v[,1], v[,2], v[,1], v[,3]))
with(tree.lines, text( x[,2], y, labels=nodes, pos=4, cex=0.75))
with(tree.lines, {
    b <- !(nodes %in% ex.align.2.k2.nj$edge[,1])
    text( x[b,2] + 0.002, y[b], ex.align.2.k2.nj$tip.label[ nodes[b] ], cex=0.8, pos=4 )
    })

total.nodes <- length(unique(as.numeric(ex.align.2.k2.nj$edge)))
leaf.n <- length( with(ex.align.2.k2.nj, setdiff(edge[,2], edge[,1])) )

## we can use this to make a species order
sp.y <- with(tree.lines, {
    b <- !(nodes %in% ex.align.2.k2.nj$edge[,1])
    lab <- ex.align.2.k2.nj$tip.label[ nodes[b] ]
    y.pos <- y[b]
    names(y.pos) <- sub("_", " ", lab)
    y.pos
})

pdf("species_intron_size_distributions.pdf", height=8.27, width=11.69,
    title='Species intron size distributions')
par(mfrow=c(4,6))
par(mar=c(2.1, 4.1, 4.1, 2.1))
invisible(lapply( names(sort(sp.y, decreasing=TRUE)), function(sn){
    plot( int.s.l.h[[sn]], main=sn, xlab='', cex.main=0.75 )
    abline(v=8, col='red', lwd=1)
}))
dev.off()

### We can fit all of these to a normal mixture model using the mixtools
### library. There are probably other tools available, but this is the first
### one I came across; described here:
### https://www.r-bloggers.com/fitting-mixture-distributions-with-the-r-package-mixtools/
###

install.packages('mixtools')
require(mixtools)
## nicely up to date!
## to fit binormal data:

tmp <- normalmixEM( int.s.l[!is.na(int.s.l[,'poecilia reticulata']), 'poecilia reticulata'] )
plot( tmp, which=1 )

plot.mixEM( tmp, which=2 )
lines(density(tmp$x), lty=2, lwd=2)
plotMix( tmp )

norm.models <- vector(mode='list', length=ncol(intron.orth.t[[1]]$l) )
names(norm.models) <- colnames(intron.orth.t[[1]]$l)
for(i in 1:length(norm.models)){
    x <- log2( intron.orth.t[[1]]$l[,i] )
    x <- x[!is.na(x)]
    x <- x[is.finite(x)]
    norm.models[[i]] <- normalmixEM(x)
}

pdf("species_intron_size_models.pdf", height=8.27, width=11.69,
    title='Species intron size distributions')
par(mfrow=c(4,5))
par(mar=c(2.1, 4.1, 4.1, 2.1))
for(sn in names(sort(sp.y, decreasing=TRUE))){
    plotMix( norm.models[[sn]], main=sn, xlab='log2 intron length', legpos=NULL, cex.main=0.75,
            obs.lwd=4, obs.lty=1, obs.col=rgb(0.5, 0, 0, 1), mod.col=rgb(0, 0, 0.5, 0.5), mod.lwd=3)
    abline(v=8, col='red', lwd=1)
}
dev.off()

## we then need to make a substitution matrix
al.size <- as.integer( 1 + max(int.s.l, na.rm=TRUE) - min(int.s.l, na.rm=TRUE))
al.offset <- 64L  ## starts at @
sub.matrix <- make.sub.matrix( al.size )
## first row and column indicate missing values, and should have maximum penalties
missing.penalty <- as.integer(max(sub.matrix) + 1)
sub.matrix[1,] <- missing.penalty
sub.matrix[,1] <- missing.penalty

int.s.enc <- encode.dist( int.s.l, offset=al.offset )
## each one has a total of 63608 intron sizes which are considered (though
## with some missing values to be expected.

## then we can infer the ancestral states.
nj.edge <- matrix( as.integer(ex.align.2.k2.nj$edge), nrow=nrow(ex.align.2.k2.nj$edge))

int.s.inf <- sankoff( nj.edge, c(total.nodes, leaf.n), sub.matrix, c(al.offset, al.size),
                      int.s.enc )

saveRDS( int.s.inf, "int_s_inf.rds")
## and it seems that we got something from that..
## we can look at the inferred states (for leaf nodes these are just the same as that
## given.

dev.set(4)

## for takifugu, tetraodon, their ancestor and the ancestor of everything
par(mfrow=c(2,2))
hist( int.s.inf[[1]]$state[,1], main=names(int.s.enc)[1] )
hist( int.s.inf[[2]]$state[,1], main=names(int.s.enc)[2] )
hist( int.s.inf[[326]]$state[,1], main=paste(names(int.s.enc)[1:2], collapse=", ") )
hist( int.s.inf[[173]]$state[,1], main="common ancestor" )

## lets look at some mammals..
hist( int.s.inf[[154]]$state[,1], main=names(int.s.enc)[154] )
hist( int.s.inf[[155]]$state[,1], main=names(int.s.enc)[155] )
hist( int.s.inf[[186]]$state[,1], main=paste(names(int.s.enc)[c(154,155)], collapse=", ") )
hist( int.s.inf[[215]]$state[,1], main="mammals" )

hist( int.s.inf[[173]]$state[,1], main="common ancestor" )

## We can define the set of differences between a child and parent
## where we ignore missing values..

delta.abs <- sapply(tree.lines$nodes, function(i){
    p.i <- int.s.inf[[i]]$edge[ int.s.inf[[i]]$edge[,2] == 0, 1]
##    p.i <- which(!(int.s.inf[[i]]$edge[,2]))
    if(!length(p.i))
        return(0)
##    p.i <- int.s.inf[[i]]$edge[p.i,1]
    b <- int.s.inf[[i]]$state != 0 & int.s.inf[[p.i]]$state != 0
    mean( abs(int.s.inf[[i]]$state[b] - int.s.inf[[p.i]]$state[b]) )
})

delta <- sapply(tree.lines$nodes, function(i){
    p.i <- int.s.inf[[i]]$edge[ int.s.inf[[i]]$edge[,2] == 0, 1 ]
    if(!length(p.i))
        return(0)
    b <- int.s.inf[[i]]$state != 0 & int.s.inf[[p.i]]$state != 0
    mean( int.s.inf[[i]]$state[b] - int.s.inf[[p.i]]$state[b] )
})


dev.set(3)
plot.new()
plot.window(xlim=usr[1:2], ylim=usr[3:4], xaxs='i', yaxs='i')
with(tree.lines, segments(x[,1], y, x[,2], y, col=hsvScale(delta.abs), lwd=2))
with(tree.lines, segments(v[,1], v[,2], v[,1], v[,3]))
with(tree.lines, text( x[,2], y, labels=nodes, pos=4, cex=0.75))
with(tree.lines, {
    b <- !(nodes %in% ex.align.2.k2.nj$edge[,1])
    text( x[b,2] + 0.002, y[b], ex.align.2.k2.nj$tip.label[ nodes[b] ], cex=0.6, pos=4 )
    })

blue.to.red <- function(x, conv=sqrt){
    rgb( ifelse(x > 0, conv( (x-0)/max(x) ), 0),
         0,
         ifelse(x < 0, conv( (x-0)/min(x)), 0) )
}

#dev.set(3)
plot.new()
plot.window(xlim=usr[1:2], ylim=usr[3:4], xaxs='i', yaxs='i')
with(tree.lines, segments(x[,1], y, x[,2], y, col=blue.to.red(delta), lwd=2))
with(tree.lines, segments(v[,1], v[,2], v[,1], v[,3]))
with(tree.lines, text( x[,2], y, labels=nodes, pos=4, cex=0.75))
with(tree.lines, {
    b <- !(nodes %in% ex.align.2.k2.nj$edge[,1])
    text( x[b,2] + 0.002, y[b], ex.align.2.k2.nj$tip.label[ nodes[b] ], cex=0.8, pos=4 )
    })

m.i <- match(1:length(int.s.inf), tree.lines$nodes )

## there are a lot of ways in which we can look at this data.
## One way would be to plot the accumulated delta for all lineages;
## One would hope that would separate the groups in some reasonable manner..
## Also nice would be to plot phylogenetic edge length vs delta or abs(delta)
## to see what outliers we have.

## The internal nodes with the biggest reduction in size is 266, which is the
## teleostei common node. This is followed by the node that joins up
## the puffer fishes and some other teleost internal nodes. I believe that this
## probably argues for a general decrease in intron size in teleosts.

### traverse the leaf nodes to get the distance. Use both the tree (for
### phylogenetic differences) and a summarised change in intron size
### We have a total of leaf.n leaf nodes..

leaf.int.d <- lapply( 1:leaf.n, function(i){
    pg.dist <- c()
    int.d <- c()
    node.i <- c()
    ascend.tree <- function(j){
        b <- ex.align.2.k2.nj$edge[,2] == j
        if(sum(b) == 1){
            pg.dist <<- c(ex.align.2.k2.nj$edge.length[b], pg.dist)
            p.j <- ex.align.2.k2.nj$edge[b,1]
            b2 <- (int.s.inf[[j]]$state[,1] != 0 & int.s.inf[[p.j]]$state[,1] != 0)
            int.d <<- c(mean(int.s.inf[[j]]$state[b2,2]), int.d)
            node.i <<- c(j, node.i)
            ascend.tree( p.j )
        }
    }
    ascend.tree(i)
    cbind('pg'=cumsum(pg.dist),
          'int.d'=cumsum(int.d),
          'node'=node.i)
})
names(leaf.int.d) <- colnames( int.s.l )

saveRDS(leaf.int.d, "leaf_int_d.rds")

## collect the states of the ancestors of a given leaf
collect.ancestor.states <- function( leaf ){
    node.i <- leaf
    if(!is.numeric(leaf))
        node.i <- which(colnames(int.s.l) == leaf)
    if(!is.numeric(node.i))
        stop("Unable to get a numeric node identifier from: ", leaf)
    m <- matrix(nrow=nrow( int.s.inf[[node.i]]$state), ncol=0)
    ascend.tree <- function(j){
        b <- ex.align.2.k2.nj$edge[,2] == j
        if(sum(b) == 1){
            m <<- cbind( int.s.inf[[j]]$state[,1], m )
            p.j <- ex.align.2.k2.nj$edge[b,1]
            ascend.tree( p.j )
        }
    }
    ascend.tree(node.i)
    m
}


draw.delta.lines <- function(m, cols=NA, ...){
    x <- c(0, m[,'pg'])
    y <- c(0, m[,'int.d'])
    tel.b <- cumsum( m[,'node'] == 266 )
    mammal.b <- cumsum( m[,'node'] == 215 ) ## does not include marsupials
    aves.b <- cumsum( m[,'node'] == 268 )
    if(is.na(cols))
        cols <- rgb(tel.b, aves.b, mammal.b)
    l <- length(x)
    segments( x[1:(l-1)], y[1:(l-1)], x[2:l], y[2:l], col=cols, ... )
    nr <- nrow(x)
    c('x'=x[l], 'y'=y[l])
}

## plot these values
dev.set(3)

xlim <- range( unlist(lapply( leaf.int.d, function(x){ x[,'pg'] })) )
xlim[1] <- 0
ylim <- range( unlist(lapply( leaf.int.d, function(x){ x[,'int.d'] })) )

plot.new()
plot.window(xlim=xlim, ylim=ylim)
abline(h=0, col='grey')
ends <- sapply(leaf.int.d, draw.delta.lines)
axis(2)

identify(ends[1,], ends[2,], labels=colnames( int.s.l) )

draw.delta.lines( leaf.int.d[['hippocampus comes']], cols=rgb(0, 0, 1, 0.7) )
draw.delta.lines( leaf.int.d[['denticeps clupeoides']], cols=rgb(0, 1, 0, 0.7) )
draw.delta.lines( leaf.int.d[['hucho hucho']], cols=rgb(0, 1, 0, 0.7), lwd=2 )
draw.delta.lines( leaf.int.d[['callorhinchus milii']], cols=rgb(1, 0.5, 0, 1), lwd=2 )
draw.delta.lines( leaf.int.d[['erpetoichthys calabaricus']], cols=rgb(1, 0.5, 0, 1), lwd=2 )
draw.delta.lines( leaf.int.d[['latimeria chalumnae']], cols=rgb(1, 0, 0.5, 1), lwd=2 )


## and that is really very nice. We see teleosts specifically going
## down, and down, and we also see the aves branch specifically
## decreasing after divergence from non-aves lineages.

## here we see a good correlation between genome size and y-position (which
## is the same thing as seeing a correlation between genome size and mean
## intron size). In particular we see that Betta splendens is close to takifugu
## although it is not _that_ closely related by phylogeny. Well to be fair,
## it is actually pretty close as it comes off the branch point above the
## takifugu / tetraodon.
## The pattern observed, suggests, a random reduction in teleost intron size
## if that is the case we should see a reasonable proportion of short introns
## that are not shared betweed betta and takifugu, or betta and gasterosteus.

plot(int.s.l[,'tetraodon nigroviridis'], int.s.l[,'betta splendens'], cex=0.5, col=rgb(0,0,0,0.3) )
plot(int.s.l[,'tetraodon nigroviridis'], int.s.l[,'takifugu rubripes'], cex=0.5, col=rgb(0,0,0,0.3) )
plot(int.s.l[,'tetraodon nigroviridis'], int.s.l[,'anabas testudineus'], cex=0.5, col=rgb(0,0,0,0.3) )
plot(int.s.l[,'tetraodon nigroviridis'], int.s.l[,'maylandia zebra'], cex=0.5, col=rgb(0,0,0,0.3) )
plot(int.s.l[,'tetraodon nigroviridis'], int.s.l[,'hippocampus comes'], cex=1, col=rgb(0,0,0,0.3) )
plot(int.s.l[,'tetraodon nigroviridis'], int.s.l[,'denticeps clupeoides'], cex=1, col=rgb(0,0,0,0.3) )
plot(int.s.l[,'tetraodon nigroviridis'], int.s.l[,'danio rerio'], cex=1, col=rgb(0,0,0,0.3) )

## We think that there are a set of long introns that cannot reduce to a certain size,
## and others. In tetraodon and takifugu, most introns which can become smalll have been.
## however, not all, as can be seen from the plot between denticeps clupeoides and tetraodon
## here there are a a good number of introns which are long in tetraodon, but short in
## denticeps. It would seem reaonable for this set of introns to include the full set
## of introns. I.e. that these would default to some sort of reasonable minimum size.
##

## We may also wish to ask, what is the distribution of minimum sizes?
## Consider only lenths above 64 bases as smaller than that looks to be artificial.

int.s.min <- apply( int.s.l, 1, function(x){ min( x[x >= 6], na.rm=TRUE ) })
## the horrid addition of sapply to get 0s instead of length 0 
int.s.min.i <- sapply( apply( int.s.l, 1, function(x){ which.min( x[x >= 6] ) }),
                      function(x){ ifelse(length(x), x, 0) })

hist(int.s.min)
sum( int.s.min < 8 ) / length( int.s.min )  ## 94% !
sum( int.s.min < 7 ) / length( int.s.min )  ## 88% !

par(mar=c(14.1, 4.1, 4.1, 2.1))
barplot( table( c(1:172, int.s.min.i)) - 1, las=2, names.arg=c('NA', colnames(int.s.l)) )

int.s.max <- apply( int.s.l, 1, function(x){ max( x[x >= 6], na.rm=TRUE ) })
## the horrid addition of sapply to get 0s instead of length 0 
int.s.max.i <- sapply( apply( int.s.l, 1, function(x){ which.max( x[x >= 6] ) }),
                      function(x){ ifelse(length(x), x, 0) })

hist(int.s.max)

par(mar=c(14.1, 4.1, 4.1, 2.1))
barplot( table( c(1:172, int.s.max.i)) - 1, las=2, names.arg=c('NA', colnames(int.s.l)) )

## The great variance in minimal intron size is somewhat surprising; but it argues for
## a continuous process of intron size minimisation in teleosts... or at least continuous
## change.

### We can also ask how do the states change through an inferred path of evolution.
### Consider

##tmp <- collect.ancestor.states( 'danio rerio' )
##tmp <- collect.ancestor.states( 'takifugu rubripes' )
## tmp <- collect.ancestor.states( 'betta splendens' )
tmp <- collect.ancestor.states( 'poecilia formosa' )
## this gives us 6 columns rather than the 7 that we could get if we consider
## the root.
all.h <- hist(tmp)
ind.h <- apply(tmp, 2, hist, breaks=all.h$breaks, plot=FALSE)

par(mfrow=c(1,1))
plot( all.h$mids, all.h$counts, type='n', ylim=range(sapply(ind.h, function(x){ x$counts })) )
cols <- hsvScale(1:ncol(tmp))
for(i in 1:ncol(tmp))
    lines(all.h$mids, ind.h[[i]]$counts, col=cols[i], lwd=2)

## can we do someting silly like do a pca on the full set here and see what happens..
## maybe exclude the acutal leaf as that has a lot of 0s in it

tmp.pca <- prcomp(tmp, scale=TRUE)

tmp.pca <- prcomp(t(tmp[,-ncol(tmp)]))
plot(tmp.pca)
plot(tmp.pca$x[,1], tmp.pca$x[,2], cex=2.5)
text(tmp.pca$x[,1], tmp.pca$x[,2])
o <- order(tmp.pca$rot[,1])
image(tmp[o,])
image(tmp[o[2000:3000], ])
image(tmp[rev(o)[10000:11000], ])

o <- order(tmp.pca$rot[,2])

### Is there a tendency for small introns to stay small
### or for large introns not to grow. We can collect
### state transitions from a given node by recursing 
### down the tree and collecting transitions
### We are using the ex.align.2.k2.nj tree in this function
### The edge component has parents in column 1 and child nodes
### in column 2

collect.transitions <- function( root, edges=ex.align.2.k2.nj$edge, inf=int.s.inf,
                                dummy.r=18){
    if(!is.numeric(root))
        stop("Unable to get a numeric node identifier from: ", root)
    tr.tables <- vector(mode='list', length=0)
    dummy.d <- -dummy.r:dummy.r
    dummy.r <- abs( dummy.d )
    dummy.t <- table( dummy.r, dummy.d )
    col.i <- as.numeric(colnames(dummy.t))
    row.i <- as.numeric(rownames(dummy.t))
    loss.b <- t( sapply( row.i, function(r){
        r + col.i <= 0 }))
    descend.tree <- function(i){
        c.i <- which( edges[,1] == i )
        if(!length(c.i))
            return
        for(j in edges[c.i,2]){
            d <- inf[[j]]$state[,1] - inf[[i]]$state[,1]
            tbl <- table( c(dummy.r, inf[[i]]$state[,1]),
                          c(dummy.d, d) )
            tbl <- tbl - dummy.t
            tbl[ loss.b ] <- 0
            tr.tables <<- c(tr.tables, list(tbl))
            descend.tree(j)
        }
    }
    descend.tree(root)
    tr.tables
}

## compare the root only to leaves
collect.transitions.2 <- function( root, edges=ex.align.2.k2.nj$edge, inf=int.s.inf,
                                dummy.r=18){
    if(!is.numeric(root))
        stop("Unable to get a numeric node identifier from: ", root)
    root.state <- inf[[root]]$state[,1]
    tr.tables <- vector(mode='list', length=0)
    tbl.names <- vector(mode='character', length=0)
    dummy.d <- -dummy.r:dummy.r
    dummy.r <- abs( dummy.d )
    dummy.t <- table( dummy.r, dummy.d )
    descend.tree <- function(i){
        c.i <- which( edges[,1] == i )
        if(!length(c.i)){
            b <- inf[[i]]$state[,1] > 0
            d <- inf[[ i ]]$state[b,1] - root.state[b]
            tbl <- table( c(dummy.r, root.state[b]),
                          c(dummy.d, d) )
            tbl <- tbl - dummy.t
            tr.tables <<- c(tr.tables, list(tbl))
            tbl.names <<- c(tbl.names, colnames(int.s.l)[i])
            return
        }
        for(j in edges[c.i,2]){
            descend.tree(j)
        }
    }
    descend.tree(root)
    names(tr.tables) <- tbl.names
    tr.tables
}


## note that the root nodes for:
## teleost: 266
## mammal: 215 (without marsupials)
## aves  : 268
## I don't quite remember how I found these, but I believe it is by plotting
## the whole tree. I need to do that as well

teleost.transitions <- collect.transitions( 266 )
mammal.transitions <- collect.transitions( 215, dummy.r=19 )

teleost.transitions.2 <- collect.transitions.2( 266 )
mammal.transitions.2 <- collect.transitions.2( 215, dummy.r=19 )

image.transitions <- function(trans, ...){
    image(x=as.numeric(colnames(trans)), y=as.numeric(rownames(trans)),  z=t(log(1 + trans)),
          ...)
}

plot.transitions <- function(transitions){
    tr.nm <- names(transitions)
    for(i in 1:length(transitions)){
        image.transitions( transitions[[i]], main=tr.nm[i] )
        d <- as.numeric(colnames(transitions[[i]]))
        s <- as.numeric(rownames(transitions[[i]]))
        counts <- sapply(1:nrow(transitions[[i]]), function(j){
            b <- d + s[j] > 0
            c('neg'= sum( transitions[[i]][j, d < 0 & b ]),
              'pos'= sum( transitions[[i]][j, d > 0 & b ]),
              'all'= sum( transitions[[i]][j, ] ))
            })
        colnames(counts) <- rownames(transitions[[i]])
        print(counts)
        print( cbind(sprintf("%.2f", counts['neg',] / counts['all',]),
                     sprintf("%.2f", counts['pos',] / counts['all',])) )
        inpt <- readline("next: ")
    }
}

plot.transitions(teleost.transitions)
plot.transitions(mammal.transitions)

plot.transitions(teleost.transitions.2)
plot.transitions(mammal.transitions.2)

## these are quite variable and not so easy to conclude things from;
## a simpler quesion might be to compare the extant nodes with a single
## suitable ancestor node.

### some formal testing of gene ontology enrichment.
### The problem we have is that genes with many introns are more likely
### to be picked. ?? There are packages that take this into control, but
### the ones I have found seem to be rather too specific to gene expression
### analysis or similar.

## Start off with simply using GOstats
## and the org.hs.db



if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("GOstats")
BiocManager::install("org.Hs.eg.db")

require("GOstats")
require("org.Hs.eg.db")


ens2eg <- as.list(org.Hs.egENSEMBL2EG)
table( sapply( ens2eg, length ))

##     1     2     3     4     5     7 
## 28813   364    22     6     5     1 

eg2ens <- as.list(org.Hs.egENSEMBL)


## I have human ids in
## intron.orth.t[[1]]$id[,'homo sapiens']
## but note that not all of the rows are present

hum.g.na <- is.na(intron.orth.t[[1]]$id[,'homo sapiens'])

sum(hum.g.na)
## [1] 303
sum(hum.g.na) / length(hum.g.na)
## [1] 0.004804338
sum(!hum.g.na)
## 62765

## this is counting many times
length(unique( intron.orth.t[[1]]$id[!hum.g.na,'homo sapiens']))
## 5708

hum.g <- unique( intron.orth.t[[1]]$id[!hum.g.na,'homo sapiens'])
sum( hum.g %in% names( ens2eg ) )
## 5695  ## this will be our universe..

hum.universe <- unlist(ens2eg[ hum.g ])

## we then have some variance estimates in
## teleost.var, non.teleost.var, all.var

do.hp.test <- function(universe, sel.ens, ont='BP'){
    sel.eg <- unlist(ens2eg[ sel.ens ])
    params <- new("GOHyperGParams",
                  geneIds=sel.eg,
                  universeGeneIds=universe,
                  annotation="org.Hs.eg.db",
                  ontology=ont,
                  conditional=FALSE,
                  testDirection="over")
    hyperGTest(params)
}

sel.ens <- unique( intron.orth.t[[1]]$id[ all.var[,'mean'] > 10 &
                                          all.var[,'sd'] < 1, 'homo sapiens' ])
sel.ens <- sel.ens[ !is.na(sel.ens) ]

tmp <- do.hp.test( hum.universe, sel.ens, ont='BP')
## and this gives me e-19, but with a reasonably small odds ratio (1.9)

mean.l.bp <- lapply( 6:15, function(x){
    sel.ens <- unique( intron.orth.t[[1]]$id[ all.var[,'mean'] > x &
                                              all.var[,'sd'] < 1, 'homo sapiens' ])
    sel.ens <- sel.ens[ !is.na(sel.ens) ]
    list('sel'=sel.ens, 'hg'=do.hp.test( hum.universe, sel.ens, ont='BP'))
})

min.l.bp <- lapply( 5:13, function(x){
    sel.ens <- unique( intron.orth.t[[1]]$id[ all.var[,'min'] > x &
                                              all.var[,'sd'] < 1, 'homo sapiens' ])
    sel.ens <- sel.ens[ !is.na(sel.ens) ]
    list('sel'=sel.ens, 'hg'=do.hp.test( hum.universe, sel.ens, ont='BP'))
})


min.l.bp.2 <- lapply( 5:13, function(x){
    sel.ens <- unique( intron.orth.t[[1]]$id[ all.var[,'min'] > x
                                              , 'homo sapiens' ])
    sel.ens <- sel.ens[ !is.na(sel.ens) ]
    list('sel'=sel.ens, 'hg'=do.hp.test( hum.universe, sel.ens, ont='BP'))
})


tel.min.l.bp <- lapply( 5:13, function(x){
    sel.ens <- unique( intron.orth.t[[1]]$id[ teleost.var[,'min'] > x
                                              , 'homo sapiens' ])
    sel.ens <- sel.ens[ !is.na(sel.ens) ]
    list('sel'=sel.ens, 'hg'=do.hp.test( hum.universe, sel.ens, ont='BP'))
})


mammal.min.l.bp <- lapply( 5:15, function(x){
    sel.ens <- unique( intron.orth.t[[1]]$id[ mammal.var[,'min'] > x
                                              , 'homo sapiens' ])
    sel.ens <- sel.ens[ !is.na(sel.ens) ]
    list('sel'=sel.ens, 'hg'=do.hp.test( hum.universe, sel.ens, ont='BP'))
})


mammal.max.l.bp <- lapply( 5:15, function(x){
    sel.ens <- unique( intron.orth.t[[1]]$id[ mammal.var[,'max'] < x
                                              , 'homo sapiens' ])
    sel.ens <- sel.ens[ !is.na(sel.ens) ]
    list('sel'=sel.ens, 'hg'=do.hp.test( hum.universe, sel.ens, ont='BP'))
})


## to get all genes annotated with a specific group:
go2eg <- as.list(org.Hs.egGO2ALLEGS)

## from :
lapply( tel.min.l.bp, function(x){ head(summary(x$hg)) })

## we have the strongest enrichment from entry #4 ( x > 8 )
head(summary(tel.min.l.bp[[4]]$hg), n=10)

##        GOBPID       Pvalue OddsRatio  ExpCount Count Size
## 1  GO:0023052 1.349236e-24  1.835704  850.7532  1023 1689
## 2  GO:0007154 1.888800e-24  1.828819  859.8198  1032 1707
## 3  GO:0007165 1.217057e-22  1.807112  795.3459   957 1579
## 4  GO:0048731 5.581517e-18  1.704816  733.8942   873 1457
## 5  GO:0050794 7.821533e-18  1.616900 1559.4623  1712 3096
## 6  GO:0048856 8.997116e-18  1.650375  890.0420  1036 1767
## 7  GO:0007275 1.128844e-17  1.666048  821.5384   964 1631
## 8  GO:0032501 5.955897e-17  1.600669 1049.7156  1197 2084
## 9  GO:0007155 8.630577e-17  2.701575  166.2218   238  330
## 10 GO:0010646 1.183684e-16  1.779536  525.3615   644 1043
##                                  Term
## 1                           signaling
## 2                  cell communication
## 3                 signal transduction
## 4                  system development
## 5      regulation of cellular process
## 6    anatomical structure development
## 7  multicellular organism development
## 8    multicellular organismal process
## 9                       cell adhesion
## 10   regulation of cell communication

## we could also try:
plot(-log10(summary(tel.min.l.bp[[4]]$hg)$Pvalue))

b <- summary(tel.min.l.bp[[4]]$hg)$Pvalue <= 1e-10
sum(b)
## [1] 35
## which is quite a lot; and maybe too many.. 

## anyway, simply do
tel.min.l.bp.sum <- summary( tel.min.l.bp[[4]]$hg )

## and for looking at the other way around, lets do
mammal.max.l.bp.sum <- summary( mammal.max.l.bp[[7]]$hg )

## we can then get all of the entrez gene identifiers associated with
## the set of these terms:

tel.min.l.bp.go.eg <- lapply( tel.min.l.bp.sum$GOBPID, function(x){ go2eg[ x ] })
sapply( tel.min.l.bp.go.eg, function(x){ length(x[[1]] ) } )

 ##  [1] 12425 12494 11421  7782 23523  9773  8693 13254  2061  5503  2075  5574
 ## [13] 28650 26560 14754 10666  5050  4889  5047 19322  3964  4606  3780  7047
 ## [25]  2471  6493  4874  1318  6167  3412  2948  6276  5046  1458  2810  7886
 ## [37]  7462  2838   880  2815  9812  2334  1465  9444  1907 11321  2451   917
 ## [49]  2164  1912  1787   919  4087  1118  4011  1443  2052  2074  1026   940
 ## [61]   971   832  1131   853   849  3589  1623  1274  3158  2241   558  1845
 ## ...

## these are generally very large groups that we see
## but they are likely to overlap very much with each other.
## Given that genes have many introns it would make sense to take the maximal
## size per gene and then visualise this..

mammal.max.l.bp.go.eg <- lapply( mammal.max.l.bp.sum$GOBPID, function(x){ go2eg[ x ] })
sapply( mammal.max.l.bp.go.eg, function(x){ length(x[[1]] ) } )
## these are generally smaller groups.. 


## give some statistics about enrichment or depletion
gene.set.association <- function(eg, intron.stat=teleost.var[,'min'],
                                 ids=intron.orth.t[[1]]$id,
                                 min.l=5, int.criteria=which.max,
                                 perm.n=10, decreasing=FALSE,
                                 entrez=TRUE){
    if(entrez){
        ens.id <- unlist( eg2ens[eg] )
        ens.id <- ens.id[!is.na(ens.id)]
    }else{
        ens.id <- eg
    }
    ##
    ## obtain the indices of the longest intron for each gene (ens.i)
    ens.i <- tapply( 1:length(intron.stat), ids[,'homo sapiens'], function(i){
        i[ int.criteria( intron.stat[i] )] })
    ## then get the length of the longest intron
    int.s <- intron.stat[ens.i]
    ## get the human ensembl ids (int stands for intron here)
    int.id <- ids[ens.i, 'homo sapiens']
    ## very short introns are probably unreasonable; ignore these
    ## removing both the ids and the intron sizes
    int.id <- int.id[int.s >= min.l]
    int.s <- int.s[int.s >= min.l] 
    ## reorder by intron size from small to big
    o <- order(int.s, decreasing=decreasing)
    int.id <- int.id[o]
    int.s <- int.s[o]
    
    ## which of the introns are in eg?
    o.b <- int.id %in% ens.id
    ## the following will give the number of genes in the set (eg)
    ## with an increase in size
    get.stats <- function(b){
        q <- cumsum( b )
        p <- phyper(q, sum(b), sum(!b), 1:length(b))
        p2 <- phyper(q-1, sum(b), sum(!b), 1:length(b), lower.tail=FALSE)
        exp <- 1:length(b) * sum(b) / length(b)
        data.frame('id'=int.id, 'l'=int.s, 'b'=b, 'q'=q, 'p'=p, 'p2'=p2, 'exp'=exp )
    }
    stats <- get.stats(o.b)
    stats.permed <- lapply( 1:perm.n, function(i){
        get.stats( sample(o.b) )
    })
    list('o'=stats, 'o.perm'=stats.permed)
}

## takes an object as returned by gene.set.association
plot.association.stats <- function(stats.l, min.sample=50,
                                   cols=c(1, 2, 4, 6), p2=FALSE,
                                   draw.permuted=FALSE,
                                   draw.min.p=FALSE, ...){
    ## only plot points where a gene is a member;
    ## makes the plot a little bit cleaner
    stats <- stats.l$o
    stats.p <- stats.l$o.perm
    b <- stats$b
    b[1:min.sample] <- FALSE
    ## plot observed vs expected
    par(mar=c(5.1, 4.1, 4.1, 8.1))
    plot( stats$l[b], stats$q[b] / stats$exp[b],
         xlab='log2 intron size', ylab='', ##  'observed / expected',
         col=cols[1], type='l', axes=FALSE, ...)
    axis(1)
    axis(2, col=cols[1], col.axis=cols[1])
    ## lets do the permuted values
    if(draw.permuted){
        invisible( lapply(stats.p, function(x){
            lines(x$l[ x$b ], x$q[ x$b ] / x$exp[ x$b ], col=rgb(0.7, 0.7, 0.7))
        }))
    }
    ## then plot p-values
    if(!p2)
        pp <- -log10( stats$p[b] )
    else
        pp <- -log10( stats$p2[b] )
##            
    plot.window( xlim=range( stats$l[b] ), ylim=range(pp) )
    lines( stats$l[b], pp, col=cols[2] )
    axis(4, line=0, col=cols[2], col.axis=cols[2], )
    if(draw.min.p){
        min.p.i <- which.max( pp )
        abline(v=(stats$l[b])[min.p.i], col=cols[2], lty=2)
    }
    ## then plot n (or rather q)
    plot.window( xlim=range( stats$l[b] ), ylim=range(stats$q[b]) )
    lines( stats$l[b], stats$q[b], type='l', col=cols[3] )
    axis(4, line=2.5, col=cols[3], col.axis=cols[3])
    plot.window( xlim=range( stats$l[b] ), ylim=c(0, nrow(stats) ))
    lines( stats$l, 1:nrow(stats), col=cols[4], lty=2, lwd=2 )
    axis(4, line=5, col=cols[4], col.axis=cols[4])
    plot.window(xlim=c(0,1), ylim=c(0,1))
    legend('bottomright', legend=c('obs/exp', '-log10 p', 'q', 'n'),
           text.col=cols, box.lty=0, lty=c(1,1,1,2), col=cols, lwd=3, inset=c(0.0, 0.1))
}

### These plots should be good after a bit of cosmetic cleanup...
### we have peaks of depletion betwen 8 and 9 for each of these.
### I should perhaps include expected / observed ratios for each?

intron.depletion <- lapply( tel.min.l.bp.go.eg, function(x){
    gene.set.association( x[[1]] ) })

## the following reversed the sort and so calculates statistics
## for the right-hand subsets; where we see enrichment of the
## same set of genes. In this case it might make more sense to
## 
intron.depletion.2 <- lapply( tel.min.l.bp.go.eg, function(x){
    gene.set.association( x[[1]], decreasing=TRUE ) })

intron.depletion.3 <- lapply( tel.min.l.bp.go.eg, function(x){
    gene.set.association( x[[1]], decreasing=TRUE, int.criteria=which.min ) })


for(i in 1:length(intron.depletion)){
    plot.association.stats( intron.depletion[[i]],
                          main=paste(tel.min.l.bp.sum[i, c('GOBPID', 'Term')],
                                   collapse=" "), min.sample=100,
                          draw.min.p=TRUE)
    inpt <- readline('next: ')
}

for(i in 1:length(intron.depletion.2)){
    plot.association.stats( intron.depletion.2[[i]],
                          main=tel.min.l.bp.sum$Term[i], min.sample=100, p2=TRUE,
                          draw.min.p=TRUE)
    inpt <- readline('next: ')
}

for(i in 1:length(intron.depletion.3)){
    plot.association.stats( intron.depletion.3[[i]],
                          main=tel.min.l.bp.sum$Term[i], min.sample=100, p2=TRUE,
                          draw.min.p=TRUE)
    inpt <- readline('next: ')
}


## this seems to do the wrong thing. I need to go through
## the function and make it work for enrichment as well as depletion
mammal.intron.stats <- lapply( mammal.max.l.bp.go.eg, function(x){
    gene.set.association( x[[1]], intron.stat=mammal.var[,'max'], int.criteria=which.min) })

## ask the same question, but use the biggest intron per gene
mammal.intron.stats.2 <- lapply( mammal.max.l.bp.go.eg, function(x){
    gene.set.association( x[[1]], intron.stat=mammal.var[,'max'], int.criteria=which.max) })
## using the largest intron loses us the enrichment. Lets try with the mean size
## as well. Should really have the median as well


for(i in 1:length(mammal.intron.stats)){
    plot.association.stats( mammal.intron.stats[[i]],
                          main=mammal.max.l.bp.sum$Term[i], min.sample=100, p2=TRUE )
    inpt <- readline('next: ')
}
## The above actually has some potentially interesting plots, with an enrichment for
## short introns for some of these

par(mfrow=c(1,2))
for(i in 1:length(mammal.intron.stats.2)){
    plot.association.stats( mammal.intron.stats[[i]],
                          main=mammal.max.l.bp.sum$Term[i], min.sample=100, p2=TRUE )
    plot.association.stats( mammal.intron.stats.2[[i]],
                          main=mammal.max.l.bp.sum$Term[i], min.sample=100, p2=TRUE )
    inpt <- readline('next: ')
}



## Ensembl has entries of regulatory regions of different types; We can thus test
## the hypothesis that introns that are long across teleosts are more likely to
## contain enhancer elements in humans. That can be controlled for length as
## the majority of introns are long in humans...


### I have now incorporated the ensembl transcript identifiers into the
### intron.orth.t tables; I can thus get the intron sequences from the different
### species and ask whether we can identify conserved regions;
##
## Use:
## teleost.var            --> identify rows of interest
## intron.orth.t[[1]]$tr  --> get transcript identities
## intron.orth.fam[[1]]   --> get the family identity
##                        construct the filename containing
## intron.f               --> read.exons( intron.f )
## For each sequence align to the danio rerio sequence using
## local.aligns() in /home/lmj/R/exon_aligneR/functions.r
## 
## We ned a substitution matrix
## I have unfortunately used the same name more than once

sub.matrix <- make.sub.matrix()

## lets use PTHR24083_SF46
## as an example; that has between 2^10 and 2^12.7
## in teleosts. Actually a bit smaller in mammals
## between 8.7 -> 13. So up to 16,384 bases.

cols <- c(rgb(0.9, 0.9, 0.9), rgb(1, 0, 1),
                    rgb(c(1,0,0,0.5), c(0,1,0,0.5),
                        c(0,0,1,0)) )
names(cols) <- c('-', 'I', 'A', 'C', 'G', 'T')

## max.l is the maximum length of the score tables
align.introns <- function(i, seed.sp='danio rerio', al.sp=seed.sp,
                          gap=c(-10,-1), min.width=15, min.score=40,
                          sm=sub.matrix, max.l=1e7){
    fam.id <- intron.orth.fam[[seed.sp]][i]
    tr <- intron.orth.t[[seed.sp]]$tr[i,]
    int.i <- intron.orth.t[[seed.sp]]$i[i,]
    
    int.seq <- read.exons( intron.f[ fam.id ] )
    int.meta <- sapply( int.seq, function(x){ c('id'=x$id, 'tr'=x$tr, 'sp'=x$sp, 'n'=length(x$l) )})
    ## select the sequences that we will use
    b <- int.meta['tr',] %in% tr
    int.seq <- int.seq[b]
    int.meta <- int.meta[,b]
    ## so we can order by species
    names(int.seq) <- int.meta['sp',]
    colnames(int.meta) <- int.meta['sp',]
    
    int.i <- int.i[ names(int.i) %in% names(int.seq) ]
    al.i <- which( names(int.i) == al.sp )
    a.seq <- int.seq[[al.sp]]$e[ int.i[al.sp] ]
    b.seq <- sapply( names(int.i)[-al.i], function(sp){
        int.seq[[sp]]$e[ int.i[sp] ] })
    b.seq <- b.seq[ !is.na(b.seq) ]
    ## we could define a reasonable order here using ex.align.2.k2
    sp.o <- names( sort( ex.align.2.k2[ sub(" ", "_", al.sp), ]))
    sp.o <- sub("_", " ", sp.o)
    sp.o <- sp.o[ sp.o %in% names(b.seq) ]
    aligns <- vector(mode='list', length=length(sp.o))
    for(j in 1:length(sp.o)){
        sp <- sp.o[j]
        if(as.double(nchar(a.seq)) * as.double(nchar(b.seq[sp])) < max.l )
            aligns[[j]] <- local.aligns( a.seq, b.seq[sp], sm$offset,
                                        sm$size, sm$sm,
                                        gap=gap, min.width=min.width, min.score=min.score )
    }
    list(sp1=al.sp, sp=sp.o, al=aligns, l1=nchar(a.seq), l2=unname(sapply(b.seq, nchar)) )
}

## we can test this with o[1088]
## --> 55047 (COUPTF)

al.test <- align.introns( 55047 )

## select a set of interesting introns;
b <- teleost.var[,'min'] > 10 & !is.na(teleost.var[,'sd']) & teleost.var[,'n'] > 10  ## 1594 introns;
## a bit too many, but, let us take the top 100 with the least
## amount of variance
b.i <- which(b)
o <- order( teleost.var[b,'sd'] )

## lets do 50 alignments:
## we get a seg fault at the 13th of these; lets try to work out why..
## PTHR24373_SF15

### This will be too slow to do here and we are starting to use stupid amounts of memory;
### Instead export the relevant information; read these into a perl script and
### call Rscript on a suitable R script that runs the alignments
### We can then use fork() to parallelize.

## or we can try using the parallel library and see how that works out

write.table(intron.orth.fam[[1]], "dr_intron_fam.txt",
            sep="\t", quote=FALSE)
write.table(intron.orth.t[[1]]$id, "dr_intron_orthology_id.txt",
            row.names=FALSE, sep="\t", quote=FALSE)
write.table(intron.orth.t[[1]]$i, "dr_intron_orthology_i.txt",
            row.names=FALSE, sep="\t", quote=FALSE)
write.table(intron.orth.t[[1]]$tr, "dr_intron_orthology_tr.txt",
            row.names=FALSE, sep="\t", quote=FALSE)

write.table(intron.orth.t[[1]]$l, "dr_intron_orthology_l.txt",
            row.names=FALSE, sep="\t", quote=FALSE)

write.table(teleost.var, "dr_teleost_var.txt",
            row.names=FALSE, sep="\t", quote=FALSE)

write.table(mammal.var, "dr_mammal_var.txt",
            row.names=FALSE, sep="\t", quote=FALSE)

write.table(sauria.var, "dr_sauria_var.txt",
            row.names=FALSE, sep="\t", quote=FALSE)

write.table(intron.f, 'intron_files.txt',
            sep="\t", quote=FALSE)

write.table(teleost.b, "teleost_b.txt", sep="\t", row.names=TRUE, quote=FALSE)
write.table(mammal.b, "mammal_b.txt", sep="\t", row.names=TRUE, quote=FALSE)
write.table(sauria.b, "sauria_b.txt", sep="\t", row.names=TRUE, quote=FALSE)

## debug(align.introns)
## align.introns( 58887 )

## intron.al.t50 <- vector(mode='list', length=50)
## for(i in 1:50){
##     j <- b.i[o[i]]
##     print(paste("aligning introns for family: ", intron.orth.fam[[1]][j]))
##     intron.al.t50[[i]] <- list('i'=j, align.introns(j) )
## }
