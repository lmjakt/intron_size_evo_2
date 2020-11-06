## First ask the question:
##
## Do introns that are long across the teleosts contain fewer repeats
## than those that have variable length.
##
## This is to some extent not an interesting question because the answer
## should follow from the expansion of intron size in D.r. from repeat
## spreading. And that has been published. But we can then go on to use
## the masking to refine the intron alignments in a way that may be better
## than the blast based one I used previously.

mask <- read.table( 'danio_rerio_introns_mask.txt', header=FALSE, sep="\t", stringsAsFactors=FALSE,
                   colClasses=c('character', 'integer', 'character', 'character', 'integer', 'integer', 'integer', 'character',
                                'character', 'character', 'character', 'character'))
colnames(mask) <- c('chr', 'strand', 'gene', 'tr', 'rank', 'start', 'end', 'reps', 'scores', 'r.start', 'r.end', 'mask')

orth <- readRDS('../R_intron_alignments_summary_2/orth.rds')
tel.var <- readRDS('../R_intron_alignments_summary_2/tel_var.rds')
mam.var <- readRDS('../R_intron_alignments_summary_2/mam_var.rds')
sp.class <- readRDS('../R_trees_distances/sp_class_3.rds')
int.s.b <-  readRDS('../R_trees_distances/int_s_b.rds')
## colnames in orth are 'genus.species'
## species names in sp.class and int.s.b are genus_species

sum( orth$id$danio.rerio %in% mask$gene ) - nrow(orth$id)
## 0 ? entries without introns?
dim(mask)

## [1] 229159     12

## we wish to subset the two:

i <- match( paste( orth$id[,'danio.rerio'], orth$i[,'danio.rerio'], sep="_" ), paste( mask$gene, mask$rank, sep="_" ) )

mask.m <- mask[i[!is.na(i)], ]
## those now seem to match properly. I can now go ahead and ask simple questions.
## Well, maybe not so simple.

## lengths of masks and positions
mask.m.l <- t(sapply(mask.m$mask, function(x){
    m <- utf8ToInt( x ) - 48
    ## 0 -> not masked
    ## 1 -> masked
    c('l'=length(m), 'ml'=sum(m))
}))
rownames(mask.m.l) <- NULL

## I had some bugs in my script, but this is now true
## after ironing those out.
all( mask.m.l[,'l'] == orth$l[,'danio.rerio'] )
## [1] TRUE

plot(mask.m.l[,'l'], mask.m.l[,'ml'])
plot(log2(mask.m.l[,'l']), log2(mask.m.l[,'ml']), cex=0.4)
abline(v=8, col='red')

plot(log2(mask.m.l[,'l']), mask.m.l[,'ml'] / mask.m.l[,'l'], cex=0.4)

hist( mask.m.l[,'ml'] / mask.m.l[,'l'] )

plot(tel.var[,'min'], mask.m.l[,'ml'] / mask.m.l[,'l'], cex=0.4)
plot(tel.var[,'10%'], mask.m.l[,'ml'] / mask.m.l[,'l'], cex=0.4)

plot(tel.var[,'10%'], log2(mask.m.l[,'l'] -  mask.m.l[,'ml'] ), cex=0.4)
abline(0, 1, col='red')

par(mfrow=c(1,2))
plot(log2(orth$l[,'danio.rerio']), log2(orth$l[,'gasterosteus.aculeatus']), cex=0.4)
abline(0,1, col='red')
plot(log2(mask.m.l[,'l'] -  mask.m.l[,'ml'] ), log2(orth$l[,'gasterosteus.aculeatus']), cex=0.4)
abline(0,1, col='red')

par(mfrow=c(1,2))
plot(log2(orth$l[,'danio.rerio']), log2(orth$l[,'astyanax.mexicanus']), cex=0.4)
abline(0,1, col='red')
plot(log2(mask.m.l[,'l'] -  mask.m.l[,'ml'] ), log2(orth$l[,'astyanax.mexicanus']), cex=0.4)
abline(0,1, col='red')

par(mfrow=c(1,2))
plot(tel.var[,'10%'], log2(orth$l[,'danio.rerio']), cex=0.4, xlim=c(6,18), ylim=c(6,18))
abline(0,1, col='red')
plot(tel.var[,'10%'], log2(mask.m.l[,'l'] -  mask.m.l[,'ml'] ), cex=0.4, xlim=c(6,18), ylim=c(6,18))
abline(0,1, col='red')
##
## This suggests that a reasonable proportion of the increas in intron length in danio rerio
## comes from the proliferation of repetitive DNA. That is not surprising, and what we would
## expect to see.

###### To make this interesting we want to have a look at the alignments we obtained and to see
###### what proportion overlap with repeat regions. That will be a bit slow.


