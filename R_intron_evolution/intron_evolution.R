## Do we really have common retention of long introns?
## My previous analysis was flawed in that the probability
## of an intron changing size from
##    a --> b 
## where b is less than a given threshold
## is clearly dependant on the size in a, and we can
## we can see this by simulation.

## Here I want to compare the results of simulation over specified windows
## of a (size in ancestor) against those observed in the real
## data.

## mainly for hsvScale
source("~/R/general_functions.R")
source("ie_functions.R")
require(parallel)

## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2

half.w <- a4.w * 85 / 210
full.w <- a4.w * 170 / 210

## lets get some genome sizes
tmp <- read.table("../R_172_genomes/genome_sizes.txt")
genome.size <- tmp[,1]
names(genome.size) <- rownames(tmp)
genome.size <- sort(genome.size)
rm(tmp)

## Some species classification:
sp.class <- readRDS( "../R_trees_distances/sp_class_3.rds")
class.col <-  readRDS( "../R_trees_distances/class_col_3.rds")
sp.col <- readRDS("../R_trees_distances/sp_col_3.rds")

## we can also get species with too many too short introns
int.s.b <- readRDS("../R_trees_distances/int_s_b.rds")

## for plotting we can use a taxonomic order obtained from the
## phylogenetic tree:
sp.y <- readRDS("../R_trees_distances/sp_y.rds")

sp.names <- sub("_", " ", names(sp.y))
substring(sp.names, 1, 1) <- toupper(substring(sp.names, 1, 1))
names( sp.names ) <- names(sp.y)

sp.labels <- sapply( strsplit(sp.names, " ", fixed=TRUE), function(x){
    paste(substring(x[1], 1, 1), " ", x[2], sep="")})


ex.align.2.k2.nj <- readRDS( "../R_172_genomes/ex_align_2_k2_nj.rds" )
ex.align.2.k2 <- readRDS("../R_172_genomes/ex_align_2_k2.rds")
## we may not need to change labels; check after for consistency
## ex.align.2.k2.nj$tip.label <- sub("_", ".", ex.align.2.k2.nj$tip.label)

int.s.inf.2 <- readRDS("../R_trees_distances_2/int_s_inf_2.rds")
tree.lines <- readRDS("../R_trees_distances/tree_lines.rds")
clade.nodes <- readRDS("../R_trees_distances/clade_nodes.rds")

## this depends on internal ordering of nodes, but seems to work ok.
nodes.col <- rep('black', length(clade.nodes$mammalia))

for(clade in names(clade.nodes))
    nodes.col[ clade.nodes[[clade]] ] <- class.col[clade]

### I probably won't need the orthology, but it may be useful
### at some point. If I do need it, uncomment the following and
### run it:
orth <- readRDS( "../R_intron_alignments_summary_2/orth.rds")
## unfortunately this uses genus.species whereas everything else here
## uses genus_species
for(nm in c('id', 'tr', 'i', 'l'))
    colnames(orth[[nm]]) <- sub(".", "_", colnames(orth[[nm]]), fixed=TRUE)

## we should test that the simulation gives us similar end distributions.
## test for 266 -> 78 (D. rerio), 3 (B. splendens), 50 (Maylandia zebra)

sim.266.78 <- simulate.is.evol(266, 78, gen.n=1, inf=int.s.inf.2, use.blurred=FALSE )
with( int.s.inf.2[[78]], hist( state[ state[,1] > 10*log2(76), 1], breaks=100 ))
hist( sim.266.78[[2]]$size, breaks=100 )
hist( sim.266.78[[3]]$size, breaks=100 )

sim.266.3 <- simulate.is.evol(266, 3, gen.n=1, inf=int.s.inf.2, use.blurred=FALSE )
with( int.s.inf.2[[3]], hist( state[ state[,1] > 10*log2(76), 1], breaks=100 ))
hist( sim.266.3[[2]]$size, breaks=100 )
hist( sim.266.3[[3]]$size, breaks=100 )

sim.266.50 <- simulate.is.evol(266, 50, gen.n=1, inf=int.s.inf.2, use.blurred=FALSE )
with( int.s.inf.2[[50]], hist( state[ state[,1] > 10*log2(76), 1], breaks=100 ))
hist( sim.266.50[[2]]$size, breaks=100 )
hist( sim.266.50[[3]]$size, breaks=100 )

### and that does seem to do exactly what it should do.

## To test common retention of long between descendants of 267 (D. rerio and others)
## against descendants of 276 (all other teleosts )
## reduce the time do not test everything against everything, but have a selected
## set of leaves from 267

## this is for: E. electricus, D. clupeoides, S. formosus, I. punctatus
## and D. rerio
leaves.267.i <- c(9, 10, 28, 30, 78)
leaves.276 <- collect.leaves(276, tree=ex.align.2.k2.nj)

## do this using a loop rather than nested lapply
## as with a loop we will at least know where we are when we hit
## an error.

ci.267 <- vector(mode='list', length=length(leaves.267.i))
for(i in 1:length(ci.267)){
    ci.267[[i]] <- vector(mode='list', length=length(leaves.276$i))
    cat(i, "\n")
    ci.267[[i]] <- mclapply( leaves.276$i, function(j){
        test.common.retention(266, leaves.267.i[i], j, inf=int.s.inf.2, sim.n=5, thr=80, breaks=seq(80, 140, 10),
                              min.l=10*log2(76))},
        mc.cores=24)
}

cairo_pdf("common_length_retained_267_276_p.pdf", width=0.9 * a4.w * pdf.m, height=0.8 * a4.h * pdf.m, onefile=TRUE)
par(mfrow=c(5,4))
for(i in 1:length(ci.267)){
    sp1 <- sp.names[ ex.align.2.k2.nj$tip.label[ leaves.267.i[i] ] ]
    for(j in 1:length(ci.267[[i]])){
        sp2 <- sp.names[ ex.align.2.k2.nj$tip.label[ leaves.276$i[j] ] ]
        plot.tested.stat( ci.267[[i]][[j]], stat.name='p', sp1, sp2, cex.main=0.7 )
    }
}
dev.off()

## We can combine independent p-values using Fisher's method.
##
## -2 * sum( log(p[i]) )
## from which we can get a p-value using:
##
## pchisq( -2 * sum( log(p) ), df=2*length(p)
##
## And this seems to do something reasonable with our data.

ci.267.2 <- vector(mode='list', length=length(leaves.267.i))
for(i in 1:length(ci.267.2)){
    print(i)
    ci.267.2[[i]] <- mclapply( leaves.276$i, function(j){
        test.common.retention(266, leaves.267.i[i], j, inf=int.s.inf.2, sim.n=5, thr=80, breaks=seq(80, 140, 1),
                              min.l=10*log2(76))},
        mc.cores=24)
}


cairo_pdf("common_length_retained_267_276_2_p.pdf", width=0.9 * a4.w * pdf.m, height=0.8 * a4.h * pdf.m, onefile=TRUE)
par(mfrow=c(5,4))
for(i in 1:length(ci.267.2)){
    sp1 <- sp.names[ ex.align.2.k2.nj$tip.label[ leaves.267.i[i] ] ]
    for(j in 1:length(ci.267.2[[i]])){
        sp2 <- sp.names[ ex.align.2.k2.nj$tip.label[ leaves.276$i[j] ] ]
        plot.tested.stat( ci.267.2[[i]][[j]], stat.name='p', sp1, sp2, cex.main=0.7 )
    }
}
dev.off()

cairo_pdf("common_length_retained_267_276_2_r.pdf", width=0.9 * a4.w * pdf.m, height=0.8 * a4.h * pdf.m, onefile=TRUE)
par(mfrow=c(5,4))
for(i in 1:length(ci.267.2)){
    sp1 <- sp.names[ ex.align.2.k2.nj$tip.label[ leaves.267.i[i] ] ]
    for(j in 1:length(ci.267.2[[i]])){
        sp2 <- sp.names[ ex.align.2.k2.nj$tip.label[ leaves.276$i[j] ] ]
        plot.tested.stat( ci.267.2[[i]][[j]], stat.name='r', sp1, sp2, cex.main=0.7 )
    }
}
dev.off()

#### The results from the above are rathe difficult to interpret. We see very clear effects
#### for two short genomes (obs / expected) ratios even when we split down to a single 10th
#### of a log 2. But for long and long, or long and short, we do not see any clear effects,
#### and in fact observe less than would be expected, which I do not know how to explain.
#### We can calculate aggregate p-values for these, and we should get the expected effects.
#### But before that we should visualise the general tendency for change.

## to get a set of plots:
## 266: Common ancestor of teleosts
## 78:  D. rerio
## 3: Betta splendens
## 2: Takifugu rubripes 


cairo_pdf("teleost_ancestor_descendant.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
par(mfrow=c(6, 4))
tel.i <- which(sp.class[,'teleostei'])
tel.i <- tel.i[ order( sp.y[names(tel.i)], decreasing=TRUE ) ]
cols <- hcl.colors(256, "YlOrRd", rev=TRUE)
teleost.is.evo <- sapply(tel.i, function(i){
    tmp <- hist.2d( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    sp <- rownames(sp.class)[i]
    substring(sp, 1, 1) <- toupper( substring(sp, 1, 1) )
    sp <- sub("_", " ", sp)
    with(tmp, plot(xv/10, yv/10, col=rgb(0,0,0,0.05), cex=0.5, xlab="Ancestor", ylab=sp))
    mtext(sp, side=3, line=1)
    abline(0, 1, col='red', lty=2)
    with(tmp, image(x/10, y/10, h, col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, sweep(h, 1, apply(h, 1, max), "/"), col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, t(scale(t(h))), col=cols, xlab="Ancestor", ylab=sp))
    tmp
})
names(teleost.is.evo) <- rownames(sp.class)[tel.i]
dev.off()

cairo_pdf("eutherian_ancestor_descendant.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
par(mfrow=c(6, 4))
euth.i <- which(sp.class[,'eutheria'])
euth.i <- euth.i[ order( sp.y[names(euth.i)], decreasing=TRUE ) ]
cols <- hcl.colors(256, "YlOrRd", rev=TRUE)
eutherian.is.evo <- sapply(euth.i, function(i){
    tmp <- hist.2d( int.s.inf.2[[215]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    sp <- rownames(sp.class)[i]
    substring(sp, 1, 1) <- toupper( substring(sp, 1, 1) )
    sp <- sub("_", " ", sp)
    with(tmp, plot(xv/10, yv/10, col=rgb(0,0,0,0.05), cex=0.5, xlab="Ancestor", ylab=sp))
    mtext(sp, side=3, line=1)
    abline(0, 1, col='red', lty=2)
    with(tmp, image(x/10, y/10, h, col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, sweep(h, 1, apply(h, 1, max), "/"), col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, t(scale(t(h))), col=cols, xlab="Ancestor", ylab=sp))
    tmp
})
names(eutherian.is.evo) <- rownames(sp.class)[euth.i]
dev.off()

cairo_pdf("mammalian_ancestor_descendant.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
par(mfrow=c(6, 4))
mam.i <- which(sp.class[,'mammalia'])
mam.i <- mam.i[ order( sp.y[names(mam.i)], decreasing=TRUE ) ]
cols <- hcl.colors(256, "YlOrRd", rev=TRUE)
mammalian.is.evo <- sapply(mam.i, function(i){
    tmp <- hist.2d( int.s.inf.2[[239]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    sp <- rownames(sp.class)[i]
    substring(sp, 1, 1) <- toupper( substring(sp, 1, 1) )
    sp <- sub("_", " ", sp)
    with(tmp, plot(xv/10, yv/10, col=rgb(0,0,0,0.05), cex=0.5, xlab="Ancestor", ylab=sp))
    mtext(sp, side=3, line=1)
    abline(0, 1, col='red', lty=2)
    with(tmp, image(x/10, y/10, h, col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, sweep(h, 1, apply(h, 1, max), "/"), col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, t(scale(t(h))), col=cols, xlab="Ancestor", ylab=sp))
    tmp
})
names(mammalian.is.evo) <- rownames(sp.class)[mam.i]
dev.off()


cairo_pdf("tetrapod_ancestor_descendant.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
par(mfrow=c(6, 4))
tpod.i <- rev( collect.leaves(248)$i )
cols <- hcl.colors(256, "YlOrRd", rev=TRUE)
tpod.is.evo <- sapply(tpod.i, function(i){
    tmp <- hist.2d( int.s.inf.2[[248]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    sp <- rownames(sp.class)[i]
    substring(sp, 1, 1) <- toupper( substring(sp, 1, 1) )
    sp <- sub("_", " ", sp)
    with(tmp, plot(xv/10, yv/10, col=rgb(0,0,0,0.05), cex=0.5, xlab="Ancestor", ylab=sp))
    mtext(sp, side=3, line=1)
    abline(0, 1, col='red', lty=2)
    with(tmp, image(x/10, y/10, h, col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, sweep(h, 1, apply(h, 1, max), "/"), col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, t(scale(t(h))), col=cols, xlab="Ancestor", ylab=sp))
    tmp
})
names(tpod.is.evo) <- rownames(sp.class)[tpod.i]
dev.off()

cairo_pdf("teleost_290_ancestor_descendant.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
par(mfrow=c(6, 4))
tel.290.i <- rev( collect.leaves(290)$i )
cols <- hcl.colors(256, "YlOrRd", rev=TRUE)
teleost.290.is.evo <- sapply(tel.290.i, function(i){
    tmp <- hist.2d( int.s.inf.2[[290]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    sp <- rownames(sp.class)[i]
    substring(sp, 1, 1) <- toupper( substring(sp, 1, 1) )
    sp <- sub("_", " ", sp)
    with(tmp, plot(xv/10, yv/10, col=rgb(0,0,0,0.05), cex=0.5, xlab="Ancestor", ylab=sp))
    mtext(sp, side=3, line=1)
    abline(0, 1, col='red', lty=2)
    with(tmp, image(x/10, y/10, h, col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, sweep(h, 1, apply(h, 1, max), "/"), col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, t(scale(t(h))), col=cols, xlab="Ancestor", ylab=sp))
    tmp
})
names(teleost.290.is.evo) <- rownames(sp.class)[tel.290.i]
dev.off()

cairo_pdf("teleost_288_ancestor_descendant.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
par(mfrow=c(6, 4))
tel.288.i <- rev( collect.leaves(288)$i )
cols <- hcl.colors(256, "YlOrRd", rev=TRUE)
teleost.288.is.evo <- sapply(tel.288.i, function(i){
    tmp <- hist.2d( int.s.inf.2[[288]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    sp <- rownames(sp.class)[i]
    substring(sp, 1, 1) <- toupper( substring(sp, 1, 1) )
    sp <- sub("_", " ", sp)
    with(tmp, plot(xv/10, yv/10, col=rgb(0,0,0,0.05), cex=0.5, xlab="Ancestor", ylab=sp))
    mtext(sp, side=3, line=1)
    abline(0, 1, col='red', lty=2)
    with(tmp, image(x/10, y/10, h, col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, sweep(h, 1, apply(h, 1, max), "/"), col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, t(scale(t(h))), col=cols, xlab="Ancestor", ylab=sp))
    tmp
})
names(teleost.288.is.evo) <- rownames(sp.class)[tel.288.i]
dev.off()

cairo_pdf("teleost_276_ancestor_descendant.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
par(mfrow=c(6, 4))
tel.276.i <- rev( collect.leaves(276)$i )
cols <- hcl.colors(256, "YlOrRd", rev=TRUE)
teleost.276.is.evo <- sapply(tel.276.i, function(i){
    tmp <- hist.2d( int.s.inf.2[[276]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    sp <- rownames(sp.class)[i]
    substring(sp, 1, 1) <- toupper( substring(sp, 1, 1) )
    sp <- sub("_", " ", sp)
    with(tmp, plot(xv/10, yv/10, col=rgb(0,0,0,0.05), cex=0.5, xlab="Ancestor", ylab=sp))
    mtext(sp, side=3, line=1)
    abline(0, 1, col='red', lty=2)
    with(tmp, image(x/10, y/10, h, col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, sweep(h, 1, apply(h, 1, max), "/"), col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, t(scale(t(h))), col=cols, xlab="Ancestor", ylab=sp))
    tmp
})
names(teleost.276.is.evo) <- rownames(sp.class)[tel.76.i]
dev.off()

cairo_pdf("teleost_267_ancestor_descendant.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
par(mfrow=c(6, 4))
tel.267.i <- rev( collect.leaves(267)$i )
cols <- hcl.colors(256, "YlOrRd", rev=TRUE)
teleost.267.is.evo <- sapply(tel.267.i, function(i){
    tmp <- hist.2d( int.s.inf.2[[267]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    sp <- rownames(sp.class)[i]
    substring(sp, 1, 1) <- toupper( substring(sp, 1, 1) )
    sp <- sub("_", " ", sp)
    with(tmp, plot(xv/10, yv/10, col=rgb(0,0,0,0.05), cex=0.5, xlab="Ancestor", ylab=sp))
    mtext(sp, side=3, line=1)
    abline(0, 1, col='red', lty=2)
    with(tmp, image(x/10, y/10, h, col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, sweep(h, 1, apply(h, 1, max), "/"), col=cols, xlab="Ancestor", ylab=sp))
    with(tmp, image(x/10, y/10, t(scale(t(h))), col=cols, xlab="Ancestor", ylab=sp))
    tmp
})
names(teleost.267.is.evo) <- rownames(sp.class)[tel.276.i]
dev.off()


tmp <- diff.2d(266, 78, inf=int.s.inf.2 )

tmp2 <- hist.2d( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[78]]$state[,1] )
with(tmp2, image(x, y, h))

tmp3 <- summarise.delta( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[78]]$state[,1] )
tmp4 <- summarise.delta( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[3]]$state[,1] )
tmp5 <- summarise.delta( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[27]]$state[,1] )

par(mfrow=c(7,3))
for(i in c(2, 3, 27, 9, 10, 78, 105)){
    tmp <- summarise.delta( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    plot(tmp[,'x'], tmp[,'mode'], main='mode', type='l')
    abline(tmp[1,'x'], -1, col='red')
    plot(tmp[,'x'], tmp[,'median'], main='median', type='l')
    abline(tmp[1,'x'], -1, col='red')
    plot(tmp[,'x'], tmp[,'mean'], main='mean', type='l')
}

## for birds and  mammalians. Start from the tetrapod common ancestor: 248
par(mfrow=c(7,3))
for(i in c(61, 75, 108, 127, 107, 118, 155)){
    tmp <- summarise.delta( int.s.inf.2[[248]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    plot(tmp[,'x'], tmp[,'mode'], main='mode', type='l')
    abline(tmp[1,'x'], -1, col='red')
    plot(tmp[,'x'], tmp[,'median'], main='median', type='l')
    abline(tmp[1,'x'], -1, col='red')
    plot(tmp[,'x'], tmp[,'mean'], main='mean', type='l')
}




##### in order to compare minimisation across the whole clade we can ask
##### what fraction of a given size of introns in an ancestor are minimised
##### in a descendant. This can be done for the whole tree, including intermediate
##### nodes. The summarise.minimisation function does something like this.

## for this it helps if we have the tree and the colours and nodes associated with
## it.
## we have already obtained thesein clade.nodes, and nodes.col
## so we can obtain the required information.

## get descendants of node 259
jawed.nodes <- collect.descendants(259)
b <- is.na(jawed.nodes$names)
jawed.nodes$leaf <- !b
jawed.nodes$names[b] <- as.character( jawed.nodes$i[b] )

jawed.minimisation <- lapply(jawed.nodes$i, function(i){
    summarise.minimisation( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[i]]$state[,1] )
})

### give names to these:
names( jawed.minimisation ) <- jawed.nodes$i

### do the same for:
### 1. All teleosts (descendants of 266)
### 2. Descendants of 276
### 3. Descendants of 267
### 4. All mammals (descendants of 239)

teleost.nodes <- collect.descendants(266)
b <- is.na(teleost.nodes$names)
teleost.nodes$leaf <- !b
teleost.nodes$names[b] <- as.character(teleost.nodes$i[b])

teleost.276.nodes <- collect.descendants(276)
b <- is.na(teleost.276.nodes$names)
teleost.276.nodes$leaf <- !b
teleost.276.nodes$names[b] <- as.character(teleost.276.nodes$i[b])

teleost.267.nodes <- collect.descendants(267)
b <- is.na(teleost.267.nodes$names)
teleost.267.nodes$leaf <- !b
teleost.267.nodes$names[b] <- as.character(teleost.267.nodes$i[b])

mammal.nodes <- collect.descendants(239)
b <- is.na(mammal.nodes$names)
mammal.nodes$leaf <- !b
mammal.nodes$names[b] <- as.character(mammal.nodes$i[b])

eutherian.nodes <- collect.descendants(215)
b <- is.na(eutherian.nodes$names)
eutherian.nodes$leaf <- !b
eutherian.nodes$names[b] <- as.character(eutherian.nodes$i[b])


teleost.minimisation <- lapply( teleost.nodes$i, function(i){
    summarise.minimisation( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[i]]$state[,1] )
})
names(teleost.minimisation) <- teleost.nodes$i

teleost.276.minimisation <- lapply( teleost.276.nodes$i, function(i){
    summarise.minimisation( int.s.inf.2[[276]]$state[,1], int.s.inf.2[[i]]$state[,1] )
})
names(teleost.276.minimisation) <- teleost.276.nodes$i

teleost.267.minimisation <- lapply( teleost.267.nodes$i, function(i){
    summarise.minimisation( int.s.inf.2[[267]]$state[,1], int.s.inf.2[[i]]$state[,1] )
})
names(teleost.267.minimisation) <- teleost.267.nodes$i

mammal.minimisation <- lapply( mammal.nodes$i, function(i){
    summarise.minimisation( int.s.inf.2[[239]]$state[,1], int.s.inf.2[[i]]$state[,1] )
})
names(mammal.minimisation) <- mammal.nodes$i

eutherian.minimisation <- lapply( eutherian.nodes$i, function(i){
    summarise.minimisation( int.s.inf.2[[215]]$state[,1], int.s.inf.2[[i]]$state[,1] )
})
names(eutherian.minimisation) <- eutherian.nodes$i

## we can try to plot the values for all of the nodes:
par(mfrow=c(1,1))
x <- jawed.minimisation[[1]]$x
y <- sapply( jawed.minimisation, function(x){ x$p } )

plot(x, y[[1]], ylim=range(unlist(y)), type='n', ylab='Proportion minimised')
for(i in 1:length(jawed.minimisation)){
    if(!is.na( jawed.nodes$names[i] ))
        with(jawed.minimisation[[i]], lines(x, p, col=nodes.col[ jawed.nodes$i[i] ] ))
}

## use the names to take some examples..

i <- c(1,2,3,4,9,10,28,30,78, 127, 155)
plot(x, y[[1]], ylim=c(0,1), type='n', ylab='Proportion minimised')
for(j in i){
    with(jawed.minimisation[[as.character(j)]], lines(x, p, type='l', col=nodes.col[j]))
}

## for somewhat nicer visualisation we can consider smoothening the curve:
dyn.load("~/R/c/convolve_ks.so")
kernel <- dnorm(0:6, sd=2)

x <- jawed.minimisation[[1]]$x
y <- sapply( jawed.minimisation, function(x){ x$p } )
i <- c(1,2,3,4,9,10,28,30,78, 127, 155)
plot(x, y[[1]], ylim=c(0,1), type='n', ylab='Proportion minimised')
for(j in i){
    with(jawed.minimisation[[as.character(j)]], {
        y <-  .Call("convolve_ks", p, kernel)
        lines(x,
              y,
              type='l', col=nodes.col[j])
        text(x[1], y[1], j, pos=2)
    })
}


vis.jawed.minim <- function(cex.xylab=1, line.xylab=2, cex.leg=1, cex.axis=1, mgp.axis=c(3,1,0), tcl.axis=-mgp.axis[2],
                            mar=c(5.1, 4.1, 2.1, 4.1)){
    i <- 1:172
    i <- i[ i %in% jawed.nodes$i ]
    i <- i[ int.s.b[ ex.align.2.k2.nj$tip.label[i] ]]
    x <- jawed.minimisation[[1]]$x / 10
    h <- with(int.s.inf.2[[259]], hist(state[ state[,1] >= 10*log2(76),1] / 10, breaks=60, plot=FALSE ))
    par(mar=mar)
    plot.new()
    plot.window(xlim=range(h$breaks), ylim=range(h$count), yaxs='i', xaxs='i')
    mtext("Proportion minimised", side=2, line=line.xylab, cex=cex.xylab)
    mtext("Ancestral size", side=1, line=line.xylab, cex=cex.xylab)
    mtext("Ancestral count", side=4, line=line.xylab, cex=cex.xylab)
    with(h, rect( breaks[-length(breaks)], 0, breaks[-1], counts, col=rgb(0.8, 0.8, 0.8), border=NA))
    axis(1, cex.axis=cex.axis, mgp=mgp.axis, tcl=tcl.axis)
    axis(4, cex.axis=cex.axis, mgp=mgp.axis, tcl=tcl.axis)
    plot.window( xlim=range(h$breaks), ylim=c(0,0.7), yaxs='i', xaxs='i' )
    axis(2, cex.axis=cex.axis, mgp=mgp.axis, tcl=tcl.axis)
    abline(v=log2(100), lty=2)
    for(j in i){
        with(jawed.minimisation[[as.character(j)]], {
            y <-  .Call("convolve_ks", p, kernel)
            lines(x/10,
                  y,
                  type='l', col=nodes.col[j])
        })
    }
    ## add lines for teleost, eutherian and avian ancestors
    anc.i <- c(266, 268, 215)
    anc.lab <- c('teleost anc', 'avian anc', 'eutherian anc')
    anc.col <- c("#FF0000FF", "#00F0F0FF", "#F000F0FF")
    for(i in 1:length(anc.i)){
        j <- anc.i[i]
        with(jawed.minimisation[[as.character(j)]], {
            y <-  .Call("convolve_ks", p, kernel)
            lines(x/10,
                  y,
                  type='l', col=anc.col[i], lwd=4, lty=1)
        })
    }
    legend('topright', names(class.col), lwd=2, col=class.col, cex=cex.leg)
}

vis.jawed.minim()



vis.mod.delta <- function(delta=jawed.delta, alpha=0.2, lwd=1, extant.only=TRUE,
                          sd=3, min.f=0.01, smooth.modes=TRUE,
                          cex.leg=1, cex.xylab=1, cex.axis=1, mgp.axis=c(1,0.2,0),
                          tcl.axis=-mgp.axis[2], xlab="", ylab=""){
## show all by default..
    sp.i <- as.numeric(names(delta))
    if(extant.only)
        delta <- delta[ sp.i <= 172 & int.s.b[sp.i] ]
    sp.i <- as.numeric(names(delta))
    cols <- nodes.col[ sp.i ]
    sp <- ex.align.2.k2.nj$tip.label[sp.i]
    all.coord <- t(do.call(cbind,
                         lapply( 1:length(delta), function(i){
                             if(length(delta[[i]])){
                                 do.call(cbind, lapply(1:length(delta[[i]]$ridges$pts), function(j){
                                     with(delta[[i]]$ridges, cbind( rbind(pts[[j]][,1],
                                                                 t(.Call("gs_smooth", pts[[j]], sd, min.f)),
                                                                 i, j ), c(NA, NA, i, j))) }))
                             }else{
                                 matrix(nrow=4, ncol=1)
                             }
                         })))
    all.coord[,1] <- all.coord[,1] / 10
    all.coord[,2] <- all.coord[,2] / 10
    dl <- all.coord[,2] - all.coord[,1]
    pts <- cbind(all.coord[,1], dl)
    plot(pts, type='n', mgp=mgp.axis, cex.lab=cex.xylab, cex.axis=cex.axis, tcl=tcl.axis,
         xlab=xlab, ylab=ylab)
##    plot(all.coord, type='n')
    cols.a <- apply(col2rgb(cols), 2, function(x){ do.call(rgb, as.list(c(x/255, alpha)) )})
    invisible(tapply( 1:nrow(all.coord), all.coord[,3],
                     function(i){
                         print(sp[all.coord[i[1],3]])
                         lines( pts[i,1:2], col=cols.a[ all.coord[i[1],3]], lwd=lwd) }))
    abline(h=0)
}

vis.mod.delta(alpha=0.2, lwd=2)

par(mfrow=c(2,2))
i <- 1:172
i <- i[ i %in% teleost.nodes$i ]
i <- i[ int.s.b[ ex.align.2.k2.nj$tip.label[i] ]]
x <- teleost.minimisation[[1]]$x
plot(x, teleost.minimisation[[1]]$p, ylim=c(0,1), type='n', ylab='Proportion minimised')
abline(v=10*log2(100), lty=2)
for(j in i){
    with(teleost.minimisation[[as.character(j)]], {
        y <-  .Call("convolve_ks", p, kernel)
        lines(x,
              y,
              type='l', col=nodes.col[j])
        text(x[1], y[1], j, pos=2)
    })
##    inpt <- readline(paste(j, ex.align.2.k2.nj$tip.label[j], " "))
}

i <- 1:172
i <- i[ i %in% teleost.276.nodes$i ]
i <- i[ int.s.b[ ex.align.2.k2.nj$tip.label[i] ]]
x <- teleost.276.minimisation[[1]]$x
plot(x, teleost.276.minimisation[[1]]$p, ylim=c(0,1), type='n', ylab='Proportion minimised')
abline(v=10*log2(100), lty=2)
for(j in i){
    with(teleost.276.minimisation[[as.character(j)]], {
        y <-  .Call("convolve_ks", p, kernel)
        lines(x,
              y,
              type='l', col=nodes.col[j])
        text(x[1], y[1], j, pos=2)
    })
    inpt <- readline(paste(j, ex.align.2.k2.nj$tip.label[j], " "))
}

i <- 1:172
i <- i[ i %in% teleost.267.nodes$i ]
i <- i[ int.s.b[ ex.align.2.k2.nj$tip.label[i] ]]
x <- teleost.267.minimisation[[1]]$x
plot(x, teleost.267.minimisation[[1]]$p, ylim=c(0,1), type='n', ylab='Proportion minimised')
abline(v=10*log2(100), lty=2)
for(j in i){
    with(teleost.267.minimisation[[as.character(j)]], {
        y <-  .Call("convolve_ks", p, kernel)
        lines(x,
              y,
              type='l', col=nodes.col[j])
        text(x[1], y[1], j, pos=2)
    })
    inpt <- readline(paste(j, ex.align.2.k2.nj$tip.label[j], " "))
}

## i <- 1:172
## i <- i[ i %in% mammal.nodes$i ]
## i <- i[ int.s.b[ ex.align.2.k2.nj$tip.label[i] ]]
## x <- mammal.minimisation[[1]]$x
## plot(x, mammal.minimisation[[1]]$p, ylim=c(0,1), type='n', ylab='Proportion minimised')
## abline(v=10*log2(100), lty=2)
## for(j in i){
##     with(mammal.minimisation[[as.character(j)]], {
##         y <-  .Call("convolve_ks", p, kernel)
##         lines(x,
##               y,
##               type='l', col=nodes.col[j])
##         text(x[1], y[1], j, pos=2)
##     })
##     inpt <- readline(paste(j, ex.align.2.k2.nj$tip.label[j], " "))
## }

i <- 1:172
i <- i[ i %in% eutherian.nodes$i ]
i <- i[ int.s.b[ ex.align.2.k2.nj$tip.label[i] ]]
x <- eutherian.minimisation[[1]]$x
##plot(x, eutherian.minimisation[[1]]$p, ylim=c(0,1), type='n', ylab='Proportion minimised')
abline(v=10*log2(100), lty=2)
for(j in i){
    with(eutherian.minimisation[[as.character(j)]], {
        y <-  .Call("convolve_ks", p, kernel)
        lines(x,
              y,
              type='l', col=nodes.col[j])
        text(x[1], y[1], j, pos=2)
    })
##    inpt <- readline(paste(j, ex.align.2.k2.nj$tip.label[j], " "))
}


## genome size in genome.size
## int.s.b, selects species with fewer very short introns.

## sum of introns larger than 128 bp but smaller than 2^10
## in the jawed vertebrate ancestor that have been minimised
## in the descendant (< 100)
jawed.minim.10 <- t(sapply(jawed.minimisation, function(x){
    row.b <- x$x > 70 & x$x <= 100
    col.b <- x$y <= 10*log2(100)
    c('t'=sum(x$tbl[row.b,]), 'm'=sum(x$tbl[row.b,col.b]), 'p'=sum(x$tbl[row.b,col.b])/sum(x$tbl[row.b,]))
}))

### then lets get some points from these.
vis.min.summary <- function(b=int.s.b, id.pts=NULL, cex.lab=1, cex.leg=1,
                            cex.xylab=1, mgp.axis=c(3,1,0), tcl.axis=-mgp.axis[2], cex.axis=1,
                            do.id=FALSE){
    i <- which( jawed.nodes$i <= 172 & int.s.b[ jawed.nodes$name ] )
    x <- log2(genome.size[ jawed.nodes$name[i] ])
    y <- jawed.minim.10[i, 'p']
##
    plot( x, y, col=nodes.col[ jawed.nodes$i[i] ], pch=19, xlab='log2 genome size', ylab='Proportion minimised',
         mgp=mgp.axis, tcl=tcl.axis, cex.lab=cex.xylab, cex.axis=cex.axis)
    legend('topright', names(class.col), pch=19, col=class.col, cex=cex.leg)
    if(is.null(id.pts) && do.id){
        id.pts <- identify( x, y, sp.labels[jawed.nodes$name[i]], pos=TRUE )
    }else{
        if(!is.null(id.pts))
            text(x[id.pts$ind], y[id.pts$ind], labels=sp.labels[(jawed.nodes$name[i])[id.pts$ind]], pos=id.pts$pos, cex=cex.lab )
    }
    invisible(id.pts)
}

pts <- vis.min.summary()

vis.min.summary(id.pts=pts)

plot( log10(x), y, col=nodes.col[ jawed.nodes$i[i] ], pch=19 )
identify( log10(x), y, jawed.nodes$name[i] )
    

i <- which( jawed.nodes$i <= 172 & int.s.b[ jawed.nodes$name ] )
x <- genome.size[ jawed.nodes$name[i] ]
y <- jawed.minim.10[i, 'p']
plot( x, y, col=nodes.col[ jawed.nodes$i[i] ], pch=19 )
identify( x, y, jawed.nodes$name[i] )

plot( log10(x), y, col=nodes.col[ jawed.nodes$i[i] ], pch=19 )
identify( log10(x), y, jawed.nodes$name[i] )


jawed.minim.10.12 <- t(sapply(jawed.minimisation, function(x){
    row.b <- x$x > 100 & x$x <= 120
    col.b <- x$y <= 10*log2(100)
    c('t'=sum(x$tbl[row.b,]), 'm'=sum(x$tbl[row.b,col.b]), 'p'=sum(x$tbl[row.b,col.b])/sum(x$tbl[row.b,]))
}))

y2 <- jawed.minim.10.12[i, 'p']
points( x, y2, col=nodes.col[ jawed.nodes$i[i] ], pch=1 )
identify( x, y2, jawed.nodes$name[i] )

### So there is a very clear difference indeed.

### two dimensional distributions with marginal distributions
### using split.screen()

d.scale <- function(m){ t(scale(t(m))) }

dist.w=0.25
dist.h=0.25
sc <- split.screen( c(1,4) )
hist.2d.2( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[78]]$state[,1], 1, dist.w=dist.w, dist.h=dist.h, main="Danio rerio", h.tform=d.scale )
hist.2d.2( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[2]]$state[,1], 2, dist.w=dist.w, dist.h=dist.h, main="Takifugu rubripes", h.tform=d.scale)
hist.2d.2( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[3]]$state[,1], 3, dist.w=dist.w, dist.h=dist.h, main="Betta splendens", h.tform=d.scale )
hist.2d.2( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[10]]$state[,1], 4, dist.w=dist.w, dist.h=dist.h, main="Gasterosteus aculeatus", h.tform=d.scale )
close.screen(all.screens=TRUE)


d.scale <- function(m){ sweep(m, 1, apply(m, 1, max), "/") }

dist.w=0.25
dist.h=0.25
sc <- split.screen( c(2,4) )
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[78]]$state[,1], 1, dist.w=dist.w, dist.h=dist.h, main="Danio rerio", h.tform=d.scale,
          image.after=abline(0,1))
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[2]]$state[,1], 2, dist.w=dist.w, dist.h=dist.h, main="Takifugu rubripes", h.tform=d.scale,
          image.after=abline(0,1))
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[56]]$state[,1], 3, dist.w=dist.w, dist.h=dist.h, main="F. heteroclitus", h.tform=d.scale,
          image.after=abline(0,1) )
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[10]]$state[,1], 4, dist.w=dist.w, dist.h=dist.h, main="D. clupeoides", h.tform=d.scale,
          image.after=abline(0,1) )
## and some non teleosts.. 
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[127]]$state[,1], 5, dist.w=dist.w, dist.h=dist.h, main="Mus m", h.tform=d.scale,
          image.after=abline(0,1) )
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[155]]$state[,1], 6, dist.w=dist.w, dist.h=dist.h, main="Homo s", h.tform=d.scale,
          image.after=abline(0,1))
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[75]]$state[,1], 7, dist.w=dist.w, dist.h=dist.h, main="T guttata", h.tform=d.scale,
          image.after=abline(0,1) )
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[61]]$state[,1], 8, dist.w=dist.w, dist.h=dist.h, main="G gallus", h.tform=d.scale,
          image.after=abline(0,1) )
close.screen(all.screens=TRUE)

dist.w=0.25
dist.h=0.25
sc <- split.screen( c(2,4) )
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[78]]$state[,1], 1, dist.w=dist.w, dist.h=dist.h, main="Danio rerio", h.tform=d.scale,
          image.after=abline(0,1))
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[2]]$state[,1], 2, dist.w=dist.w, dist.h=dist.h, main="Takifugu rubripes", h.tform=d.scale,
          image.after=abline(0,1))
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[10]]$state[,1], 3, dist.w=dist.w, dist.h=dist.h, main="D. clupeoides", h.tform=d.scale,
          image.after=abline(0,1) )
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[49]]$state[,1], 4, dist.w=dist.w, dist.h=dist.h, main="L. oculatus", h.tform=d.scale,
          image.after=abline(0,1) )
## and some non teleosts.. 
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[169]]$state[,1], 5, dist.w=dist.w, dist.h=dist.h, main="E. calabricus", h.tform=d.scale,
          image.after=abline(0,1) )
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[53]]$state[,1], 6, dist.w=dist.w, dist.h=dist.h, main="C. mili", h.tform=d.scale,
          image.after=abline(0,1))
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[134]]$state[,1], 7, dist.w=dist.w, dist.h=dist.h, main="L. chalumnae", h.tform=d.scale,
          image.after=abline(0,1) )
hist.2d.2( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[93]]$state[,1], 8, dist.w=dist.w, dist.h=dist.h, main="C. abingdonii", h.tform=d.scale,
          image.after=abline(0,1) )
close.screen(all.screens=TRUE)

### Does the variance between ancestor and descendant increase or decrease
### with ancestral size?

tmp <- with(jawed.minimisation[['2']], table.to.var(tbl, x, y))
plot(tmp$rv, sqrt(tmp$v), cex=log(tmp$n)/log(4))

tmp <- with(jawed.minimisation[['78']], table.to.var(tbl, x, y))
plot(tmp$rv, sqrt(tmp$v), cex=log(tmp$n)/log(4))

tmp2 <- with(jawed.minimisation[['78']], table.to.var.2(tbl, x, y, 200))
plot(tmp2[,'rv'], sqrt(tmp2[,'v']), cex=0.5 * log(tmp2[,'n']/log(10)), type='b')
### this turns out to be difficult to interpret.

### these suggest that calculating a statistic capturing the difference
### between the total extant marginal distribution and for ancestral
### size slices might be useful.
## we can recapture a given set of discretised values using un.table
## and then use ks.test. This will not give a p-value, as we will
## have lots of ties, but we can use it to obtain the D statistic
## the bigger, the more difference..

sp.i <- c(2, 3, 9, 10, 49, 169, 53, 134, 61, 75, 127, 155)
par(mfrow=c(3,4))
for(i in sp.i){
    tmp <- with(jawed.minimisation[[as.character(i)]], table.to.ks(tbl, x, y, 1000))
    plot(tmp[,'rv'], tmp[,'d'], main=sp.names[ ex.align.2.k2.nj$tip.label[i] ])
}


## the following plays around with lowess lines a bit. Turns out to be not terribly useful.
## p
sp.i <- c(2, 3, 9, 10, 49, 169, 53, 134, 61, 75, 127, 155)
par(mfrow=c(3,4))
for(i in sp.i){
    tmp <- hist.2d( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    with(tmp, image(x, y, scale(h), col=hcl.colors(256, "YlOrRd", rev=TRUE), main=sp.names[ex.align.2.k2.nj$tip.label[i]]))
    abline(0,1, lty=2)
##    tmp2 <- scale(tmp$h)
    tmp2 <- sweep(tmp$h, 2, apply(tmp$h, 2, max), "/")
##    tmp2 <- tmp2 - min(tmp2, na.rm=TRUE)
    tmp2[is.na(tmp2)] <- 0
    tmp2 <- tmp2 * (sum(tmp$h) / sum(tmp2))

    un.t <- lapply( 1:nrow(tmp2), function(j){
        y <- un.table(tmp2[j,], tmp$y)
        x <- rep( tmp$x[j], length(y) )
        list(x=x, y=y)
    })

    un.t.x <- unlist( sapply(un.t, function(x){ x$x }))
    un.t.y <- unlist( sapply(un.t, function(x){ x$y }))
    lw <- lowess(un.t.y, un.t.x)
    plot(un.t.y, un.t.x, cex=0.24, col=rgb(0,0,0,0.2))
    lines(lw)
}

### How about finding ridges within the data?
### is that at all useful?

with(jawed.minimisation[["78"]], image(x, y, t(scale(t(tbl)))))
with(jawed.minimisation[["78"]], image(x, y, scale(tbl)))

kernel <- dnorm(0:20, sd=6)

with(jawed.minimisation[["78"]], image( x, y, blur.matrix(tbl, kernel) ))
with(jawed.minimisation[["78"]], image( x, y, blur.matrix(t(scale(t(tbl))), kernel), col=im.cols ))
with(jawed.minimisation[["2"]], image( x, y, blur.matrix(t(scale(t(tbl))), kernel), col=cols ))

sp.i <- c(2, 3, 9, 78, 49, 169, 53, 134, 61, 75, 127, 155)
par(mfrow=c(3,4))
for(i in sp.i){
    with(jawed.minimisation[[as.character(i)]],
         image( x, y, blur.matrix(t(scale(t(tbl))), kernel), col=cols,
               main=sp.names[ex.align.2.k2.nj$tip.label[i]]))
    abline(0,1,col='blue')
}

dyn.load("~/R/c/gs_smooth.so")

## We can trace the modes of the 2-dimensional distributions to see how these compare
sp.i <- c(2, 3, 9, 78, 49, 169, 53, 134, 61, 75, 127, 155)
par(mfrow=c(3,4))
for(i in sp.i){
    tmp <- hist.2d( int.s.inf.2[[259]]$state[,1], int.s.inf.2[[i]]$state[,1] )
    with(tmp, image( x, y, blur.matrix(t(scale(t(h))), kernel), col=cols,
                    main=sp.names[ex.align.2.k2.nj$tip.label[i]]))
    tmp.wd <- with(tmp, windowed.density(xv, yv, 1000, 50, max.d=0.1*diff(range(y))))
    wsum <- sapply(tmp.wd$ridges, function(x){
        sum(tmp$h[ 1 + cbind( as.integer(x[,1])-tmp$x[1], as.integer(x[,2])-tmp$y[1] ) ])
    })
    o <- order(wsum, decreasing=TRUE)
    sd <- 3
    min.f <- 0.01
    for(j in 1:4){
        lines( tmp.wd$ridges[[o[j]]][,1], .Call("gs_smooth", tmp.wd$ridges[[o[j]]], sd, min.f), col=j, lwd=3)
    }
}

sd <- 3
min.f <- 0.01
navy <- rgb(0.6, 0.6, 1)

h.tform <- function(x, k=kernel){
    blur.matrix(t(scale(t(x))), k)
}

h.log.blur <- function(x, k=kernel){
    blur.matrix(log(x + 1), k)
}

sp.i <- jawed.nodes$i[ !is.na(jawed.nodes$names) ]
names(sp.i) <- jawed.nodes$names[ !is.na(jawed.nodes$names) ]
## these will be in the correct order as well..
cairo_pdf("jawed_intron_size_change.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
jawed.delta <- vector(mode='list', length=length(sp.i))
names(jawed.delta) <- as.character(rev(sp.i))
##par(mfrow=c(6,4))
count <- 0
for(h in 1:length(sp.i)){
    i <- rev(sp.i)[h]
    i.n <- names( rev(sp.i) )[h]
    if(count %% 24 == 0){
        close.screen(all.screens=TRUE)
        split.screen( c(6,4) )
    }
    jawed.delta[[h]] <- hist.2d.2(int.s.inf.2[[259]]$state[,1], int.s.inf.2[[i]]$state[,1], 1 + count %% 24,
                                  dist.w=0.25, dist.h=0.3, h.tform=h.tform, mar=c(2.6,2.6,1.1,1.1), ax.cex=0.5, ax.mgp=c(3,0.2,0),
                                  xlab="259", main=ifelse( is.na(sp.names[i.n]), i.n, sp.names[i.n] ),
                                  trace.modes=TRUE, x.div=10, y.div=10)
    count <- count + 1
}
dev.off()

sp.i <- teleost.nodes$i[ !is.na(teleost.nodes$names) ]
names(sp.i) <- teleost.nodes$names[ !is.na(teleost.nodes$names) ]
## these will be in the correct order as well..
cairo_pdf("teleost_intron_size_change.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
##par(mfrow=c(6,4))
count <- 0
for(h in 1:length(sp.i)){
    i <- rev(sp.i)[h]
    i.n <- names(rev(sp.i))[h]
    if(count %% 24 == 0){
        close.screen(all.screens=TRUE)
        split.screen( c(6,4) )
    }
    tmp <- hist.2d.2(int.s.inf.2[[266]]$state[,1], int.s.inf.2[[i]]$state[,1], 1 + count %% 24,
                     dist.w=0.25, dist.h=0.3, h.tform=h.tform, mar=c(2.6,2.6,1.1,1.1), ax.cex=0.5, ax.mgp=c(3,0.2,0),
                     xlab="266", main=ifelse( is.na(sp.names[i.n]), i.n, sp.names[i.n] ),
                     trace.modes=TRUE, x.div=10, y.div=10)
    count <- count + 1
##     tmp <- hist.2d( int.s.inf.2[[266]]$state[,1], int.s.inf.2[[i]]$state[,1] )
##     tmp.bl <- with(tmp, blur.matrix(t(scale(t(h))), kernel))
##     tmp.pk <- apply( tmp.bl, 1, function(x){ get.peaks( tmp$y, x, 3, 0.5) })
##     with(tmp, image(x/10, y/10, tmp.bl, col=cols,
##                     main=sp.names[ex.align.2.k2.nj$tip.label[i]],
##                     xlab='Ancestor', ylab='Descendant'))
##     abline(0,1,col=navy, lty=2, lwd=2)
##     b.i <- which(sapply(tmp.pk, function(x){ length(x) > 0 }))
##     pts <- with(tmp, lapply(b.i, function(i){
##         cbind( x[i], tmp.pk[[i]])}))
##     ridges <- trace.lines(pts, diff(range(tmp$y))*0.1)
##     ridge.w <- sapply( ridges, function(x){ sum( tmp.bl[ cbind(1+x[,1]-tmp$x[1], 1+x[,2]-tmp$y[1]) ])})
##     b <- (ridge.w / max(ridge.w)) > 0.05
##     for(i in which(b)){
## ##        lines(ridges[[i]], col='blue', lwd=3)
##         lines(ridges[[i]][,1]/10, .Call("gs_smooth", ridges[[i]], sd, min.f)/10, col=navy, lwd=3)
##     }
}
dev.off()


sd <- 3
min.f <- 0.01
navy <- rgb(0.6, 0.6, 1)
### we can also do the same thing but look at specific pairs of species:
sp.i <- c(2, 10, 26, 76, 78, 61, 127, 155)
## I have redone this with marginal distributions below.. 
## cairo_pdf("intron_size_example_pairs.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
## par(mfrow=c(6,4))
## ##par(mfrow=c(4,7))
## for(i in 1:length(sp.i)){
##     if(i == length(sp.i))
##         break
##     for(j in (i+1):length(sp.i)){
##         tmp.sp <- c(sp.i[j], sp.i[i])
##         for(m in 1:2){
##             tmp.sp <- rev(tmp.sp)
##             ## print(paste(i, ",", j, ":", sp.i[i], "vs", sp.i[j], ":",
##             ##             sp.names[ex.align.2.k2.nj$tip.label[tmp.sp[1]], sp.names[ex.align.2.k2.nj$tip.label[tmp.sp[2]]
##             tmp <- hist.2d( int.s.inf.2[[tmp.sp[1]]]$state[,1], int.s.inf.2[[tmp.sp[2]]]$state[,1] )
##             tmp.bl <- with(tmp, blur.matrix(t(scale(t(h))), kernel))
##             tmp.pk <- apply( tmp.bl, 1, function(x){ get.peaks( tmp$y, x, 3, 0.5) })
##             with(tmp, image(x/10, y/10, tmp.bl, col=cols,
##                             xlab=sp.names[ex.align.2.k2.nj$tip.label[tmp.sp[1]]],
##                             ylab=sp.names[ex.align.2.k2.nj$tip.label[tmp.sp[2]]]))
##             abline(0,1,col=navy, lty=2, lwd=2)
##             b.i <- which(sapply(tmp.pk, function(x){ length(x) > 0 }))
##             pts <- with(tmp, lapply(b.i, function(i){
##                 cbind( x[i], tmp.pk[[i]])}))
##             ridges <- trace.lines(pts, diff(range(tmp$y))*0.1)
##             ridge.w <- sapply( ridges, function(x){ sum( tmp.bl[ cbind(1+x[,1]-tmp$x[1], 1+x[,2]-tmp$y[1]) ])})
##             b <- (ridge.w / max(ridge.w)) > 0.05
##             for(k in which(b)){
##                 ##        lines(ridges[[i]], col='blue', lwd=3)
##                 lines(ridges[[k]][,1]/10, .Call("gs_smooth", ridges[[k]], sd, min.f)/10, col=navy, lwd=3)
##             }
##         }
##     }
## }
## dev.off()

    


x11(xpos=6000)
split.screen( c(2,2) )

hist.2d.2( int.s.inf.2[[sp.i[1]]]$state[,1], int.s.inf.2[[sp.i[2]]]$state[,1], 1, dist.h=0.2, h.tform=h.tform)
hist.2d.2( int.s.inf.2[[sp.i[2]]]$state[,1], int.s.inf.2[[sp.i[1]]]$state[,1], 2, dist.h=0.2, h.tform=h.tform)

hist.2d.2( int.s.inf.2[[sp.i[1]]]$state[,1], int.s.inf.2[[sp.i[5]]]$state[,1], 3, dist.h=0.2, h.tform=h.tform)
hist.2d.2( int.s.inf.2[[sp.i[5]]]$state[,1], int.s.inf.2[[sp.i[1]]]$state[,1], 4, dist.h=0.2, h.tform=h.tform)

sp.i <- c(49, 2, 10, 26, 78, 77, 169, 61, 127, 155)
cairo_pdf("intron_size_example_pairs.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
dists <- vector(mode='list', length=sum( 2:length(sp.i) - 1 ) * 2)
dists.i <- matrix(nrow=length(sp.i), ncol=length(sp.i))
count <- 0
for(i in 1:length(sp.i)){
    if(i == length(sp.i))
        break
    for(j in (i+1):length(sp.i)){
        tmp.sp <- c(sp.i[j], sp.i[i])
        dists.i[i,j] <- count + 1
        dists.i[j,i] <- count + 2
        for(m in 1:2){
            tmp.sp <- rev(tmp.sp)
            if(count %% 24 == 0){
                close.screen(all.screens=TRUE)
                split.screen( c(6,4) )
            }
            screen( 1 + count %% 24 )
            dists[[count+1]] <- hist.2d.2( int.s.inf.2[[tmp.sp[1]]]$state[,1], int.s.inf.2[[tmp.sp[2]]]$state[,1], 1 + count %% 24,
                                         dist.w=0.25, dist.h=0.3,
                                         h.tform=h.tform, mar=c(2.6,2.6,1.1,1.1), ax.cex=0.5, ax.mgp=c(3,0.2,0),
                                         xlab=sp.names[ex.align.2.k2.nj$tip.label[tmp.sp[1]]],
                                         ylab=sp.names[ex.align.2.k2.nj$tip.label[tmp.sp[2]]],
                                         trace.modes=TRUE, x.div=10, y.div=10
                                         )
            dists[[count+1]]$sp.i <- tmp.sp
            dists[[count+1]]$sp <- sp.names[ ex.align.2.k2.nj$tip.label[ tmp.sp ]]
            count <- count + 1
        }
    }
}
dev.off()

## The following plots can be used to indicate the proportion of non-overlapping
## minimisation (i.e. different in different species). But it is probably better
## to concentrate on the teleosts for this as the evolutionary distances get a bit
## large otherwise.
par(mfrow=c(3,4))
for(i in 1:nrow(dists.i)){
    j <- dists.i[i, !is.na(dists.i[i,])]
    xy <- lapply( j, function(k){
        list(x=dists[[k]]$h$x, y=dists[[k]]$h$y)
    })
    ## plot.new()
    ## plot.window( xlim=range(unlist(lapply(xy, function(x){ x$x }))),
    ##              ylim=range(unlist(lapply(xy, function(x){ x$y }))))
    ## axis(1)
    ## axis(2)
    min.lines <- lapply( j, function(k){
        with(dists[[k]]$h, {
             x <- sweep( h, 1, rowSums(h), "/" )
             x[ !is.finite(x) ] <- 0
             x <- rowSums(x[,1:5])
             x
        })
##        dists[[j]]$ht[,1]
        ##lines( x$pts[[1]] )
##        tmp <- sapply( x$pts[x$b], lines )
    })
    sp <- lapply(dists[j], function(x){
        x$sp
    })
    plot(xy[[1]]$x, min.lines[[1]], ylim=range(unlist(min.lines)), type='n',
         xlab=sp[[1]][1])
    for(j in 1:length(min.lines)){
        lines(xy[[j]]$x, min.lines[[j]], col=j, lwd=2)
    }
    legend('topright', legend=sapply(sp, function(x){x[2]}), col=1:length(sp), lty=1, lwd=2)
}

### We can also look at the proportion of pairs of sizes that are within a given range
sp.i <- c(49, 2, 10, 26, 78, 77, 169, 61, 127, 155)
is.cons <- vector(mode='list', length=sum( 2:length(sp.i) - 1 ) * 2)
is.cons.i <- matrix(nrow=length(sp.i), ncol=length(sp.i))
colnames(is.cons.i) <- sp.i
rownames(is.cons.i) <- sp.i
count <- 0
for(i in 1:length(sp.i)){
    if(i == length(sp.i))
        break
    for(j in (i+1):length(sp.i)){
        tmp.sp <- c(sp.i[c(j,i)])
        is.cons.i[i,j] <- count + 1
        is.cons.i[j,i] <- count + 2
        for(m in 1:2){
            tmp.sp <- rev(tmp.sp)
            ## use discretized and log transformed data.. 
            A <- int.s.inf.2[[tmp.sp[1]]]$state[,1]
            B <- int.s.inf.2[[tmp.sp[2]]]$state[,1]
            min.v <- 10 * log2(76)
            b <- A >= min.v & B >= min.v
            A <- A[b]
            B <- B[b]
            ## this will make for a noisy signal, but it can be cleaned up later
            count <- count + 1
            is.cons[[count]] <- list(sp.i=tmp.sp,
                                     sp=sp.names[ ex.align.2.k2.nj$tip.label[tmp.sp] ],
                                     dt=tapply(1:length(A), A, function(k){
                                         table( abs( A[k] - B[k] ) )}) )
        }
    }
}

i <- 33
y <- with( is.cons[[i]], sapply(dt, function(x){
    d <- as.numeric(names(x))
    sum( x[d < 10]) / sum(x)
}))

### We can also look at the proportion of pairs of sizes that are within a given range
sp.i <- c(259, 266, 49, 1:172)
is.corr <- matrix(nrow=sum( 2:length(sp.i) - 1 ), ncol=9)
is.corr.i <- matrix(nrow=length(sp.i), ncol=length(sp.i))
sp.p <- matrix(nrow=sum( 2:length(sp.i) - 1 ), ncol=2)
colnames(is.corr.i) <- sp.i
rownames(is.corr.i) <- sp.i
count <- 0
for(i in 1:length(sp.i)){
    if(i == length(sp.i))
        break
    for(j in (i+1):length(sp.i)){
        tmp.sp <- c(sp.i[c(i,j)])
        is.corr.i[i,j] <- count + 1
        is.corr.i[j,i] <- count + 1
##        for(m in 1:2){
##        tmp.sp <- rev(tmp.sp)
        ## use discretized and log transformed data.. 
        A <- int.s.inf.2[[tmp.sp[1]]]$state[,1]
        B <- int.s.inf.2[[tmp.sp[2]]]$state[,1]
        min.v <- 10 * log2(76)
        b <- A >= min.v & B >= min.v
        A <- A[b]
        B <- B[b]
        thr <- 95
        b1 <- A <= thr ## & B <= thr
        b2 <- B <= thr ## A > thr & B > thr
        ## this will make for a noisy signal, but it can be cleaned up later
        count <- count + 1
        is.corr[count,] <- c( cor( A[b1], B[b1] ), cor(A[!b1], B[!b1]),
                             cor( A[b2], B[b2] ), cor(A[!b2], B[!b2] ),
                             cor( A[b1 & b2], B[b1 & b2] ),
                             cor( A[!(b1 | b2)], B[!(b1 | b2)] ), cor(A, B),
                             tmp.sp)
        sp.p[count, ] <- ifelse( is.na(ex.align.2.k2.nj$tip.label[ tmp.sp ]),
                                as.character(tmp.sp), sp.labels[ex.align.2.k2.nj$tip.label[ tmp.sp ]] )
        ##        }
    }
}

## this is for 259 against all
r <- 1:174
##r <- which(is.corr[,8] == 266)
d.266 <- sapply( is.corr[r,9], nodes.sep, B=266, tree=ex.align.2.k2.nj )
d.259 <- sapply( is.corr[r,9], nodes.sep, B=259, tree=ex.align.2.k2.nj )

plot( is.corr[ r, 1 ], is.corr[ r, 2 ], col=nodes.col[is.corr[r,9]], pch=19 )
abline(0,1)
identify( is.corr[ r, 1 ], is.corr[ r, 2 ], labels=sp.p[r,2] )

plot( is.corr[ r, 3 ], is.corr[ r, 4 ], col=nodes.col[is.corr[r,9]], pch=19 )
abline(0,1)
identify( is.corr[ r, 3 ], is.corr[ r, 4 ], labels=sp.p[r,2] )

plot( is.corr[ r, 1 ], is.corr[ r, 2 ], col=nodes.col[is.corr[r,9]], pch=19, type='n' )
##text( is.corr[ r, 5 ], is.corr[ r, 6 ], d.266)
text( is.corr[ r, 1 ], is.corr[ r, 2 ], d.259, col=nodes.col[is.corr[r,9]])
points( is.corr[ r, 1 ], is.corr[ r, 2 ], col=nodes.col[is.corr[r,9]], pch=1, cex=3, lwd=3 )
abline(0,1)

plot( is.corr[ r, 3 ], is.corr[ r, 4 ], col=nodes.col[is.corr[r,9]], pch=19, type='n' )
##text( is.corr[ r, 5 ], is.corr[ r, 6 ], d.266)
text( is.corr[ r, 3 ], is.corr[ r, 4 ], d.259, col=nodes.col[is.corr[r,9]])
points( is.corr[ r, 3 ], is.corr[ r, 4 ], col=nodes.col[is.corr[r,9]], pch=1, cex=3, lwd=3 )
abline(0,1)
identify( is.corr[ r, 3 ], is.corr[ r, 4 ], labels=sp.p[r,2] )

plot( is.corr[ r, 3 ], is.corr[ r, 7 ], col=nodes.col[is.corr[r,9]], pch=19, type='n' )
##text( is.corr[ r, 5 ], is.corr[ r, 6 ], d.266)
text( is.corr[ r, 3 ], is.corr[ r, 7 ], d.259, col=nodes.col[is.corr[r,9]])
points( is.corr[ r, 3 ], is.corr[ r, 7 ], col=nodes.col[is.corr[r,9]], pch=1, cex=3, lwd=3 )
abline(0,1)

plot( is.corr[ r, 5 ], is.corr[ r, 6 ], col=nodes.col[is.corr[r,9]], pch=19, type='n' )
##text( is.corr[ r, 5 ], is.corr[ r, 6 ], d.266)
text( is.corr[ r, 5 ], is.corr[ r, 6 ], d.259)
points( is.corr[ r, 5 ], is.corr[ r, 6 ], col=nodes.col[is.corr[r,9]], pch=1, cex=3, lwd=3 )
abline(0,1)

identify( is.corr[ r, 5 ], is.corr[ r, 6 ], labels=sp.p[r,2] )


plot(r, is.corr[r,7], col=nodes.col[is.corr[r,9]], pch=19 )
identify( r, is.corr[ r, 7 ], labels=sp.p[r,2] )

plot(r, is.corr[r,1], col=nodes.col[is.corr[r,9]], pch=19 )
identify( r, is.corr[ r, 1 ], labels=sp.p[r,2] )

par(mfrow=c(2,2))
plot(is.corr[,1], is.corr[,2], col=nodes.col[is.corr[,8]], pch=1, cex=2.2, lwd=3)
points(is.corr[,1], is.corr[,2], col=nodes.col[is.corr[,9]], pch=19, cex=1.2)
abline(0,1)
identify( is.corr[,1], is.corr[,2], labels=paste(sp.p[,1], sp.p[,2]), cex=0.75 )

plot(is.corr[,3], is.corr[,4], col=nodes.col[is.corr[,8]], pch=1, cex=2.2, lwd=3)
points(is.corr[,3], is.corr[,4], col=nodes.col[is.corr[,9]], pch=19, cex=1.2)
abline(0,1)
identify( is.corr[,3], is.corr[,4], labels=paste(sp.p[,1], sp.p[,2]), cex=0.75 )

plot(is.corr[,5], is.corr[,6], col=nodes.col[is.corr[,8]], pch=1, cex=2.2, lwd=3)
points(is.corr[,5], is.corr[,6], col=nodes.col[is.corr[,9]], pch=19, cex=1.2)
abline(0,1)
identify( is.corr[,5], is.corr[,6], labels=paste(sp.p[,1], sp.p[,2]), cex=0.75 )

plot(1:nrow(is.corr), is.corr[,7], col=nodes.col[is.corr[,8]], pch=1, cex=2.2, lwd=3)
points(1:nrow(is.corr), is.corr[,7], col=nodes.col[is.corr[,9]], pch=19, cex=1.2)
identify(1:nrow(is.corr), is.corr[,7], labels=paste(sp.p[,1], sp.p[,2]), cex=0.75 )


i <- 33
y <- with( is.corr[[i]], sapply(dt, function(x){
    d <- as.numeric(names(x))
    sum( x[d < 10]) / sum(x)
}))


plot(as.numeric(names(y)), y, xlab=is.cons[[i]]$sp[1], main=is.cons[[i]]$sp[2])

sp.i <- c(49, 2, 10, 26, 78, 77, 169, 61, 127, 155)
ks.d <- vector(mode='list', length=sum( 2:length(sp.i) - 1 ) * 2)
count <- 0
par(mfrow=c(6,4))
for(i in 1:length(sp.i)){
    if(i == length(sp.i))
        break
    for(j in (i+1):length(sp.i)){
        tmp.sp <- c(sp.i[j], sp.i[i])
        for(m in 1:2){
            tmp.sp <- rev(tmp.sp)
            sp.1 <- ex.align.2.k2.nj$tip.label[tmp.sp[1]]
            sp.2 <- ex.align.2.k2.nj$tip.label[tmp.sp[2]]
            count <- count + 1
            ks.d[[count]] <- with(orth, column.ks( log2(l[,sp.1]), log2(l[,sp.2]), 10, rand.n=100 ))
            ks.d[[count]]$sp.i <- tmp.sp
            ks.d[[count]]$sp <- sp.names[ex.align.2.k2.nj$tip.label[tmp.sp]]
        }
    }
}


cairo_pdf("intron_size_example_pairs_ks.pdf", width=0.95*a4.w*pdf.m, height=0.95*a4.h*pdf.m, onefile=TRUE)
par(mfrow=c(6,4))
for(ks in ks.d){
    with(ks, {
        plot.new()
        plot.window(xlim=range(breaks), ylim=range(n))
        with(par(), rect(usr[1], usr[3], usr[2], usr[4]))
        mtext(sp[1], 1, line=3, cex=0.75)
        mtext(paste(sp[2], 'KS D' ), 2, line=2.5, cex=0.75)
        rect(breaks[-length(breaks)], 0, breaks[-1], n, col=rgb(0.9, 0.9, 0.9), border='white', lwd=0.1)
        axis(4, col.axis=rgb(0.4, 0.5, 0.5))
        plot.window(xlim=range(breaks), ylim=c(0, max(D)))
        dr <- apply( rD, 2, function(x){ quantile(x, probs=c(0.05, 0.5, 0.95)) })
        ##        arrows( mids, dr['5%',], mids, dr['95%',], angle=90, code=3, length=0.05, col=rgb(0,0,0,0.25) )
        bw <- diff(breaks)[1]
        segments( mids-0.25*bw, dr['5%',], mids+0.25*bw, dr['5%',], col=rgb(0,0,0,0.25) )
        segments( mids-0.2*bw, dr['50%',], mids+0.2*bw, dr['50%',], col=rgb(0,0,0,0.25) )
        segments( mids-0.25*bw, dr['95%',], mids+0.25*bw, dr['95%',], col=rgb(0,0,0,0.25) )
        segments( mids, dr['5%',], mids, dr['95%',], col=rgb(0,0,0,0.25) )
        points(mids, D, type='b')
        axis(2)
        axis(1)
    })
}
dev.off()


## calculate mutual information of the full and sampled distribution
## this doesn't actually work. So ignore the following.
sp.i <- c(259, 266, 49, 2, 10, 26, 78, 77, 253, 251, 249, 244, 169, 134, 61, 127, 155)
pair.d <- vector(mode='list', length=sum( 2:length(sp.i) - 1 ) * 2)
pair.d.i <- matrix(nrow=length(sp.i), ncol=length(sp.i))
colnames(pair.d.i) <- sp.i
rownames(pair.d.i) <- sp.i
count <- 0
par(mfrow=c(6,4))
for(i in 1:length(sp.i)){
    if(i == length(sp.i))
        break
    for(j in (i+1):length(sp.i)){
        tmp.sp <- c(sp.i[j], sp.i[i])
        pair.d.i[i,j] <- count + 1
        pair.d.i[j,i] <- count + 2
        for(m in 1:2){
            tmp.sp <- rev(tmp.sp)
            sp.1 <- ifelse(tmp.sp[1] <= 172, ex.align.2.k2.nj$tip.label[tmp.sp[1]], as.character(tmp.sp[1]) )
            sp.2 <- ifelse(tmp.sp[2] <= 172, ex.align.2.k2.nj$tip.label[tmp.sp[2]], as.character(tmp.sp[2]) )
            count <- count + 1
            pair.d[[count]] <- column.mean.dist(int.s.inf.2[[tmp.sp[1]]]$state[,1]/10,
                                                int.s.inf.2[[tmp.sp[2]]]$state[,1]/10,
                                                10, h.breaks=20 )
##            pair.d[[count]] <- with(orth, column.mean.dist( log2(l[,sp.1]), log2(l[,sp.2]), 10, h.breaks=20 ))
            pair.d[[count]]$sp.i <- tmp.sp
            pair.d[[count]]$sp <- c(sp.1, sp.2)
##            pair.d[[count]]$sp <- sp.names[ex.align.2.k2.nj$tip.label[tmp.sp]]
        }
    }
}

par(mfrow=c(4,6))
count <- 0
for(tmp in pair.d){
    with(tmp, {
        r <- 1:(length(mids) - 1)
##        r <- 1:length(mids)
        plot(mids[r], d[r], type='b', xlab=paste(sp[1], "size"), main=paste(sp[1], sp[2], sep=", "), ylab='')
        plot.window(xlim=range(mids[r]), ylim=range(v[r]))
        lines(mids[r], v[r], type='b', col='red')
        ## plot.window(xlim=range(mids[r]), ylim=range(ks.d[r]))
        ## lines(mids[r], ks.d[r], col='blue')
        ## plot.window(xlim=range(mids[r]), ylim=range(kl.d[r]))
        ## lines(mids[r], kl.d[r], col='green')
        ## ## plot.window(xlim=range(mids[r]), ylim=range(H[r]))
        ## ## lines(mids[r], H[r], col='purple', type='b')
        ## y <- (e.lh - lh) / n
        ## plot.window(xlim=range(mids[r]), ylim=range(y[r]))
        ## lines(mids[r], y[r], col='purple', type='b')
        ## axis(4, col.axis='purple')
        ## y <- cor
        ## plot.window(xlim=range(mids[r]), ylim=range(y[r]))
        ## lines(mids[r], y[r], col='red', type='l', lwd=2)
        axis(4, line=0, col.axis='red')
        plot.window(xlim=range(mids[r]), ylim=range(delta[r]/n[r]))
        lines(mids[r], delta[r]/n[r], col='blue')
        axis(4, line=2, col.axis='blue')
        x <- all.h$mids
        y <- sapply(h, function(x){ x$density })
        plot(x, y[,1], type='n', xlab=paste(sp[2], 'size'), main=sp[2], ylab='', ylim=range(y))
        cols <- hsvScale( 1:ncol(y) ) ## , alpha=0.75)
        with(par(), rect(breaks[-length(breaks)], usr[4]-diff(usr)[3]*0.05, breaks[-1], usr[4], col=cols, border=NA))
        for(i in 1:ncol(y))
            lines(x, y[,i], col=cols[i], lwd=2)
        lines(x, all.h$density, col='grey', lwd=2)
        with(par(), text(usr[1]+diff(usr)[1]*0.1, usr[4]-diff(usr)[3]*0.1, H.all, adj=c(0,1)))
    })
    count <- count + 2
    if(count %% 24 == 0)
        inpt <- readline("next: ")
}
## these numbers are buggered up the by the long tails..     

sp.i <- c(49, 2, 3, 4, 10, 26, 78, 77, 169, 61, 127, 155)
pair.d2 <- vector(mode='list', length=sum( 2:length(sp.i) - 1 ) * 2)
pair.d2.i <- matrix(nrow=length(sp.i), ncol=length(sp.i))
count <- 0
par(mfrow=c(6,4))
for(i in 1:length(sp.i)){
    if(i == length(sp.i))
        break
    for(j in (i+1):length(sp.i)){
        tmp.sp <- c(sp.i[j], sp.i[i])
        pair.d2.i[i, j] <- count + 1
        pair.d2.i[j, i] <- count + 2
        for(m in 1:2){
            tmp.sp <- rev(tmp.sp)
            sp.1 <- ex.align.2.k2.nj$tip.label[tmp.sp[1]]
            sp.2 <- ex.align.2.k2.nj$tip.label[tmp.sp[2]]
            count <- count + 1
            pair.d2[[count]] <- with(orth, column.mean.dist.2( log2(l[,sp.1]), log2(l[,sp.2]), 10 ))
            pair.d2[[count]]$sp.i <- tmp.sp
            pair.d2[[count]]$sp <- sp.names[ex.align.2.k2.nj$tip.label[tmp.sp]]
        }
    }
}

par(mfrow=c(4,6))
count <- 1
for(tmp in pair.d2){
    with(tmp, {
        plot(x.m, d, type='b', xlab=sp[1], ylab='', main=sp[2])
        ## plot.window(xlim=range(x.m), ylim=range(v))
        ## lines(x.m, v, type='b', col='red')
        plot.window(xlim=range(x.m), ylim=range(xy.c))
        lines(x.m, xy.c, type='b', col='blue')
        axis(4)
        plot.window(xlim=range(x.m), ylim=range(kl.d))
        lines(x.m, kl.d, type='b', col='red')
        axis(4, line=2)
    })
    if(count %% 24 == 0)
        inpt <- readline("next: ")
    count <- count + 1
}

## all with takifugu
par(mfrow=c(2,2))

tak.i <- pair.d2.i[2,]
tak.i <- tak.i[ !is.na(tak.i) ]
tak.stats <- lapply(tak.i, function(x){
    with( pair.d2[[x]], cbind(n, d, v, x.m, y.m, xy.c, xv, yv, kl.d, chi.d, ks.d))
})
names(tak.stats) <- sapply(pair.d2[tak.i], function(x){ x$sp[2] })

tak.x <- sapply(tak.stats, function(x){ x[,'x.m'] })
tak.d <- sapply(tak.stats, function(x){ x[,'d'] })
tak.v <- sapply(tak.stats, function(x){ x[,'v'] })
tak.kl.d <- sapply(tak.stats, function(x){ x[,'kl.d'] })
tak.chi.d <- sapply(tak.stats, function(x){ x[,'chi.d'] })
tak.ks.d <- sapply(tak.stats, function(x){ x[,'ks.d'] })

plot(tak.x[,1], tak.d[,1], xlim=range(tak.x), ylim=range(log2(tak.d)), type='n')
for(i in 1:length(tak.stats)){
    lines(tak.x[,i], log2(tak.d[,i]), col=i, lwd=2)
}
legend( 'topright', legend=names(tak.stats), lty=1, col=1:length(tak.stats))

plot(tak.x[,1], tak.v[,1], xlim=range(tak.x), ylim=range(sqrt(tak.v)), type='n')
for(i in 1:length(tak.stats)){
    lines(tak.x[,i], sqrt(tak.v[,i]), col=i, lwd=2)
}
legend( 'topright', legend=names(tak.stats), lty=1, col=1:length(tak.stats))

vis.min.pts <- vis.min.summary()

## common minimisation:
tmp <- common.minimisation.by.ancestral.size(anc=266, d1=10, d2=3, disc.n=1, max.s2=110)
with(tmp$df, plot(as, p1, type='l'))
with(tmp$df, lines(as, p2, type='l', col='red'))
with(tmp$df, lines(as, p1.2, type='l', col='blue'))
with(tmp$df, lines(as, p1 * p2, type='l', col='blue', lty=2))
with(tmp$df, plot(as, n, type='l'))
with(tmp$df, plot(as, -log10(p), type='l'))

with(tmp$df, plot(as, log2(p1.2 / (p1 * p2)), type='l', lwd=2))
with(tmp$df, plot.window( xlim=range(as), ylim=range(-log10(p)) ))
with(tmp$df, lines(as, -log10(p), col='blue'))
axis(4)

## lets look at a simple test for common minimisation..
leaf.276 <- collect.descendants(276, leaf.only=TRUE)

com.d.min.10.276 <- lapply( leaf.276$i, function(d2){
    common.minimisation.by.ancestral.size(anc=266, d1=10, d2=d2, disc.n=1, max.s2=110)
})

b <- rep(TRUE, length(leaf.276$i)) ## int.s.b[ leaf.276$names ] ## & !(leaf.276$i %in% c(21, 6) )
pp <- -log10(sapply(com.d.min.10.276[b], function(x){ x$comb['p'] }))
pp2 <- -log10(sapply(com.d.min.10.276[b], function(x){ x$comb['p2'] }))
oe <- sapply( com.d.min.10.276[b], function(x){ x$comb['obs'] / x$comb['exp'] })
gs <- log2( genome.size[ leaf.276$names[b] ])
dd <- ex.align.2.k2[10, leaf.276$i[b]]

plot( gs, pp )
plot( gs, pp2 )
identify(gs, pp2, leaf.276$i[b])

with(com.d.min.10.276[[16]]$df, plot(as, -log10(p), type='l'))

plot(gs, oe )

plot(dd, pp2)
plot(dd, oe)


plot.commonality <- function(com, vl=c(80,100), oe.range=NULL){
    with(com, {
        plot(as, n, yaxt='n', col='grey', type='l' )
        plot.window(xlim=range(as), ylim=c(0,1))
        lines(as, p1, col='red')
        lines(as, p2, col='blue')
        lines(as, p1.2, col='purple')
        if(is.null(oe.range))
            oe.range <- range( log2( p1.2 / (p1 * p2)), ra.rm=TRUE, finite=TRUE )
        plot.window(xlim=range(as), ylim=oe.range)
        axis(2)
        abline(h=0, lty=2)
        lines(as, log2( p1.2 / (p1 * p2) ), col=rgb(0, 0.8, 0))
        max.oe <- log2( pmax(p1, p2) / (p1 * p2) )
        lines(as, max.oe, col='brown')
        pv <- -log10(com$p)
        plot.window(xlim=range(as), ylim=c(range(pv)))
        axis(4)
        lines(as, pv, col=rgb(0, 0.5, 0.5), lwd=3)
        abline(h=0, v=vl, lty=3)
        
    })
}


## sp1 should be a single species; species two a vector of species
plot.correlations.distances <- function(sp1, sp2, int.l=log2(orth$l), sp.d=ex.align.2.k2, min.l=log2(76),
                                        gen.size=genome.size, sp.b=int.s.b[ sp2 ],
                                        do.id=FALSE){
    ## enforce expectation. Should really warn, but.. 
    sp1 <- sp1[1]
    ## remove any not fitting criteria
    sp2 <- sp2[ sp.b ]
    ## use sapply in order to specify common introns for individual
    ## pairs. Return also, the number of introns used for each pair
    b1 <- int.l[,sp1] >= min.l & is.finite(int.l[,sp1])
    pair.cor <- sapply( sp2, function(sp){
        b <- b1 & (int.l[,sp] >= min.l) & is.finite(int.l[,sp])
        c( 'n'=sum(b), 'c'=cor( int.l[b, sp1], int.l[b, sp]))
    })
    gs <- log2(gen.size[sp2])
    cex <- 1 + gs - min(gs)
    plot( sp.d[ sp1, sp2 ], pair.cor['c',], cex=cex,
         pch=19)
    mod1 <- lm( pair.cor['c',] ~ sp.d[sp1, sp2] )
    mod2 <- lm( pair.cor['c',] ~ sp.d[sp1, sp2] + gs )
    mod3 <- lm( pair.cor['c',] ~ gs )
    lab <- NULL
    if(do.id){
        lab <- identify( sp.d[ sp1, sp2 ], pair.cor['c',], sp2, pos=TRUE )
    }
    invisible(list(lab=lab, mod1=mod1, mod2=mod2, mod3=mod3,
                   d=sp.d[ sp1, sp2 ], cor=pair.cor['c',], gs=gs,
                   sp1=sp1, sp2=sp2))
}

## test the above function and look at the models:
sp1 <- 'gasterosteus_aculeatus'
sp2 <- rownames(sp.class)[ sp.class[,'teleostei'] ]
sp2 <- sp2[-4] ## remove gasterosteus itself

par(mfrow=c(2,2))
tmp <- plot.correlations.distances(sp1, sp2, do.id=TRUE)
with(tmp, plot( mod1$residuals, mod2$residuals, cex=1 + gs - min(gs), pch=19,
               col=1 + as.numeric(int.s.b[sp2])))
abline(0,1, lty=2)
## small genomes benefit from the additional model component (genome size)
## whereas large genomes have higher residuals in the two component model
## (mod1 = cor ~ distance
##  mod2 = cor ~ distance + genome_size

with(tmp, plot( gs, mod1$residuals, pch=19, cex=1 + gs - min(gs),
               ylim=range(c(mod1$residuals, mod2$residuals)),
               col=1 + as.numeric(int.s.b[sp2])))
abline(h=0, lty=2)
## there is a clear relationship between the residual and genome size
with(tmp, plot( gs, mod2$residuals, pch=19,
                 cex=1 + gs - min(gs), col=1 + as.numeric(int.s.b[sp2]),
               ylim=range(c(mod1$residuals, mod2$residuals)) ))
abline(h=0, lty=2)

with(tmp, sum( mod1$residuals^2 ))
## 0.04126273
with(tmp, sum( mod2$residuals^2 ))
## 0.0251926
with(tmp, sum( mod3$residuals^2 ))
## 0.4556272
## but note that mod3 still has a p-value of: 0.001317

with(tmp, summary( lm( mod1$residuals ~ gs )))
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  1.126919   0.215864   5.220 3.62e-06 ***
## gs          -0.038190   0.007315  -5.221 3.61e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.02326 on 49 degrees of freedom
## Multiple R-squared:  0.3575,	Adjusted R-squared:  0.3443 
## F-statistic: 27.26 on 1 and 49 DF,  p-value: 3.609e-06

## this should take into account the additional model
with(tmp, summary( lm( mod2$residuals ~ gs )))
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)
## (Intercept) -2.464e-17  2.104e-01       0        1
## gs           8.278e-19  7.130e-03       0        1

## Residual standard error: 0.02267 on 49 degrees of freedom
## Multiple R-squared:  4.684e-34,	Adjusted R-squared:  -0.02041 
## F-statistic: 2.295e-32 on 1 and 49 DF,  p-value: 1

## and lets try the same thing with B splendens
sp1 <- 'betta_splendens'
sp2 <- rownames(sp.class)[ sp.class[,'teleostei'] ]
sp2 <- sp2[-3] ## remove B splendens

par(mfrow=c(2,2))
tmp <- plot.correlations.distances(sp1, sp2)
with(tmp, points(d, mod1$fitted.values, cex=1 + gs - min(gs), col='red', lwd=3 ))
with(tmp, points(d, mod2$fitted.values, cex=1 + gs - min(gs), col='blue', lwd=3 ))
with(tmp, plot( mod1$residuals, mod2$residuals, cex=1 + gs - min(gs), pch=19,
               col=1 + as.numeric(int.s.b[sp2])))
abline(0,1, lty=2)
## mod1 = cor ~ distance
## mod2 = cor ~ distance + genome_size

with(tmp, plot( gs, mod1$residuals, pch=19, cex=1 + gs - min(gs),
               ylim=range(c(mod1$residuals, mod2$residuals)),
               col=1 + as.numeric(int.s.b[sp2])))
abline(h=0, lty=2)
## there is a clear relationship between the residual and genome size
with(tmp, plot( gs, mod2$residuals, pch=19,
                 cex=1 + gs - min(gs), col=1 + as.numeric(int.s.b[sp2]),
               ylim=range(c(mod1$residuals, mod2$residuals)) ))
abline(h=0, lty=2)

with(tmp, sum( mod1$residuals^2 ))
## 0.05008009
with(tmp, sum( mod2$residuals^2 ))
## 0.0297754
with(tmp, sum( mod3$residuals^2 ))
## 0.3197196

summary(tmp$mod1)  ## p <= 2.2e-16
summary(tmp$mod2)  ## sp.d <= 2e-16, gs 6.66e-07
summary(tmp$mod3)  ## p <= 0.0002725 (negatively associated)

with(tmp, summary( lm( mod1$residuals ~ gs )))
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  1.268639   0.236709   5.359 2.23e-06 ***
## gs          -0.042991   0.008021  -5.360 2.23e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.02538 on 49 degrees of freedom
## Multiple R-squared:  0.3696,	Adjusted R-squared:  0.3568 
## F-statistic: 28.73 on 1 and 49 DF,  p-value: 2.226e-06

## this though has no significance. Which means that the
## model fits both bit and small genomes equally badly
with(tmp, summary( lm( mod1$residuals^2 ~ gs )))

## this gives nothing significance as it should.
with(tmp, summary( lm( mod2$residuals ~ gs )))

## and lets try the same thing with D clupeoides
sp1 <- 'denticeps_clupeoides'
sp2 <- rownames(sp.class)[ sp.class[,'teleostei'] ]
sp2 <- sp2[-10] ## remove D clupeoides

par(mfrow=c(2,2))
tmp <- plot.correlations.distances(sp1, sp2, do.id=TRUE, sp.b=TRUE)
with(tmp, points(d, mod1$fitted.values, cex=1 + gs - min(gs), col='red', lwd=3 ))
with(tmp, points(d, mod2$fitted.values, cex=1 + gs - min(gs), col='blue', lwd=3 ))
with(tmp, plot( mod1$residuals, mod2$residuals, cex=1 + gs - min(gs), pch=19,
               col=1 + as.numeric(int.s.b[sp2])))
abline(0,1, lty=2)

## we can also try to show these correlations along a tree..
par(mfrow=c(1,2))
td.266 <- draw.tree( ex.align.2.k2.nj, tree.lines, root=266, label.nodes=FALSE, rm=0.2, leaf.cex=0.75, leaf.font=3 )
plot.new()
plot.window( xlim=range(tmp$cor^2), ylim=td.266$usr[3:4], yaxs='i' )
x0 <- with(par(), usr[1])
with(td.266, segments(x0, lines$y[ nodes$leaf.b ], tmp$cor[ nodes$sp ]^2, col=hsvScale(tmp$gs[nodes$sp]), lwd=4) )

## (mod1 = cor ~ distance
##  mod2 = cor ~ distance + genome_size
## The summary of mod2 suggests that the genome size is a slighly worse
## predictor of intron size correlation than phylogenetic distance

summary(tmp$mod1) ## p <= 0.0007833
summary(tmp$mod2) ## p <= 1e-07, sp.d 6.72e-7, gs 4.93e-6
summary(tmp$mod3) ## p <= 0.006618

with(tmp, plot( gs, mod1$residuals, pch=19, cex=1 + gs - min(gs),
               ylim=range(c(mod1$residuals, mod2$residuals)),
               col=1 + as.numeric(int.s.b[sp2])))
abline(h=0, lty=2)

with(tmp, plot( gs, mod1$residuals, pch=19,
                 cex=1 + gs - min(gs), col=1 + as.numeric(int.s.b[sp2]),
               ylim=range(c(mod1$residuals, mod2$residuals)) ))
abline(h=0, lty=2)
## there is a clear relationship between the residual and genome size
## with model 1,
with(tmp, plot( gs, mod2$residuals, pch=19,
                 cex=1 + gs - min(gs), col=1 + as.numeric(int.s.b[sp2]),
               ylim=range(c(mod1$residuals, mod2$residuals)) ))
abline(h=0, lty=2)
## but not with model 2?

with(tmp, sum( mod1$residuals^2 ))
## 0.05113176
with(tmp, sum( mod2$residuals^2 ))
## 0.03296174
with(tmp, sum( mod3$residuals^2 ))
## 0.05541899
## Genome size is almost as good as distance in predicting
## the correlation.

with(tmp, summary( lm( mod1$residuals ~ gs )))
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  1.176883   0.242113   4.861 1.24e-05 ***
## gs          -0.039892   0.008206  -4.861 1.24e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.02653 on 49 degrees of freedom
## Multiple R-squared:  0.3254,	Adjusted R-squared:  0.3116 
## F-statistic: 23.63 on 1 and 49 DF,  p-value: 1.242e-05

## p-value is 1 from this... 
with(tmp, summary( lm( mod2$residuals ~ gs )))

## and lets try the same thing with D rerio
sp1 <- 'danio_rerio'
sp2 <- rownames(sp.class)[ sp.class[,'teleostei'] ]
sp2 <- sp2[sp2 != 'danio_rerio'] ## remove D rerio

par(mfrow=c(2,2))
tmp <- plot.correlations.distances(sp1, sp2, do.id=TRUE)
with(tmp, points(d, mod1$fitted.values, cex=1 + gs - min(gs), col='red', lwd=3 ))
with(tmp, points(d, mod2$fitted.values, cex=1 + gs - min(gs), col='blue', lwd=3 ))
with(tmp, plot( mod1$residuals, mod2$residuals, cex=1 + gs - min(gs), pch=19,
               col=1 + as.numeric(int.s.b[sp2])))
abline(0,1, lty=2)

## (mod1 = cor ~ distance
##  mod2 = cor ~ distance + genome_size
## The summary of mod2 suggests that the genome size is a slighly worse
## predictor of intron size correlation than phylogenetic distance
summary(tmp$mod1)  ## p <= 2.0e-5
summary(tmp$mod2)  ## p <= 2.4e-7, sp.d 1.3e-7, gs 0.000427
summary(tmp$mod3)  ## p <= 0.1249 hopeless on its own.

with(tmp, plot( gs, mod1$residuals, pch=19, cex=1 + gs - min(gs),
               ylim=range(c(mod1$residuals, mod2$residuals)),
               col=1 + as.numeric(int.s.b[sp2])))
abline(h=0, lty=2)
## there is a somewhat relationship between the residual and genome size
with(tmp, plot( gs, mod2$residuals, pch=19,
                 cex=1 + gs - min(gs), col=1 + as.numeric(int.s.b[sp2]),
               ylim=range(c(mod1$residuals, mod2$residuals)) ))
abline(h=0, lty=2)  ## does not look that different from the above.. 

with(tmp, sum( mod1$residuals^2 ))
## 0.03413675
with(tmp, sum( mod2$residuals^2 ))
## 0.02629006
with(tmp, sum( mod3$residuals^2 ))
## 0.04723896
## but these are not that different from each other.. ?

with(tmp, summary( lm( mod1$residuals ~ gs )))
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.794845   0.220060   3.612 0.000714 ***
## gs          -0.026965   0.007465  -3.612 0.000713 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.02346 on 49 degrees of freedom
## Multiple R-squared:  0.2103,	Adjusted R-squared:  0.1942 
## F-statistic: 13.05 on 1 and 49 DF,  p-value: 0.0007134
## p-value is 1 from this... 

with(tmp, summary( lm( mod2$residuals ~ gs )))
## p = 1

## and lets try the same thing with E. electricus
sp1 <- 'electrophorus_electricus'
sp2 <- rownames(sp.class)[ sp.class[,'teleostei'] ]
sp2 <- sp2[sp2 != 'electrophorus_electricus'] ## remove D rerio

par(mfrow=c(2,2))
tmp <- plot.correlations.distances(sp1, sp2, do.id=TRUE)
## with(tmp, points(d, mod1$fitted.values, cex=1 + gs - min(gs), col='red', lwd=3 ))
## with(tmp, points(d, mod2$fitted.values, cex=1 + gs - min(gs), col='blue', lwd=3 ))
with(tmp, plot( mod1$residuals, mod2$residuals, cex=1 + gs - min(gs), pch=19,
               col=1 + as.numeric(int.s.b[sp2])))
abline(0,1, lty=2)

with(tmp, plot(gs, cor, pch=19, cex=2))

## (mod1 = cor ~ distance
##  mod2 = cor ~ distance + genome_size
## The summary of mod2 suggests that the genome size is a slighly worse
## predictor of intron size correlation than phylogenetic distance
summary(tmp$mod1) ## p <= 1.2e-9
summary(tmp$mod2) ## p <= 6.3e-14, sp.d 1e-14, gs 9.7e-7
summary(tmp$mod3) ## 0.476

with(tmp, plot( gs, mod1$residuals, pch=19, cex=1 + gs - min(gs),
               ylim=range(c(mod1$residuals, mod2$residuals)),
               col=1 + as.numeric(int.s.b[sp2])))
abline(h=0, lty=2)
## there is a clear relationship between the residual and genome size
with(tmp, plot( gs, mod2$residuals, pch=19,
                 cex=1 + gs - min(gs), col=1 + as.numeric(int.s.b[sp2]),
               ylim=range(c(mod1$residuals, mod2$residuals)) ))
abline(h=0, lty=2)


with(tmp, sum( mod1$residuals^2 ))
## 0.05647369
with(tmp, sum( mod2$residuals^2 ))
## 0.03409357
with(tmp, sum( mod3$residuals^2 ))
## 0.1196849

with(tmp, summary( lm( mod1$residuals ~ gs )))
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.794845   0.220060   3.612 0.000714 ***
## gs          -0.026965   0.007465  -3.612 0.000713 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.02346 on 49 degrees of freedom
## Multiple R-squared:  0.2103,	Adjusted R-squared:  0.1942 
## F-statistic: 13.05 on 1 and 49 DF,  p-value: 0.0007134
## p-value is 1 from this... 

with(tmp, summary( lm( mod2$residuals ~ gs )))

com.min.10.276 <- t(sapply( leaf.276$i, function(d2){
    common.minimisation(266, 10, d2, anc.range=c(10*log2(256), 10*log2(1024)))}))

## 21 and 6 are P. magnuspinnatus and H. comes
## they disturb some of the observations below..
## It's not really OK to remove them, but for the com.long below
## it makes a clearer picture.
b <- int.s.b[ leaf.276$names ] ## & !(leaf.276$i %in% c(21, 6) )
plot( log2(genome.size[ leaf.276$names[b] ]), com.min.10.276[b,'oe'] )
summary( lm( com.min.10.276[b,'oe'] ~ log2(genome.size[ leaf.276$names[b] ])) )

plot( log2(genome.size[ leaf.276$names[b] ]), -log10(com.min.10.276[b,'p']) )
summary( lm( -log10(com.min.10.276[b,'p']) ~ log2(genome.size[ leaf.276$names[b] ])) )

## No relationship between genome size and common minimisation.
## We can compare to doing the same thing for common maintenance
com.long.10.276 <- t(sapply( leaf.276$i, function(d2){
    common.long(266, 10, d2, min.s2=10*log2(256), anc.range=c(10*log2(1024), 10*log2(2048)))}))

b <- int.s.b[ leaf.276$names ] ## & !(leaf.276$i %in% c(21, 6))
plot( log2(genome.size[ leaf.276$names[b] ]), com.long.10.276[b,'oe'] )
identify( log2(genome.size[ leaf.276$names[b] ]), com.long.10.276[b,'oe'], leaf.276$names[b] )
plot( ex.align.2.k2[10, leaf.276$i[b]], com.long.10.276[b,'oe'] )

summary( lm( com.long.10.276[b,'oe'] ~ log2(genome.size[ leaf.276$names[b] ])) ) ## p <= 0.0002213
summary( lm( com.long.10.276[b,'oe'] ~ ex.align.2.k2[10, leaf.276$i[b]] )) ## p <= 0.567
summary( lm( com.long.10.276[b,'oe'] ~ ex.align.2.k2[10, leaf.276$i[b]] + log2(genome.size[ leaf.276$names[b] ]) )) ## p <= 0.0008648

plot( log2(genome.size[ leaf.276$names[b] ]), -log10(com.long.10.276[b,'p']) )
summary( lm( -log10(com.long.10.276[b,'p']) ~ log2(genome.size[ leaf.276$names[b] ])) ) ## p <= 0.012

## and use a discretised function that recombines the p-values into one.
com.d.long.10.276.256 <- lapply( leaf.276$i, function(d2){
    common.long.by.ancestral.size( 266, 10, d2, min.s2=10*log2(256), disc.n=1)
})

## what if we use different sizes.. ?
com.d.long.10.276.128 <- lapply( leaf.276$i, function(d2){
    common.long.by.ancestral.size( 266, 10, d2, min.s2=10*log2(128), disc.n=1)
})

com.d.long.10.276.512 <- lapply( leaf.276$i, function(d2){
    common.long.by.ancestral.size( 266, 10, d2, min.s2=10*log2(512), disc.n=1)
})

com.d.long.10.276.1024 <- lapply( leaf.276$i, function(d2){
    common.long.by.ancestral.size( 266, 10, d2, min.s2=10*log2(1024), disc.n=1)
})

## And some plotting of summary statistics.. 
b <- rep(TRUE, length(leaf.276$i))  ## int.s.b[ leaf.276$names ] & !(leaf.276$i %in% c(21, 6))
pp <- -log10(sapply(com.d.long.10.276.256[b], function(x){ x$comb['p'] }))
names(pp) <- leaf.276$names[b]
pp2 <- -log10(sapply(com.d.long.10.276.256[b], function(x){ x$comb['p2'] }))
names(pp2) <- leaf.276$names[b]
oe <- sapply( com.d.long.10.276.256[b], function(x){ x$comb['obs'] / x$comb['exp'] })
names(oe) <-  leaf.276$names[b]
gs <- log2(genome.size[ leaf.276$names[b] ])
dd <- ex.align.2.k2[10, leaf.276$i[b]]

par(mfrow=c(1,2))
par(mar=c(4.1, 2.1, 2.1, 0.1))
td.276 <- draw.tree( ex.align.2.k2.nj, tree.lines, root=276, label.nodes=FALSE, rm=0.2, leaf.cex=1, leaf.font=3 )
par(mar=c(4.1, 0.1, 2.1, 2.1))
plot.new()
plot.window( xlim=c(0, max(pp2)), ylim=td.276$usr[3:4], yaxs='i' )
with(td.276, segments( 0, lines$y[ nodes$leaf.b ], pp2[ nodes$sp ], col=hsvScale(gs[nodes$sp]), lwd=4 ))
axis(1)


## let us have a look the NCBI taxonomies
ncbi.lin <- readRDS("../R_172_genomes/ncbi_lineages.rds")
names(ncbi.lin) <- sub(" ", "_", names(ncbi.lin))
## one entry for each leaf. Ordered by genome size, but with appropriate
## names.

## to get the taxons represented in a given tree we can do something like:
leaf.276.taxons <- lapply( ncbi.lin[ leaf.276$names ], function(x){ x[,c('id', 'x', 'name')] })

length(leaf.276.taxons)
tmp <- sort( table( do.call( rbind, leaf.276.taxons )[,'id'] ) )
## 25 ids are all 42..

leaf.276.taxons <- lapply( leaf.276.taxons, function(x){
    cbind( tmp[ as.character(x[,'id']) ], x, stringsAsFactors=FALSE)
})

all(sapply( leaf.276.taxons, function(x){ "Percomorphaceae" %in% x$name }))
## TRUE:
## All (42) descendants of 276 are all members of Percomorphaceae
## no species are Percomorphaceae and do not descend from 276



KL.266 <- clade.KL(266)
KL.276 <- clade.KL(276)  ## all Percomorphaceae (see above line 2132)
KL.247 <- clade.KL(247) ## amphibians, reptiles, birds
KL.253 <- clade.KL(253) ## almost everything that is not teleost, jawless or gar
KL.259 <- clade.KL(259) ## all jawed vertebrates


im.cols <- hcl.colors(256, "YlOrRd", rev=TRUE)

cairo_pdf("KL_253_root.pdf", width=8, height=12, onefile=TRUE)
par(mfrow=c(3,2))
for(sp in names(KL.253$root)){
    with(KL.253$root[[sp]], {
        image(b1, b2, f2, main=sp, col=im.cols)
        image(b1, b2, LR, main=sp, col=im.cols)
        plot(b1[-1] - diff(b1)/2, rowSums(KL), type='l', main=sp)
        plot(b2[-1] - diff(b2)/2, colSums(KL), type='l', main=sp)
        plot(b1[-1] - diff(b1)/2, rowSums(chi.d, na.rm=TRUE), type='l', main=sp)
        plot(b2[-1] - diff(b2)/2, colSums(chi.d, na.rm=TRUE), type='l', main=sp)
    })
##    input <- readline("next: ")
}
dev.off()

## Example plots
par(mfrow=c(3,2))
##sp1 <- 'danio_rerio'
sp1 <- 'denticeps_clupeoides'
for(sp2 in names(KL.266$leaves[[sp1]])){
    with(KL.266$leaves[[sp1]][[sp2]], {
        image(b1, b2, f2, main=sp2, col=im.cols)
        image(b1, b2, LR, main=sp2, col=im.cols)
        x1 <- b1[-1] - diff(b1)
        x2 <- b2[-1] - diff(b2)
        plot(x1, rowSums(KL), type='l', main=sp2)
        plot.window(xlim=range(x1), ylim=range(xf))
        lines(x1, xf, col='red')
        plot(x2, colSums(KL), type='l', main=sp2)
        plot.window(xlim=range(x2), ylim=range(yf))
        lines(x2, yf, col='red')
        plot(x1, rowSums(chi.d), type='l')
        plot(x2, colSums(chi.d), type='l')
    })
    input <- readline("next: ")
}


## take mutual information information and plot in some
## reasonable manner to summarise.
## plot extant pairwise mutual information..
plot.x.mi <- function(mi, sp1='denticeps_clupeoides', exclude.sp=NULL, by.sp1=TRUE,
                      sp1.h=by.sp1, min.s=76, breaks.n=20, plot.tp=c('KL', 'chi.d', 'LR')[1],
                      col.by=c("gs", "k2", "sp")[1], xlab=NA, ylab=NA, leg.cex=0.5, leg.pos='topright',
                      ...){
    b2mids <- function(x){ x[-1] - diff(x)/2 }
    sp2 <- names(mi$leaves[[sp1]])
    sp2 <- sp2[ !(sp2 %in% c(sp1, exclude.sp)) ]
    if(by.sp1){
        b.id <- 'b1'
        sum.f <- rowSums
        if(is.na(xlab))
            xlab <- paste(sp1, "size")
    }else{
        b.id <- 'b2'
        sum.f <- colSums
        if(is.na(xlab))
            xlab <- "size"
    }
    x <- sapply(mi$leaves[[sp1]][sp2], function(z){
        b2mids( z[[b.id]] ) })
    if(plot.tp == 'KL'){
        x.mi <- sapply(mi$leaves[[sp1]][sp2], function(z){
            sum.f(z$KL, na.rm=TRUE) })
    }
    if(plot.tp == 'chi.d'){
        x.mi <- sapply(mi$leaves[[sp1]][sp2], function(z){
        sum.f(z$chi.d, na.rm=TRUE) })
    }
    if(plot.tp == 'LR'){
        x.mi <- sapply(mi$leaves[[sp1]][sp2], function(z){
        sum.f(z$LR, na.rm=TRUE) })
    }
    plot.new()
    h <- NULL
    if(by.sp1 && sp1.h){
        s <- log2(orth$l[,sp1])
        s <- s[ !is.na(s) & s >= log2(min.s) ]
        h <- hist(s, breaks=breaks.n, plot=FALSE)
        plot.window(xlim=range(x), ylim=range(h$counts))
        rect( h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts, col=rgb(0.8, 0.8, 0.8), border=NA )
    }
    par('new'=TRUE)
    plot(1,1, type='n', xlim=range(x), ylim=range(x.mi), xlab=xlab, ylab=ylab, frame.plot=FALSE, ...)
    par('new'=FALSE)
    cols.v <- log2(genome.size[ colnames(x) ])
    cols <- hsvScale(cols.v)
    if(col.by[1] == 'k2'){
        cols.v <- ex.align.2.k2[sp1, sp2]
        cols <- hsvScale( cols.v )
    }
    if(col.by[1] == 'sp')
        cols <- sp.col[sp2]
    for(i in 1:ncol(x))
        lines(x[,i], x.mi[,i], col=cols[i])
    ## make a legend. If by species then we use class.col
    ## otherwise the cols.v and make a heatmap
    if(col.by[1] == 'sp'){
        legend(legpos, legend=names(class.col), lty=1, col=class.col, cex=leg.cex)
    }else{
        plot.window(xlim=c(0,1), ylim=c(0,1))
        v <- seq(min(cols.v), max(cols.v), length.out=100)
        vc <- hsvScale(v)
        leg.y <- seq(0.5, 0.9, length.out=length(v))
        leg.yd <- diff(leg.y)[1]
        rect(0.95, leg.y, 1.0, leg.y+leg.yd, col=vc, border=NA)
        digits <- floor( 2 - log10(max(abs(v))) )
        digits <- ifelse(digits < 0, 0, digits)
        tv <- unique( round( seq(min(v), max(v), length.out=5), digits=digits ))
        ty <- leg.y[1] + leg.yd/2 + diff(range(leg.y)) * (tv - min(v)) / diff(range(v))
        segments(0.925, ty, 0.95, ty)
        text( 0.92, ty, tv, adj=c(1,0.5), cex=leg.cex )
    }
    invisible(list(sp1=sp1, sp2=sp2, x=x, mi=x.mi, h=h, col=cols))
}

leaves.267 <- collect.leaves(267)

## some examples
plot.x.mi( KL.247, sp1='gallus_gallus',  exclude.sp=c(names(int.s.b[ !int.s.b ])), breaks.n=100 )
plot.x.mi( KL.247, sp1='gopherus_agassizii',  exclude.sp=c(names(int.s.b[ !int.s.b ])), breaks.n=100 )
plot.x.mi( KL.247, sp1='crocodylus_porosus',  exclude.sp=c(names(int.s.b[ !int.s.b ])), breaks.n=100 )
plot.x.mi( KL.247, sp1='anolis_carolinensis',  exclude.sp=c(names(int.s.b[ !int.s.b ])), breaks.n=100 )
plot.x.mi( KL.247, sp1='pogona_vitticeps',  exclude.sp=c(names(int.s.b[ !int.s.b ])), breaks.n=100 )
plot.x.mi( KL.247, sp1='sphenodon_punctatus',  exclude.sp=c(names(int.s.b[ !int.s.b ])), breaks.n=100 )


### Simple stats discretised by the level of each stat:
### calulates the correlation and mutual information of the
### total sets or by the discretised sets.
cor.mi.259.hs.2 <- stats.discretise( 'homo_sapiens', 259, brk.l=3 )
cor.mi.259.dr.2 <- stats.discretise( 'danio_rerio', 259, brk.l=3 )
cor.mi.259.dc.2 <- stats.discretise( 'denticeps_clupeoides', 259, brk.l=3 )
cor.mi.259.bs.2 <- stats.discretise( 'betta_splendens', 259, brk.l=3 )

## convert to a more convenient table.. 
hs.stats.2 <- t(sapply( cor.mi.259.hs.2, function(x){
    c(x$cor.x, x$cor.y, x$mi.x, x$mi.y, x$cor.all, x$mi.all) }))

dr.stats.2 <- t(sapply( cor.mi.259.dr.2, function(x){
    c(x$cor.x, x$cor.y, x$mi.x, x$mi.y, x$cor.all, x$mi.all) }))

dc.stats.2 <- t(sapply( cor.mi.259.dc.2, function(x){
    c(x$cor.x, x$cor.y, x$mi.x, x$mi.y, x$cor.all, x$mi.all) }))

bs.stats.2 <- t(sapply( cor.mi.259.bs.2, function(x){
    c(x$cor.x, x$cor.y, x$mi.x, x$mi.y, x$cor.all, x$mi.all) }))

## these are explored further down in this file and compard to data
## generated using a random model of intron size evolution.

make.summary.figure <- function(fname=NULL){
    adj=1
    padj=0.16
    line=1.5
    if(!is.null(fname)){
        ## for 4 rows:
        ##cairo_pdf( fname, width=full.w*pdf.m, height=full.w*pdf.m*1.25 )
        ## for 3 rows:
        cairo_pdf( fname, width=full.w*pdf.m, height=full.w*pdf.m*1 )
    }
    ## make a figure with 4 columns and 4 rows and see how that works out..
    ## we have to use split.screen since some underlying functions use that.
    close.screen(all.screens=TRUE)
##    panels <- split.screen( c(4,3) )
    panels <- split.screen( c(3,3) )
    h2.i <- c(78, 2, 127) ##,
        ## 266, 253, 248,
        ## 276, 21, 6,
        ## 282)
    h2.names <- c('D rerio', 'T rubripes', 'M musculus') ##,
        ## 'teleost root', 'super-tetrapod root', 'tetrapod root',
        ## '276', 'P magnuspinnatus', 'H comes', '282'
        ## )
    h2.xi <- c(rep(259, 3))
    h2.xn <- c(rep('jawed root', 3))
    ## first show discretised distributions for 259 and D. rerio
    ax.cex <- 0.5
    ax.mgp <- c(1,0.2,0)
    ax.tcl <- -ax.mgp[2]
    cex.xylab <- 0.75
    cex.leg <- 0.65
    panel <- panels[1]
    mar <- c(2.6, 2.6, 1.1, 2.6)
    with( pair.d[[ pair.d.i['259','78'] ]], {
        screen(panel)
        par(mar=mar)
        x <- all.h$mids
        y <- sapply(h, function(x){ x$density })
        sp.lab <- ifelse( is.na(sp.names[sp]), sp, sp.names[sp] )
        plot(x, y[,1], type='n', xlab=paste(sp.lab[2], 'size'), main="", ylab='density', ylim=range(y),
             xlim=range(c(breaks, x)),
             mgp=ax.mgp, tcl=ax.tcl, cex.lab=cex.xylab, cex.axis=ax.cex)
        cols <- hsvScale( 1:ncol(y), alpha=1)
        with(par(), rect(breaks[-length(breaks)], usr[4]-diff(usr)[3]*0.05, breaks[-1], usr[4], col=cols, border=NA))
        for(i in 1:ncol(y))
            lines(x, y[,i], col=cols[i], lwd=2)
        lines(x, all.h$density, col='grey', lwd=2)
        with(par(), mtext(LETTERS[panel], side=2, at=usr[4], line=line, cex=mt.cex, las=1, padj=padj, adj=adj))
    })
    panel <- panel + 1
    ## kernel <- dnorm(0:20, sd=6)  ## used as default by h.tform (but needs to be defined outside of this function)
    for(i in 1:length(h2.i)){
        tmp <- hist.2d.2(int.s.inf.2[[h2.xi[i]]]$state[,1], int.s.inf.2[[h2.i[i]]]$state[,1], panel,
                         dist.w=0.25, dist.h=0.3, h.tform=h.tform, mar=mar, ax.cex=0.5, ax.mgp=c(3,0.2,0),
                         xlab=h2.xn[i], main="", ylab=h2.names[i],
                         trace.modes=TRUE, x.div=10, y.div=10, panel.label=LETTERS[panel], panel.padj=padj, panel.adj=adj,
                         panel.line=1)
        panel <- panel+1
    }
    screen(panel)
    par(mar=mar)
    vis.mod.delta(alpha=0.2, lwd=2, cex.leg=cex.leg, cex.xylab=cex.xylab, cex.axis=ax.cex, mgp.axis=ax.mgp, xlab='Jawed ancestor size',
                  ylab='Extant - ancestral size')
    with(par(), mtext(LETTERS[panel], side=2, at=usr[4], line=line, cex=mt.cex, las=1, adj=adj, padj=padj))
    panel <- panel+1
    screen(panel)
    par(mar=mar)
    vis.jawed.minim(cex.axis=ax.cex, mgp.axis=ax.mgp, cex.xylab=cex.xylab, line.xylab=1, cex.leg=cex.leg,
                    mar=mar)
    with(par(), mtext(LETTERS[panel], side=2, at=usr[4], line=line, cex=mt.cex, las=1, padj=padj, adj=adj))
    panel <- panel+1
    screen(panel)
    par(mar=mar)
    vis.min.summary(id.pts=NULL, cex.lab=cex.xylab, cex.xylab=cex.xylab, cex.leg=cex.leg, cex.axis=ax.cex, mgp.axis=ax.mgp)
    ##
    with(par(), mtext(LETTERS[panel], side=2, at=usr[4], line=line, cex=mt.cex, las=1, adj=adj, padj=padj))
    panel <- panel+1
    screen(panel)
    sp1 <- 'danio_rerio'
    sp2 <- 'oryzias_latipes'
    with(orth, {
        x <- log2(l[,sp1])
        y <- log2(l[,sp2])
        b <- !(is.na(x) | is.na(y)) & x >= log2(76) & y >= log2(76)
        par(mar=mar)
        plot(x[b], y[b], cex=0.4, col=rgb(0,0,0,0.05), cex.lab=cex.xylab, cex.axis=ax.cex, mgp=ax.mgp, tcl=ax.tcl,
             xlab='D. rerio', ylab='O. latipes')
        with(par(), mtext(LETTERS[panel], side=2, at=usr[4], line=line, cex=mt.cex, las=1, adj=adj, padj=padj))
    })
    ## show the regions used to compare correlations:
    dr.l <- cor.mi.259.dr.2[['oryzias_latipes']]$x.brk
    with(par(), rect(dr.l[-length(dr.l)], usr[3], dr.l[-1], usr[4], col=rgb(0,0,0,0.05)))
    panel <- panel+1
    screen(panel)
    par(mar=mar)
    ## show the increase in correlation for long introns
    plot(ex.align.2.k2['danio_rerio', rownames(dr.stats.2)],
         dr.stats.2[,2] / dr.stats.2[,1], col=sp.col[rownames(dr.stats.2)], pch=19,
         cex.lab=cex.xylab, cex.axis=ax.cex, mgp=ax.mgp, tcl=ax.tcl, xlab='D. rerio Kimura 2-factor',
         ylab='correlation ratio (long / short)')
        with(par(), mtext(LETTERS[panel], side=2, at=usr[4], line=line, cex=mt.cex,
                          las=1, adj=adj, padj=padj))
    legend('topleft', legend=names(class.col), col=class.col, pch=19, cex=cex.leg)
    ## The commented out code creates beautiful plots arguing that long introns contribute
    ## to the mutual information in intron lengths between species. Unfortunately, it
    ## appears that it is possible that these patterns could arise from random but restricted
    ## evolution of size, and I do not feel comfortable to use them as evidence.
    ## sp1.l <- mk.species.label(sp1)
    ## sp2.l <- mk.species.label(sp2)
    ## with( KL.266$leaves[[sp1]][[sp2]], image.marginal(b1, b2, LR, dist.w=0.25, dist.h=0.3,
    ##                                                   mar=mar, ax.cex=ax.cex, ax.mgp=c(3,0.2,0),
    ##                                                   panel.label=LETTERS[panel], panel.padj=padj,
    ##                                                   panel.adj=adj, panel.line=1,
    ##                                                   xlab=sp1.l, ylab=sp2.l))
##     panel <- panel+1
##     screen(panel)
##     par(mar=mar)
##     with( KL.266$leaves[[sp1]][[sp2]], image.marginal(b1, b2, KL, dist.w=0.25, dist.h=0.3,
##                                                       mar=mar, ax.cex=ax.cex, ax.mgp=c(3,0.2,0),
##                                                       panel.label=LETTERS[panel], panel.padj=padj,
##                                                       panel.adj=adj, panel.line=1,
##                                                       xlab=sp1.l, ylab=sp2.l))
## ###
##     panel <- panel + 1
##     screen(panel)
##     par(mar=mar)
##     plot.x.mi( KL.266, sp1='danio_rerio', exclude.sp=c(names(int.s.b[!int.s.b])), breaks.n=100, col.by='k2',
##               xlab="D. rerio size", ylab="sum KLD",  cex.lab=cex.xylab, cex.axis=ax.cex, mgp=ax.mgp, tcl=ax.tcl)
##     with(par(), mtext(LETTERS[panel], side=2, at=usr[4], line=line, cex=mt.cex, las=1))
##     panel <- panel + 1;
##     screen(panel)
##     par(mar=mar)
##     plot.x.mi( KL.259, sp1='homo_sapiens', exclude.sp=c(names(int.s.b[!int.s.b])), breaks.n=100, col.by='k2',
##               xlab="H. sapiens size", ylab="sum KLD", cex.lab=cex.xylab, cex.axis=ax.cex, mgp=ax.mgp, tcl=ax.tcl)
##     with(par(), mtext(LETTERS[panel], side=2, at=usr[4], line=line, cex=mt.cex, las=1))
    if(!is.null(fname))
        dev.off()
}

## seems that I need to call this here.. 
kernel <- dnorm(0:20, sd=6)  ## used as default by h.tform         
make.summary.figure('intron_size_evolution_summary.pdf')


## example plots showing the realtionship between divergence and the difference
## in the correlation for long and short introns

plot( ex.align.2.k2['betta_splendens', rownames(bs.stats.2)], bs.stats.2[,4] / bs.stats.2[,3],
     col=sp.col[ rownames(bs.stats.2) ],pch=19)

plot( ex.align.2.k2['betta_splendens', rownames(bs.stats.2)], bs.stats.2[,6] / bs.stats.2[,5],
     col=sp.col[ rownames(bs.stats.2) ],pch=19)

plot( ex.align.2.k2['betta_splendens', rownames(bs.stats.2)], bs.stats.2[,8] / bs.stats.2[,7],
     col=sp.col[ rownames(bs.stats.2) ],pch=19)

plot( ex.align.2.k2['danio_rerio', rownames(dr.stats.2)], dr.stats.2[,2] / dr.stats.2[,1],
     col=sp.col[ rownames(dr.stats.2) ],pch=19)

plot( ex.align.2.k2['denticeps_clupeoides', rownames(dc.stats.2)], dc.stats.2[,2] / dc.stats.2[,1],
     col=sp.col[ rownames(dc.stats.2) ],pch=19)

plot( ex.align.2.k2['homo_sapiens', rownames(hs.stats.2)], hs.stats.2[,2] / hs.stats.2[,1],
     col=sp.col[ rownames(hs.stats.2) ],pch=19)

plot( ex.align.2.k2['betta_splendens', rownames(bs.stats.2)], bs.stats.2[,2] / bs.stats.2[,1],
     col=sp.col[ rownames(bs.stats.2) ],pch=19)

hist( bs.stats.2[,1] )
hist( bs.stats.2[,2] )
hist( dr.stats.2[,1] )
hist( dr.stats.2[,2] )

plot( log2(genome.size[rownames(hs.stats.2)]), hs.stats.2[,2] / hs.stats.2[,1],
     col=sp.col[ rownames(hs.stats.2) ],pch=19)

## Both the discretised correlation and the mutual information seems to be informative
## in that we see more correlation and more mutual information for the longer
## sets of intron sizes. Still, as ever, we do need to ask if this could arise
## from neutral evolution of sizes:

## The following code is a rather inefficient means to generate randomly
## evolved sets of intron sizes and comparing them to each other in the same
## way as was done to get the stats for the comparison of extant species.

## evolve two lineages from a common ancestor
evol.size <- function(anc, n, anc2=anc,
                      ev.exp.1=quote(rnorm(length(anc), -1, 2)),
                      ev.exp.2=quote(rnorm(length(anc), 1, 1)),
                      min.s=log2(76)){
    d1 <- matrix(nrow=length(anc), ncol=n+1)
    d1[,1] <- anc
    d2 <- d1
    d2[,1] <- anc2
    for(i in 2:ncol(d1)){
        delta1 <- eval(ev.exp.1)
        delta2 <- eval(ev.exp.2)
        ## mutations which give rise to introns that are too small
        ## should be strongly selected against, so we will 
        delta1[ d1[,i-1] + delta1 < min.s ] <- 0
        delta2[ d2[,i-1] + delta2 < min.s ] <- 0
        d1[,i] <- d1[,i-1] + delta1
        d2[,i] <- d2[,i-1] + delta2
    }
    list(d1=d1, d2=d2)
}

anc <- int.s.inf.2[[266]]$state[,1] / 10
anc <- anc[ anc >= log2(76) ]


desc <- lapply( 1:1000, function(i){
    evol.size(anc, 20, anc,
              ev.exp.1=quote(rnorm(length(anc), -0.1, 0.5)),
              ev.exp.2=quote(rnorm(length(anc), -0.1, 0.5)))
})

## do fewer, but with a higher rate of minimisation
## then sample across desc2 objects.. 
desc.2 <- lapply( 1:100, function(i){
    evol.size(anc, 20, anc,
              ev.exp.1=quote(rnorm(length(anc), -0.15, 0.75)),
              ev.exp.2=quote(rnorm(length(anc), -0.15, 0.75)))
})


## evolve from two real distributions instead
anc.2 <- log2( orth$l[,'homo_sapiens'] )
anc.2 <- anc.2[ anc.2 > log2(76) & !is.na(anc.2) ]
anc.3 <- log2( orth$l[,'danio_rerio'] )
anc.3 <- anc.3[ anc.3 > log2(76) & !is.na(anc.3) ]

desc.3 <- lapply( 1:100, function(i){
    evol.size(anc.2, 20, anc.2,
              ev.exp.1=quote(rnorm(length(anc), -0.15, 0.75)),
              ev.exp.2=quote(rnorm(length(anc), -0.15, 0.75)))
})

desc.4 <- lapply( 1:100, function(i){
    evol.size(anc.3, 20, anc.3,
              ev.exp.1=quote(rnorm(length(anc), -0.15, 0.75)),
              ev.exp.2=quote(rnorm(length(anc), -0.15, 0.75)))
})

## this has an issue
d1.KL <- lapply(2:ncol(d1$d1), function(i){ KL.matrix(d1$d1[,i], d1$d2[,i]) })
with(d1.KL[[20]], image.marginal( b1, b2, KL, dist.w=0.25, dist.h=0.25 ))

desc.stats <- lapply(desc, function(d1){
    ## how about correlations and the others?
    t(sapply( 2:21, function(i){
        x <- d1$d1[,i]
        y <- d1$d2[,i]
        ## 
        brk.l <- 4
        x.brk <- seq(min(x), max(x), length.out=brk.l)
        y.brk <- seq(min(y), max(y), length.out=brk.l)
        x.l <- cut(x, breaks=x.brk, include.lowest=TRUE)
        y.l <- cut(y, breaks=y.brk, include.lowest=TRUE)
        cor.x <- tapply( 1:length(x), x.l, function(i){ cor(x[i], y[i]) })
        cor.y <- tapply( 1:length(y), y.l, function(i){ cor(x[i], y[i]) })
        mi.x <-  tapply( 1:length(x), x.l, function(i){ mutual.info(x[i], y[i], nbins1=5) })
        mi.y <-  tapply( 1:length(y), y.l, function(i){ mutual.info(x[i], y[i], nbins1=5) })
        c(cor.x, cor.y, mi.x, mi.y)
    }))
})

desc.stats.2 <- lapply(desc, function(d1){
    ## how about correlations and the others?
    t(sapply( 2:21, function(i){
        x <- d1$d1[,i]
        y <- d1$d2[,i]
        ##
        b <- x <= quantile(x, 0.99) & y <= quantile(y, 0.99)
        x <- x[b]
        y <- y[b]
        brk.l <- 3
        x.brk <- seq(min(x), max(x), length.out=brk.l)
        y.brk <- seq(min(y), max(y), length.out=brk.l)
        x.l <- cut(x, breaks=x.brk, include.lowest=TRUE)
        y.l <- cut(y, breaks=y.brk, include.lowest=TRUE)
        cor.x <- tapply( 1:length(x), x.l, function(i){ cor(x[i], y[i]) })
        cor.y <- tapply( 1:length(y), y.l, function(i){ cor(x[i], y[i]) })
        mi.x <-  tapply( 1:length(x), x.l, function(i){ mutual.info(x[i], y[i], nbins1=10) })
        mi.y <-  tapply( 1:length(y), y.l, function(i){ mutual.info(x[i], y[i], nbins1=10) })
        c(cor.x, cor.y, mi.x, mi.y)
    }))
})

## take a random descendant from the two different pools
## of
desc.2.stats <- lapply(1:length(desc.2), function(i){
    d1 <- desc.2[[i]]$d1
    x.i <- sample(10:ncol(d1), size=1)
    y.i <- sample(1:length(desc), size=10)
    x <- desc.2[[i]]$d1[,x.i]
    sapply( y.i, function(j){
        ## here take only the most evolved one..
        y <- desc[[j]]$d1[,21]
                b <- x <= quantile(x, 0.99) & y <= quantile(y, 0.99)
        x <- x[b]
        y <- y[b]
        brk.l <- 3
        x.brk <- seq(min(x), max(x), length.out=brk.l)
        y.brk <- seq(min(y), max(y), length.out=brk.l)
        x.l <- cut(x, breaks=x.brk, include.lowest=TRUE)
        y.l <- cut(y, breaks=y.brk, include.lowest=TRUE)
        cor.all <- cor(x, y)
        mi.all <- mutual.info(x, y, nbins1=10)
        cor.x <- tapply( 1:length(x), x.l, function(i){ cor(x[i], y[i]) })
        cor.y <- tapply( 1:length(y), y.l, function(i){ cor(x[i], y[i]) })
        mi.x <-  tapply( 1:length(x), x.l, function(i){ mutual.info(x[i], y[i], nbins1=10) })
        mi.y <-  tapply( 1:length(y), y.l, function(i){ mutual.info(x[i], y[i], nbins1=10) })
        c(cor.x, cor.y, mi.x, mi.y, cor.all, mi.all)
    })
})

d2.stats <- t(do.call(cbind, desc.2.stats))

## take a random descendant from the two different pools
## of
desc.3.stats <- lapply(1:length(desc.3), function(i){
    d1 <- desc.3[[i]]$d1
    x.i <- sample(10:ncol(d1), size=1)
    y.i <- sample(1:length(desc.3), size=10)
    x <- desc.3[[i]]$d1[,x.i]
    sapply( y.i, function(j){
        ## here take only the most evolved one..
        y <- desc.3[[j]]$d1[,21]
        b <- x <= quantile(x, 0.99) & y <= quantile(y, 0.99)
        x <- x[b]
        y <- y[b]
        brk.l <- 3
        x.brk <- seq(min(x), max(x), length.out=brk.l)
        y.brk <- seq(min(y), max(y), length.out=brk.l)
        x.l <- cut(x, breaks=x.brk, include.lowest=TRUE)
        y.l <- cut(y, breaks=y.brk, include.lowest=TRUE)
        cor.all <- cor(x, y)
        mi.all <- mutual.info(x, y, nbins1=10)
        cor.x <- tapply( 1:length(x), x.l, function(i){ cor(x[i], y[i]) })
        cor.y <- tapply( 1:length(y), y.l, function(i){ cor(x[i], y[i]) })
        mi.x <-  tapply( 1:length(x), x.l, function(i){ mutual.info(x[i], y[i], nbins1=10) })
        mi.y <-  tapply( 1:length(y), y.l, function(i){ mutual.info(x[i], y[i], nbins1=10) })
        c(cor.x, cor.y, mi.x, mi.y, cor.all, mi.all)
    })
})

d3.stats <- t(do.call(cbind, desc.3.stats))

desc.4.stats <- lapply(1:length(desc.4), function(i){
    d1 <- desc.4[[i]]$d1
    x.i <- sample(10:ncol(d1), size=1)
    y.i <- sample(1:length(desc.4), size=10)
    x <- desc.4[[i]]$d1[,x.i]
    sapply( y.i, function(j){
        ## here take only the most evolved one..
        y <- desc.4[[j]]$d1[,21]
        b <- x <= quantile(x, 0.99) & y <= quantile(y, 0.99)
        x <- x[b]
        y <- y[b]
        brk.l <- 3
        x.brk <- seq(min(x), max(x), length.out=brk.l)
        y.brk <- seq(min(y), max(y), length.out=brk.l)
        x.l <- cut(x, breaks=x.brk, include.lowest=TRUE)
        y.l <- cut(y, breaks=y.brk, include.lowest=TRUE)
        cor.all <- cor(x, y)
        mi.all <- mutual.info(x, y, nbins1=10)
        cor.x <- tapply( 1:length(x), x.l, function(i){ cor(x[i], y[i]) })
        cor.y <- tapply( 1:length(y), y.l, function(i){ cor(x[i], y[i]) })
        mi.x <-  tapply( 1:length(x), x.l, function(i){ mutual.info(x[i], y[i], nbins1=10) })
        mi.y <-  tapply( 1:length(y), y.l, function(i){ mutual.info(x[i], y[i], nbins1=10) })
        c(cor.x, cor.y, mi.x, mi.y, cor.all, mi.all)
    })
})

d4.stats <- t(do.call(cbind, desc.4.stats))

par(mfrow=c(2,4))
hist( d3.stats[,1], breaks=100 )
hist( d3.stats[,2], breaks=100 )
hist( d3.stats[,3], breaks=100 )
hist( d3.stats[,4], breaks=100 )

hist( dr.stats.2[,1] )
hist( dr.stats.2[,2] )
hist( dr.stats.2[,3] )
hist( dr.stats.2[,4] )

hist( d4.stats[,2] / d4.stats[,1] )
hist( d4.stats[,4] / d4.stats[,3] )
hist( d4.stats[,6] / d4.stats[,5] )
hist( d4.stats[,8] / d4.stats[,7] )

hist( dr.stats.2[,2] / dr.stats.2[,1] )
hist( dr.stats.2[,4] / dr.stats.2[,3] )
hist( dr.stats.2[,6] / dr.stats.2[,5] )
hist( dr.stats.2[,8] / dr.stats.2[,7] )

## And the ratios:
plot(d2.stats[,9], d2.stats[,2] / d2.stats[,1])
plot(d2.stats[,9], d2.stats[,4] / d2.stats[,3])
plot(d2.stats[,9], d2.stats[,6] / d2.stats[,5])
plot(d2.stats[,9], d2.stats[,8] / d2.stats[,7])

plot(dc.stats.2[,9], dc.stats.2[,2] / dc.stats.2[,1])
plot(dc.stats.2[,9], dc.stats.2[,4] / dc.stats.2[,3])
plot(dc.stats.2[,9], dc.stats.2[,6] / dc.stats.2[,5])
plot(dc.stats.2[,9], dc.stats.2[,8] / dc.stats.2[,7])

plot(bs.stats.2[,9], bs.stats.2[,2] / bs.stats.2[,1])
plot(bs.stats.2[,9], bs.stats.2[,4] / bs.stats.2[,3])
plot(bs.stats.2[,9], bs.stats.2[,6] / bs.stats.2[,5])
plot(bs.stats.2[,9], bs.stats.2[,8] / bs.stats.2[,7])

plot(dr.stats.2[,9], dr.stats.2[,2] / dr.stats.2[,1])
plot(dr.stats.2[,9], dr.stats.2[,4] / dr.stats.2[,3])
plot(dr.stats.2[,9], dr.stats.2[,6] / dr.stats.2[,5])
plot(dr.stats.2[,9], dr.stats.2[,8] / dr.stats.2[,7])

plot(d3.stats[,9], d3.stats[,2] / d3.stats[,1], xlim=range(dr.stats.2[,9]))
plot(d3.stats[,9], d3.stats[,4] / d3.stats[,3], xlim=range(dr.stats.2[,9]))
plot(d3.stats[,9], d3.stats[,6] / d3.stats[,5], xlim=range(dr.stats.2[,9]))
plot(d3.stats[,9], d3.stats[,8] / d3.stats[,7], xlim=range(dr.stats.2[,9]))

plot(dr.stats.2[,9], dr.stats.2[,2] / dr.stats.2[,1])
plot(dr.stats.2[,9], dr.stats.2[,4] / dr.stats.2[,3])
plot(dr.stats.2[,9], dr.stats.2[,6] / dr.stats.2[,5])
plot(dr.stats.2[,9], dr.stats.2[,8] / dr.stats.2[,7])

plot(d4.stats[,9], d4.stats[,2] / d4.stats[,1], xlim=range(dr.stats.2[,9]))
plot(d4.stats[,9], d4.stats[,4] / d4.stats[,3], xlim=range(dr.stats.2[,9]))
plot(d4.stats[,9], d4.stats[,6] / d4.stats[,5], xlim=range(dr.stats.2[,9]))
plot(d4.stats[,9], d4.stats[,8] / d4.stats[,7], xlim=range(dr.stats.2[,9]))

plot(dr.stats.2[,9], dr.stats.2[,2] / dr.stats.2[,1])
plot(dr.stats.2[,9], dr.stats.2[,4] / dr.stats.2[,3])
plot(dr.stats.2[,9], dr.stats.2[,6] / dr.stats.2[,5])
plot(dr.stats.2[,9], dr.stats.2[,8] / dr.stats.2[,7])

## these are massive differences, so I think we can have some confidence
## in the data.
