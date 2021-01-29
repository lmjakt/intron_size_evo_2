## read in the tables for species counts and membership of vertebrates
## and teleosts.

act.s <- read.table( "Actinopterygii_order_sizes.tsv", stringsAsFactors=FALSE, sep="\t" )
chord.s <- read.table( "Chordata_class_sizes.tsv", stringsAsFactors=FALSE, sep="\t" )

colnames(act.s) <- c('id', 'taxon', 'rank', 'n')
colnames(chord.s) <- c('id', 'taxon', 'rank', 'n')

tel.m <- read.table( "teleost_membership.tsv", stringsAsFactors=FALSE, sep="\t" )
vert.m <- read.table( "vertebrate_membership.tsv", stringsAsFactors=FALSE, sep="\t" )

colnames(tel.m) <- c('taxon', 'member')
colnames(vert.m) <- c('taxon', 'member')

all(act.s$taxon == tel.m$taxon) ## TRUE
all(chord.s$taxon == vert.m$taxon) ## FALSE, but only because of Mammalia Linneaus..

## we can merge the membership data into these
act.s <- cbind(act.s, 'teleost'=tel.m$member)
chord.s <- cbind(chord.s, 'vertebrate'=vert.m$member)

## We have two unknown clades (teleost = -1)
## Cetomimiformes       33
##   In the catalogue of life (COL) it's descendants are:
##   Barbourisiidae, Cetomimidae, Rondeletiidae
##   These are all Teleosts, so we should add 33 to the total count
##   of teleosts
## Saccopharyngiformes  28
##   Descendants in COL:  Cyematidae,  Eurypharyngidae, Monognathidae, Saccopharyngidae
##   These are also teleosts, so I should correct the values in the merged
##   table.
##
## This also suggests how the script can be improved to define membership
## automatically.

act.s[ act.s$taxon == 'Cetomimiformes', 'teleost' ] <- 1
act.s[ act.s$taxon == 'Saccopharyngiformes', 'teleost' ] <- 1

tel.total <- with(act.s, sum(n * (teleost == 1)))
vert.total <- with(chord.s, sum(n * (vertebrate == 1)))

tel.total / vert.total ## 0.474223

## i.e. almost half..
## We can format this as a markdown table?

make.table <- function( df, sep.width=2 ){
    head.widths <- nchar( colnames( df ))
    col.widths <- apply( df, 2, function(x){ max(nchar(x))} )
    col.widths <- pmax( head.widths, col.widths ) + 1
    separator <- paste( rep(" ", sep.width), collapse="" )
    blanks <- sapply(col.widths, function(x){ paste(rep(" ", x), collapse="") })
    sep.row <- paste( sapply(col.widths, function(x){ paste(rep("-", x), collapse="") }),
                     collapse=separator )
    header <- sapply( 1:ncol(df), function(i){
        x <- colnames(df)[i]
        substring( blanks[i], 1, nchar(x) ) <- x
        blanks[i]
    })
    header <- paste(header, collapse=separator)
    rows <- apply(df, 1, function(row){
        y <- sapply( 1:ncol(df), function(i){
        x <- row[i]
        substring( blanks[i], 1, nchar(x) ) <- x
        blanks[i]
        })
        paste(y, collapse=separator)
    })
    c(header,  sep.row, rows)
}

chord.md <- make.table( chord.s[ order(chord.s$n, decreasing=TRUE), ] )
act.md <- make.table( act.s[ order(act.s$n, decreasing=TRUE), ] )

writeLines( chord.md, "chordate_vertebrate_counts.md" )
writeLines( act.md, "actinopterygii_teleost_counts.md" )
