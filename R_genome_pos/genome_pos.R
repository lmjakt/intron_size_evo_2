## Check whether or not long / short introns are positioned differentely across the genome
## This is because it has been shown that long introns are found at regions of low meiotic
## recombination in Drososphila. Recombation frequencies have been calculated in D. rerio,
## though I have not been able to find the raw data for these. However, rates increase
## towards the end of chromosomes, so that I can have a look to see whether there is any
## particular tendencies.

require('entropy')
source('~/R/general_functions.R')

pos.df <- read.table("../exon_intron_pos.csv", sep="\t", stringsAsFactor=FALSE, header=FALSE)

colnames(pos.df) <- c('species', 'database', 'gene', 'tr', 'chr', 'strand', 'range',
                      'e.start', 'e.size', 'i.start', 'i.size')

## let us divide this up into species

pos.sp <- tapply(1:nrow(pos.df), pos.df$species, function(i){ pos.df[i,] })
names(pos.sp) <- sapply(pos.sp, function(x){ x$species[1] })

## we are probably mostly interested in danio.rerio
## but lets make lists of intron sizes and positions..

extract.intron <- function(x){
    tmp <- tapply( 1:nrow(x), x$chr, function(i){
        pos <- as.numeric(unlist( strsplit( x$i.start[i], "," ) ))
        size <- as.numeric(unlist( strsplit( x$i.size[i], "," ) ))
        list('chr'=x$chr[i[1]], pos=pos, size=size)
    })
    names(tmp) <- sapply(tmp, function(x){ x$chr })
    o <- order( sapply(tmp, function(x){ length(x$pos) }), decreasing=TRUE )
    tmp[o]
}

pos.intron <- lapply( pos.sp, extract.intron )

v2col <- function(v, cm){
    vv= v - min(v, na.rm=TRUE)
    val <- (cm + vv)/(cm + max(vv, na.rm=TRUE))
    hsvScale(v, val=val)
}

## OK. We will make a plot with these types of values:
## (cm is a colour moderator that moderates the colour value)
## x and y are the breaks of the rectangles, (x[-1] == x2)
## (inner) margins are bottom, left, top, right
rect_image <- function(m, x=0:ncol(m), y=0:nrow(m),
                       col.f=v2col, col.args=NULL, margins=rep(0, 4),
                       border=NA, clear=TRUE, axis=TRUE,
                       xlab=NULL, ylab=NULL, lab.cex=1.5){
    x1 <- x[-length(x)]
    x2 <- x[-1]
    y1 <- y[-length(y)]
    y2 <- y[-1]
    if(clear)
        plot.new()
    xr <- range(x)
    yr <- range(y)
    xr <- xr + diff(xr) * c(-margins[2], margins[4])
    yr <- yr + diff(yr) * c(-margins[1], margins[3])
    plot.window(xlim=xr, ylim=yr, xaxs='i', yaxs='i')
    ## and then make x1, x2, y1, and y2 the appropriate matrices
    x1 <- matrix(rep(x1, nrow(m)), nrow=nrow(m), byrow=TRUE)
    x2 <- matrix(rep(x2, nrow(m)), nrow=nrow(m), byrow=TRUE)
    y1 <- matrix(rep(y1, ncol(m)), nrow=nrow(m))
    y2 <- matrix(rep(y2, ncol(m)), nrow=nrow(m))
    if(is.null(col.args))
        cols <- col.f(m)
    else
        cols <- do.call(col.f, c(list(v=m), col.args))
    rect(x1, y1, x2, y2, col=cols, border=border)
    if(axis){
        axis(1)
        axis(2)
    }
    if(!is.null(xlab))
        mtext(xlab, side=1, line=2.5, cex=lab.cex)
    if(!is.null(ylab))
        mtext(ylab, side=2, line=2.5, cex=lab.cex)
}


plot.intron.chr <- function( x, tform=log2 ){
    for(chr in names(x)){
        par(mfrow=c(2,1))
        plot( x[[chr]]$pos, tform(x[[chr]]$size), main=chr, xaxs='i')
        tmp <- discretize2d( tform(x[[chr]]$size), x[[chr]]$pos, numBins1=25, numBins2=100 )
        rect_image( scale(tmp), col.args=list(cm=0.3) )
        inpt <- readline( paste( chr, ":"))
        if(inpt == 'q')
            break
    }
}


plot.intron.chr( pos.intron$danio_rerio )
plot.intron.chr( pos.intron$gasterosteus_aculeatus)
plot.intron.chr( pos.intron$oryzias_latipes)
