## a set of local functions.

## gets a simple classification into teleosts, mammals,
## sauria.
get.sp.taxonomy <- function(){
    teleost.b <- read.table( "../R_172_genomes/teleost_b.txt", header=TRUE, sep="\t" )
    mammal.b <-  read.table( "../R_172_genomes/mammal_b.txt", header=TRUE, sep="\t" )
    sauria.b <-  read.table( "../R_172_genomes/sauria_b.txt", header=TRUE, sep="\t" )
##
    sp.class <- cbind(teleost.b, mammal.b, sauria.b)
    rownames(sp.class) <- sub(" ", ".", rownames(teleost.b))
    colnames(sp.class) <- c('tel', 'mam', 'sau')
    sp.class
}


df2double <- function(x){
    cn <- colnames(x)
    x <- matrix( as.double(x), ncol=ncol(x))
    colnames(x) <- cn
    x
}


## extracts a score table from a single set of alignments
extract.scores <- function(x, max.l=max.l, sp.b=NULL){
    if(!is.null(sp.b))
        sp.b <- sp.b[ x$sp ]
    else
        sp.b <- rep(TRUE, length(x$sp))
    b.i <- which(as.double(x$l1) * as.double(x$l2) <= max.l & sp.b)
    do.call(rbind, lapply(b.i, function(i){
        if(is.null( x$al[[i]]$pos )) return(NULL)
        suppressWarnings( cbind( 'l1'=x$l1, 'l2'=x$l2[i], x$al[[i]]$pos[,'score'] ))
    }))
}

## extracts a score table from a single set of alignments
extract.top.scores <- function(x, max.l=max.l, sp.b=NULL){
    if(!is.null(sp.b))
        sp.b <- sp.b[ x$sp ]
    else
        sp.b <- rep(TRUE, length(x$sp))
    b.i <- which(as.double(x$l1) * as.double(x$l2) <= max.l & sp.b)
    do.call(rbind, lapply(b.i, function(i){
        if(is.null( x$al[[i]]$pos )) return(NULL)
        suppressWarnings( cbind( 'l1'=x$l1, 'l2'=x$l2[i], x$al[[i]]$pos[1,'score'] ))
    }))
}


## from a list; 
extract.scores.l <- function(al.list, max.l, sp.b=NULL ){
    scores <- do.call( rbind,
                      lapply(1:length(al.list), function(i){
                          tmp <- extract.scores(al.list[[i]], max.l=max.l, sp.b=sp.b)
                          if(is.null(tmp))
                              return(NULL)
                          cbind(i, tmp) }))
    scores <- df2double(scores)
    colnames(scores) <- c('i', 'l1', 'l2', 'score')
    cbind(scores, 'mn'=log2(scores[,'l1'] * scores[,'l2']))
}

## from a list
extract.top.scores.l <- function(al.list, max.l, sp.b=NULL ){
    scores <- do.call( rbind,
                      lapply(1:length(al.list), function(i){
                          tmp <- extract.top.scores(al.list[[i]], max.l=max.l, sp.b=sp.b)
                          if(is.null(tmp))
                              return(NULL)
                          cbind(i, tmp) }))
    scores <- df2double(scores)
    colnames(scores) <- c('i', 'l1', 'l2', 'score')
    cbind(scores, 'mn'=log2(scores[,'l1'] * scores[,'l2']))
}

extract.top.al <- function(x, max.l=max.l, sp.b=NULL){
    null.al <- list('coords'=c('score'=0),
                    'seq'=c("",""))
    if(!is.null(sp.b))
        sp.b <- sp.b[ x$sp ]
    else
        sp.b <- rep(TRUE, length(x$sp))
    b.i <- which( as.double(x$l1) * as.double(x$l2) <= max.l & sp.b)
    tmp <- lapply( b.i, function(i){
        if(is.null(x$al[[i]]$pos )) return(null.al)
        list( 'coords'= x$al[[i]]$pos[1,],
              'sp1' = x$sp1, 'sp2'=x$sp[i],
              'seq'= x$al[[i]]$seq[[1]],
             'i'=x$i)
    })
    names(tmp) <- x$sp[b.i]
    o <- order( sapply(tmp, function(y){ y$coords['score'] }),
                decreasing=TRUE )
    tmp[o]
}

## takes two sequences with gaps
degap.seq <- function(seq){
    gsub("-", "", seq)
}

extract.top.al.l <- function(al.list, max.l, sp.b=NULL ){
    als <- lapply(al.list, function(x){
        extract.top.al( x, max.l=max.l, sp.b=sp.b )
    })
}

mnq.class <- function(scores, qnt){
    qnt <- qnt[-length(qnt)]
    colSums( sapply( scores[,'mn'], function(x){ x >= qnt } ))
}

sliced.h <- function(x, slices, breaks){
    tapply(x, slices, hist, breaks=breaks, plot=FALSE)
}

## should really make a function; extract named.. 
plot.counts <- function(h, dens=FALSE, log=TRUE,
                        cols=hsv(1, 0.8, 1:length(h) / length(h)), ...){
    if(dens)
        counts <- sapply(h, function(x){ x$density })
    else
        counts <- sapply(h, function(x){ x$counts })
    if(log)
        counts <- log2( counts )
    mids <- h[[1]]$mids
    plot(mids, counts[,1], type='l', lwd=2, ylim=range(counts[is.finite(counts)], na.rm=TRUE),
         col=cols[1], ...)
    invisible( sapply( 2:ncol(counts[,-1]), function(i){
        lines(mids, counts[,i], lwd=2, col=cols[i])}))
    invisible(cols)
}

score.col <- function(scores, min.sc=min(scores), max.sc=max(scores)){
    cv <- (scores-min.sc) / (max.sc - min.sc)
    cols <- hsv(1, s=cv, v=0.8)
    list('cols'=cols, 'range'=c(min.sc, max.sc))
}
    

plot.alignments <- function(al, col.pid=FALSE, text.labels=FALSE, sp.col=sp.col.3,
                            col.f=score.col, col.tform=eval
                            ){
    ## find out how many of the species have alignments
    ## and order by the highest scoring alignment
    m.scores <- sapply(al$al, function(x){
        ifelse( is.null( nrow(x$pos) ), 0, x$pos[1,'score'] )
    })
    o <- order(m.scores, decreasing=FALSE)
    coords <- sapply(1:length(o), function(i){
        j <- o[i]
        m <- matrix(nrow=0, ncol=4)
        colnames(m) <- c('a_beg', 'a_end', 'score', 'y')
        if(!is.null(al$al[[j]]$pos))
            m <- cbind( al$al[[j]]$pos[,c('a_beg', 'a_end', 'score'), drop=FALSE],
                       'y'=rep(i, nrow(al$al[[j]]$pos)) )
        m
    })
    coords <- do.call(rbind, coords)
    ## and then we can simply call rect with the appropriate parameters
    plot.new()
    y.min <- sum( m.scores == 0 ) + 1
    xlim <- c(-al$l1 * 0.01, al$l1)
    plot.window(xlim=xlim, ylim=c(y.min, 2 + length(m.scores)), xaxs='i' )
##    cv <- coords[,'score'] / max(m.scores)
    if(col.pid)
        cols <- col.f( col.tform(unlist(sapply( o, function(i){ p.id( al$al[[i]]$seq ) }))) )
    else
        cols <- col.f( col.tform(coords[,'score']) )
##        cv <- unlist(sapply( o, function(i){ p.id( al$al[[i]]$seq ) }))
##    cols <- hsv(1, s=cv, v=0.8)
    rect( coords[,'a_beg'], coords[,'y'], coords[,'a_end'], coords[,'y']+0.8,
         col=cols$cols )
    b <- m.scores[o] > 0
    if(text.labels){
        axis(2, at=y.min:length(al$sp) + 0.5, labels=al$sp[o[y.min:length(al$sp)]], las=2)
    }else{
        with(par(), {
            y1 <- y.min:length(al$sp)
            rect(usr[1], y1, xlim[1]/4, y1+1, col=sp.col[ al$sp[o[y.min:length(al$sp)]] ], border=NA )
        })
    }
    axis(1)
    invisible( list( 'coord'=coords, 'm.scores'=m.scores[o], 'sp'=al$sp[o], 'l2'=al$l2[o], 'o'=o,
                    'scores'=coords[,'score'], 'cols'=cols))
}

capitalise <- function(lab){
    substr(lab, 1, 1) <- toupper( substr(lab, 1, 1))
    lab
}
