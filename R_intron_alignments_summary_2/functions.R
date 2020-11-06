
sp.name <- function(x, from=".", to="_", toCap=TRUE, toLow=FALSE){
    x <- sub(from, to, x, fixed=TRUE)
    if(toCap)
        substr(x, 1,1) <- toupper( substr(x, 1,1))
    if(toLow)
        substr(x, 1,1) <- tolower( substr(x, 1,1))
    x
}

to.sp <- function(x){
    x <- sub("[._]", " ", x)
    substr(x,1,1) <- toupper( substr(x,1,1) )
    x
}

read.ortho <- function(path="../R_172_genomes/", conv.sp.names=FALSE){
    files <- paste(path, c('dr_intron_fam.txt', 'dr_intron_orthology_id.txt',
                           'dr_intron_orthology_tr.txt', 'dr_intron_orthology_i.txt',
                           'dr_intron_orthology_l.txt'), sep="")
    orth <- lapply(files, read.table, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    if(conv.sp.names)
        for(i in 2:length(orth))
            colnames(orth[[i]]) <- sp.name( colnames(orth[[i]]) )
    
    names(orth) <- c('fam', 'id', 'tr', 'i', 'l')
    orth
}

## log transform with pseudocount
log.tf <- function(x, pc=1){ log2(x + pc) }

score.col <- function(scores, min.sc=min(scores), max.sc=max(scores)){
    cv <- (scores-min.sc) / (max.sc - min.sc)
    cols <- hsv(1, s=cv, v=0.8)
    list('cols'=cols, 'range'=c(min.sc, max.sc))
}
    
## 
plot.alignments <- function(al, col.pid=FALSE, text.labels=FALSE, sp.col=sp.col.3,
                            col.f=score.col, col.tform=eval, rect.border=NA
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
    plot.window(xlim=xlim, ylim=c(y.min, 2 + length(m.scores)), xaxs='i', yaxs='i' )
##    cv <- coords[,'score'] / max(m.scores)
    if(col.pid)
        cols <- col.f( col.tform(unlist(sapply( o, function(i){ p.id( al$al[[i]]$seq ) }))) )
    else
        cols <- col.f( col.tform(coords[,'score']) )
##        cv <- unlist(sapply( o, function(i){ p.id( al$al[[i]]$seq ) }))
##    cols <- hsv(1, s=cv, v=0.8)
    rect( coords[,'a_beg'], coords[,'y'], coords[,'a_end'], coords[,'y']+0.8,
         col=cols$cols, border=rect.border )
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


score.tform.l2 <- function(x, inverse=FALSE){
    if(inverse)
        return(2^x)
    log2(x)
}

## this uses functions defined in ../R_intron_alignments_2/functions.R
## the dependancies here are a complete mess.
alignment.figure <- function(al, sp.c=sp.col, class.c=class.col,
                             c.tform=score.tform.l2, c.f=score.col,
                             id.tbl=orth$id, tr.tbl=orth$tr,
                             i.tbl=orth$i, l.tbl=orth$l, class.cex=1.25,
                             layout.w=c(0.12, 0.88)){
    id <- id.tbl[ al$i, al$sp1 ]
    tr <- tr.tbl[ al$i, al$sp1 ]
    int.i <- i.tbl[ al$i, al$sp1 ]
    int.l <- l.tbl[ al$i, al$sp1 ]
    layout(matrix(2:1, nrow=1), widths=layout.w)
    par(mar=c(5.1, 0.1, 4.1, 2.1))
    tmp <- plot.alignments(al, col.f=c.f, col.tform=c.tform, sp.col=sp.c)
    ##
    mtext( sprintf("%s : %s intron no: %d length: %d", id, tr, int.i, int.l), cex=1.2, line=1 )
    ## with(par(), text(usr[1], usr[3], sprintf("%d alignments above max size (%1.1e)",
    ##                                          sum( as.double(al$l1) * as.double(al$l2) > max.l ), max.l ),
    ##                  adj=c(0, 0)))
    mtext(sprintf("%d alignments above max size (%1.1e)",
                  sum( as.double(al$l1) * as.double(al$l2) > max.l ), max.l ),
          side=1, line=2.5, at=0 )
    if(nrow(tmp$coord) == 0){
        return()
    }
    par(mar=c(5.1, 1.1, 4.1, 0.1))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,100))
    text(0.05, 100 - 1.5 * strheight("A", cex=class.cex) * 1:length(class.c),
         capitalise(names(class.c)), col=class.c, adj=c(0,0), font=2, cex=class.cex)
    hm.v <- seq( tmp$cols$range[1], tmp$cols$range[2], length.out=100 )
    hm.y <- seq( 0, 25, length.out=length(hm.v))
    hm.yd <- diff(hm.y)[1]
    rect( 0.85, hm.y, 1.0, hm.y + hm.yd, col=c.f(hm.v)$cols, border=NA )
    lab.i <- seq(1, length(hm.v), length.out=5)
    text( 0.82, hm.y[ lab.i ] + hm.yd/2, sprintf("%.1e", c.tform(hm.v[lab.i], inverse=TRUE)), adj=c(1,0.5),
         cex=0.85)
}

capitalise <- function(lab){
    substr(lab, 1, 1) <- toupper( substr(lab, 1, 1))
    lab
}

read.fa.ids <- function(fn){
    tmp <- readLines(fn)
    des <- tmp[ seq(1, length(tmp), 2) ]
    ids <- unname(sapply(des, function(x){ sub("^>([^ ]+).+", "\\1", x) }))
    introns <- as.numeric(sapply( strsplit(ids, '_'), function(x){ x[2] }))
    genes <- unname(sapply( strsplit(ids, '_'), function(x){ x[1] }))
    des.terms <- strsplit( des, " " )
    tr <- unname( sapply( des.terms, function(x){ x[2] }))
    sp1 <- unname( sapply( des.terms, function(x){ x[5] }))
    sp2 <- unname( sapply( des.terms, function(x){ x[6] }))
    coords <- t(sapply( des.terms, function(x){ as.numeric(x[7:12]) }))
    colnames(coords) <- c('a_beg', 'a_end', 'b_beg', 'b_end', 'score', 'length')
    ## we should not need the sequences.. 
    seq <- tmp[ seq(2, length(tmp), 2) ]
    list(id=ids, gene=genes, intron=introns, des=des,
         seq=seq, tr=tr, sp1=sp1, sp2=sp2, coords=coords)
}

## summarise the alignment information. This is similar to what I made for
## the extract.top.scores etc, but somewhat more reasonable I hope

## takes the 
al.summarise <- function( al, sp.class ){
    x <- lapply( colnames(sp.class), function(cl){
        b <- sp.class[ al$sp, cl ]
        if(!sum(b))
            return(NULL)
        null.al.b <- sapply( al$al, function(x){ is.null(x) || is.null(x$pos) } )
        b <- b & (!null.al.b)
        if(!sum(b))
            return(NULL)
        pos <- t( sapply( al$al[b], function(p){ p$pos[1,] }) )
        sp1 <- al$sp1
        sp2 <- al$sp[b]
        seq <- lapply( al$al[b], function(p){ p$seq[[1]] })
        l2 <- al$l2[b]
        o <- order( pos[,'score'], decreasing=TRUE )
        j <- ifelse(is.null(al$j), al$i, al$j)
        list('i'=al$i, 'j'=j, 'n'=sum(b),
             'l1'=al$l1, 'l2'=l2[o], 'sp1'=sp1, 'sp2'=sp2[o], 'pos'=pos[o,,drop=FALSE], 'seq'=seq[o] )
    })
    names(x) <- colnames(sp.class)
    x
}

## al: a list of top alignment structures obtained using al.summarise
## cl: the taxonomic class specified
## prefix: a file name prefix
export.al.seq <- function( al, cl, prefix, dir=".", orth ){
    fname <- paste(dir, "/", prefix, "_", cl, ".fasta", sep="")
    lines <- vector(mode='list', length=length(al))
    for(i in 1:length(lines)){
        if(is.null(al[[i]]) || is.null( al[[i]][[cl]]))
            next
        x <- al[[i]][[cl]]
        lines[[i]] <- paste(">", orth$tr[ x$j, x$sp1 ], "_", x$j, "_", i, "\t",
                            x$sp1, " ", x$sp2[1], "\n", degap.seq( x$seq[[1]][1]), sep="")
    }
    writeLines( do.call(c, lines), fname )
    fname
}

## bl is the blast data obtained from sequences in al.top
## al.top is a list of structures obtained using al.summarise
## cl is the species class
query.covs <- function(bl, al.top, cl, exclude.sub=NULL){
    covs <- vector(mode='list', length=length(al.top))
    q.i <- sub(".+?_[^_]+_([^_]+)$", "\\1", bl$qseqid )
    bl.j <- tapply( 1:length(q.i), q.i, eval)
    ## a tapply is probably the fastest here; but I prefer to have the
    ## the resulting data in the specified order
    ## so we do it the slow way.
    for(i in 1:length(al.top)){
        ql <- 0
        if(!is.null(al.top[[i]]) && !is.null(al.top[[i]][[cl]]))
            ql <- nchar(degap.seq( al.top[[i]][[cl]]$seq[[1]][1] ))
        j <- bl.j[[ as.character(i) ]]
        if(!is.null( exclude.sub ))
            j <- j[ !grepl(exclude.sub, bl[j, 'sseqid']) ]
        if(any( bl[j, 'qlen'] != ql ))
            stop(paste(i, "Multiple query lengths for a single row"))
        covs[[i]] <- rep(0, ql)
        for(k in j){
            r <- as.integer(bl[k, c('qstart', 'qend')])
            covs[[i]][ r[1]:r[2] ] <- covs[[i]][ r[1]:r[2] ] + 1
        }
    }
    covs
}

get.align.pars <- function(al.top, sm=NULL, cl, orth=orth.sc){
    if(!is.null(sm))
        al <- al.top[[sm]]
    else
        al <- al.top
    null.al <- rep(NA, 13)
    tmp <- t(sapply( al, function(x){
        if(is.null(x[[cl]]))
            return(null.al)
        y <- x[[cl]]
        as.double( c(y$i, y$j, y$pos[1,], y$l1, y$l2[1], sum(y$l2), orth$score[ y$j, y$sp2[1] ],
                     which(colnames(orth$score) == y$sp2[1]) ))
    }))
    colnames(tmp) <- c('i', 'j', 'a_beg', 'a_end', 'b_beg', 'b_end', 'score', 'length',
                       'l1', 'l2', 'l2.s', 'i.score', 'sp2')
    tmp
}

plot.al.scores <- function(scores, bl.cov.q, cl, sm, qnt='50%', do.plot=TRUE, min.al.score=NULL,
                           ...){
    q <- bl.cov.q[[cl]][[sm]][,qnt]
    b <- q > 0 & q < 3
    na <- is.na(scores[[sm]][[cl]][,'i'])
    sc <- scores[[sm]][[cl]]
    if(!is.null(min.al.score))
        b <- b & sc[,'i.score'] > min.al.score
    x <- log2(sc[,'l1'] * sc[,'l2.s'])
    y <- log2(sc[,'score'])
    if(do.plot)
        plot(x[b], y[b],  ...)
    invisible(list(x=x, y=y, b=b, na=na))
}

lm.res <- function(pts, mod, res.q=seq(0, 1, 0.05), use.b=TRUE ){
    m.y <- mod$coefficients[1] + pts$x * mod$coefficients[2]
    res <- pts$y - m.y
    if(use.b)
        res.qnt <- quantile(res[pts$b], probs=res.q, na.rm=TRUE)
    else
        res.qnt <- quantile(res, probs=res.q, na.rm=TRUE)
    list(pts=pts, my=m.y, res=res, qnt=res.qnt)
}

## cl and sm refer to class and sample respectively
plot.lm <- function(cls, sms, al.pars, bl.cov.q, qnt.pc='50%', min.al.score=NULL){
    pl.lm <- lapply( cls, function(cl){
        tmp <- lapply( sms, function(sm){
            pts <- plot.al.scores( al.pars, bl.cov.q, cl, sm, qnt=qnt.pc, min.al.score=min.al.score,
                                  xlab='log2 mn', ylab='score', main=paste(sm, cl))
            mod <- with(pts, lm( y[b] ~ x[b] ))
            list('lm'=mod, 'pts'=pts)
        })
        names(tmp) <- sms
        tmp
    })
    names(pl.lm) <- cls
    pl.lm
}



plot.points <- function(pts, model, model.col='red', model.lwd=2, ...){
    b <- pts$b
    plot(pts$x[b], pts$y[b], ...)
    abline(model, col=model.col, lwd=model.lwd)
}

parse.bl.names <- function(fn){
    bn <- sub(".*?([^/]+)$", "\\1",  fn)
    sm <- sub("([^_]+)_.+", "\\1", bn)
    cl <- sub("[^_]+_([^_]+)_.+", "\\1", bn)
    sp <- sub(".+?_([^_]+_[^_]+)\\.bl", "\\1", bn)
    substr( sp, 1, 1 ) <- tolower( substr( sp, 1, 1 ))
    sp <- sub("_", ".", sp)
    data.frame(fn=fn, bn=bn, sm=sm, cl=cl, sp=sp, stringsAsFactors=FALSE)
}

### Two horrible functions for plotting and returning coverage statistics.
### These make use of lots of global variables and should be rewritten.
plot.transcript.region <- function(orth.row, species, seed.coords, bl, query.id, alignments, max.nm=1e8,
                                   ex.height=0.4, bl.height=0.3){
    tr <- orth$tr[ orth.row, species ]
    int.i <- orth$i[ orth.row, species ]
    db <- dbs[ species ]
    tr.exons <- get.transcript.exons( db.cred, db, tr )
    tr.description <- get.transcript.description( db.cred, db, tr )
    if(!nrow(tr.exons)){
        plot.new()
        plot.window(xlim=c(0,1), ylim=c(0,1))
        mtext(species)
        plot.new()
        return()
    }
    strand = tr.exons[1,'strand']
    int.i.beg <- ifelse(strand == 1, tr.exons[ int.i, 'ex_end'], -tr.exons[ int.i, 'ex_start'] )
    int.i.end <- ifelse(strand == 1, tr.exons[ int.i+1, 'ex_start'], -tr.exons[ int.i+1, 'ex_end'] )
    xlim <- range( strand * tr.exons[,c('ex_start', 'ex_end')] )
    plot.new()
    plot.window(xlim=xlim, ylim=c(0,4))
    lines( xlim, c(1,1) )
    rect( strand * tr.exons[,'ex_start'], 1-ex.height, strand * tr.exons[,'ex_end'], 1 + ex.height )
    segments( int.i.beg, 1, int.i.end, 1, lwd=4 )
    ## lets draw the blast matches that fit within this directory
    bl <- bl[ bl$qseqid == query.id & sub("gb\\|([^|]+).+", "\\1", bl$sseqid) == tr.exons$name[1], ]
    if(nrow(bl))
        rect( strand * bl$sstart, 2-bl.height, strand * bl$send, 2+bl.height, density=20 )
    axis(1)
    ## lets try to draw the alignment that we have from the initial alignments.
    ## I don't really need to have a special case for the seed species here; but
    ## 
    species.al <- tolower(sub("_", ".", species))
    align <- alignments[[ as.character(orth.row) ]]
    sp1 <-align$sp1 ## sub(".", "_", align$sp1, fixed=TRUE)
    ## substr(sp1, 1,1) <- toupper(substr(sp1, 1,1))
    if(length(align)){
        if(align$sp1 == species.al){
            segments( int.i.beg + seed.coords[1], 1, int.i.beg + seed.coords[2], 1, col='red', lwd=4 ) 
        }else{
            al.i <- which( align$sp == species.al )
            if(length(al.i) == 1){
                al.pos <- align$al[[al.i]]$pos[1,]
                if( as.double(align$l1) * as.double(align$l2[al.i]) > max.nm )
                    segments( int.i.beg, 1, int.i.end, 1, lwd=4, col=rgb(0.6, 0.6, 0.6) )
                if(!is.null( al.pos ))
                    segments( int.i.beg + al.pos[3], 1, int.i.beg + al.pos[4], 1, col='red', lwd=4 )
            }
        }
    }
    mtext(species, font=2)
    with(par(), text(usr[1], usr[4],
                     paste( tr.description[,c('gene.id', 'gene', 'description')], collapse=" "),
                     adj=c(0,2) ))
    ## we are now repeating ourselves a bit, but, let's have a go at this
    plot.new()
    if(is.na( orth$l[orth.row, species] ) || is.na(orth$l[orth.row, sp1])){
        plot.window(xlim=c(0,1), ylim=c(0,1))
        return()
    }
    plot.window( xlim=c(0, orth$l[orth.row, sp1]), ylim=c(0, orth$l[orth.row, species]), asp=1 )
    al.i <- which( align$sp == species.al )
    segments( 1, 1, c(orth$l[orth.row, sp1], 1), c(1, orth$l[orth.row, species]), lwd=2 )
    lines( seed.coords[1:2], c(1,1), lwd=2, col='red' )
    axis(1)
    axis(2)
    if(length(al.i) == 1){
        al.pos <- align$al[[al.i]]$pos
        if(!is.null( al.pos ))
            segments( al.pos[,'a_beg'], al.pos[,'b_beg'], al.pos[,'a_end'], al.pos[,'b_end'], col='black', lwd=1 )
    }
}


## take a set of fasta data entry as an argument and then open the appropriate blast files
## and all sorts of global variables
plot.coverages <- function(fa.index, sp.b, alignments=aligns$long, start.from=1,
                           interactive=TRUE, do.plot=TRUE){
    ## a function that returns a depth of coverage vector
    al.cov <- function(bl, id, l, max.e=1e-5, min.pid=65, exclude.sub=NULL,
                       tform=log.tf){
        bl <- bl[ bl$qseqid == id & bl$qlen == l, ]
        if(!is.null(exclude.sub))
            bl <- bl[ !grepl(exclude.sub, bl$sseqid), ]
        cov <- rep(0, l)
        if(!nrow(bl))
            return(cov)
        for(i in 1:nrow(bl)){
            r <- bl[i,'qstart']:bl[i,'qend']
            cov[r] <- cov[r] + 1
        }
        tform(cov)
    }
    plot.cov <- function(cov, main=NULL, poly.col=rgb(0.7, 0.7, 0.7), ...){
        pts.polygon <- function(y, y.min=0){
            x <- c(1, 1:length(y), length(y))
            y <- c(0, y, 0)
            cbind(x=x, y=y)
        }
        ylim=c(0, max(cov)) * 2
        if(ylim[2] == 0)
            ylim[2] <- 1
        if(is.null(nrow(cov))){
            plot(1:length(cov), cov, type='n', ylim=ylim, main=main, xaxs='i', yaxs='i',
                 frame.plot=FALSE, yaxt='n', ...)
            polygon(pts.polygon(cov), col=poly.col)
        }else{
            plot(1:nrow(cov), 1:nrow(cov), type='n', ylim=ylim, main=main, xaxs='i', yaxs='i',
                 frame.plot=FALSE, yaxt='n', ...)
            for(i in 1:ncol(cov))
                polygon( pts.polygon(cov), col=poly.col)
        }
        axis(side=2, at=seq(0, max(cov), ceiling( max(cov) / 4 )))
    }
##
    cov.qnts.l <- vector(mode='list', length=length(fa.index))
    cov.qnts.i <- 0
    for(fa.i in fa.index){
        cov.qnts.i <- cov.qnts.i + 1
        fa.d <- query.seq[[fa.i]]
        fa.cl <- fa.files.cl[fa.i]
        fa.sm <- fa.files.sm[fa.i]
        bl.b <- bl.files.cl == fa.cl & bl.files.sm == fa.sm
        bl.i1 <- which(bl.b &  bl.files.sp == "danio.rerio")
        bl.i2 <- which(bl.b &  bl.files.sp != "danio.rerio" &
                       sp.b[ bl.files.sp ])
        qnt.probs <- seq(0, 1, 0.1)
        cov.qnts <- matrix(0, nrow=length(fa.d$id), ncol=length(qnt.probs))
        colnames(cov.qnts) <- names( quantile(1:100, qnt.probs))
        for(i in start.from:length( fa.d$id )){
            sp <- fa.d$sp2[i]
            id <- fa.d$id[i]
            tr <- fa.d$tr[i]
            al.i <- i
            while( (alignments[[al.i]]$seq.meta['tr','danio.rerio'] != tr || orth$i[ alignments[[al.i]]$i, 'danio.rerio'] != fa.d$intron[i]) &&
                   al.i <= length(fa.d$id) )
                al.i <- al.i + 1
            if( al.i > length(fa.d$id)){
                if(do.plot){
                    plot.new()
                    plot.new()
                    plot.new()
                }
                next
            }
            ## this is so ugly; and because I did not include this information
            ## when I exported the data.
            orth.row <- alignments[[al.i]]$i
            tr.2 <- orth$tr[orth.row, bl.files.sp[ bl.i2 ]]
            tr.2.i <- orth$i[orth.row, bl.files.sp[ bl.i2 ]]
            q.l <- nchar(query.seq[[fa.i]]$seq[i])
            cov.1 <- al.cov( bl.data[[bl.i1]], id, q.l, exclude.sub='ALT', tform=eval)
            cov.qnts[i, ] <- quantile( cov.1, qnt.probs, na.rm=TRUE )
            if(!do.plot)
                next
            cov.1 <- log.tf(cov.1)
            cov.2 <- sapply( bl.data[bl.i2], al.cov, id=id, l=q.l )
            par(oma=c(2.1, 2.1, 2.1, 2.1))
            par(mar=c(2.1, 3.1, 2.1, 2.1))
            layout.nr <- 1 + ncol(cov.2)
            layout.m <- cbind( 1:layout.nr, matrix( layout.nr + 1:(2*layout.nr), ncol=2, byrow=TRUE ))
            layout(layout.m, widths=c(0.2, 0.8, 0.2))
            plot.cov(cov.1, main='')
            for(j in 1:ncol(cov.2))
                plot.cov( cov.2[,j], main='') ## bl.files.sp[ bl.i2[j] ] )
            par(mar=c(2.1, 1.1, 2.1, 2.1))
            cat(i, ": ")
            plot.transcript.region( orth.row, 'danio.rerio', fa.d$coords[i,], bl.data[[bl.i1]], id, alignments )
            for(j in bl.i2){
                cat(j, ", ")
                plot.transcript.region( orth.row, bl.files.sp[j], fa.d$coords[i,], bl.data[[j]], id, alignments )
            }
            cat("\n")
            if(interactive){
                inpt <- readline(paste(i, ":"))
                if(inpt == 'q')
                    break
            }
        }
        cov.qnts.l[[ cov.qnts.i ]] <- cov.qnts
    }
    invisible( do.call( rbind, cov.qnts.l ) )
}


