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

to.sp <- function(x){
    x <- sub("[._]", " ", x)
    substr(x,1,1) <- toupper( substr(x,1,1) )
    x
}
