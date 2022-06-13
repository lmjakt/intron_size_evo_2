source("../R_common_functions/functions.R")
source("~/R/drawinR/drawing_functions.R")
require("RMySQL")

##
## we will use a series of database queries to try and determine the set of transcripts that
## are supported by RNA sequencing from each species.

## read in the database credentials.

db.cred <- read.credentials("db.cred")
## and set the password in the console.
## would be better to encrypt this, but rather too troublesome.

## Some species classification:
sp.class <- readRDS( "../R_trees_distances/sp_class_3.rds")
class.col <-  readRDS( "../R_trees_distances/class_col_3.rds")
sp.col <- readRDS("../R_trees_distances/sp_col_3.rds")
names(sp.col) <- sub("_", ".", names(sp.col))
rownames(sp.class) <- sub("_", ".", rownames(sp.class))

## for plotting we can use a taxonomic order obtained from the
## phylogenetic tree:
sp.y <- readRDS("../R_trees_distances/sp_y.rds")
## unforunately we are using the "." seperation for species names here
names(sp.y) <- sub("_", ".", names(sp.y))

sp.names <- sub("\\.", " ", names(sp.y))
substring(sp.names, 1, 1) <- toupper(substring(sp.names, 1, 1))
names( sp.names ) <- names(sp.y)

names.sp <- names(sp.names)
names(names.sp) <- sp.names

## to also map from the "_"..
sp.names.2 <- sp.names
names(sp.names.2) <- sub("\\.", "_", names(sp.names.2))

## get the orthology; 
orth <- readRDS( "../R_intron_alignments_summary_2/orth.rds")
## this doesn't give use the correct database names.
## but we can get them from:

fam.members <- read.table("../family_members/vertebrate_family_members_1_130.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
## the colnames of fam.members is the set of databases which we want to use.
core.dbs <- colnames(fam.members)
names(core.dbs) <- sub("([^_]+)_([^_]+).+", "\\1.\\2", core.dbs)
all( names(core.dbs) %in% colnames(orth$id) ) ## TRUE

## OK. Then we can start to query the information.
get.transcript.support <- function(db, cred=db.cred){
    query = paste("select distinct a.stable_id, d.logic_name from transcript a",
                  "inner join transcript_intron_supporting_evidence b on a.transcript_id=b.transcript_id",
                  "inner join intron_supporting_evidence c on b.intron_supporting_evidence_id=c.intron_supporting_evidence_id",
                  "inner join analysis d on c.analysis_id=d.analysis_id;")
    db <- ensembl.connect(cred, db)
    tbl <- dbGetQuery(db, query)
    dbDisconnect(db)
    tbl
}

transcript.cred <- lapply(core.dbs, get.transcript.support, cred=db.cred)

## we can make a plot of the total number of transcripts with supporting information
## for these.
par(mar=c(5.1, 8.1, 4.1, 2.1))
sp <- names(sp.y)
y <- barplot( sapply(transcript.cred[sp], nrow), las=2, names.arg="", col=sp.col[sp], xaxs='i', horiz=T )
mtext( sp, side=2, at=y, cex=0.5, las=2)

## So far that is looking good. However, we should actually have a look at what the logic_names are
logic.names <- lapply(transcript.cred, function(x){ table( x[,'logic_name'] ) })

## a manual inspection of a subset of these suggest that all these indicate RNA-seq based evidence for
## introns.

## We now want to estimate the number of transcripts defined by the orthology that have evidence for them

orth.tr.cred <- t(sapply(names(sp.y), function(sp){
    tr <- unique(orth$tr[ !is.na( orth$tr[,sp] ), sp])
    m <- c('n'=length(tr),
           'seq'=sum( tr %in% transcript.cred[[sp]][,'stable_id'] ))
    m <- c(m, 'r'=m[2] / m[1])
}))
rownames(orth.tr.cred) <- names(sp.y)

### 
par(mar=c(5.1, 8.1, 4.1, 2.1))
sp <- names(sp.y)
y <- barplot( orth.tr.cred[,'r.seq'], las=2, names.arg="", col=sp.col[sp], xaxs='i', horiz=T )
mtext( sp, side=2, at=y, cex=0.5, las=2)

### That gives us support from RNA-seq data; as far as I Can tell it is always 
### RNA-seq from the same species that is being annotated.

## We can also make use of the:

## transcript_supporting_feature
## this has links to either
## : dna_align_feature
## : protein_align_feature
## so we need to make two different queries.
## for dna_align_feature it seems that it is predominantly from within species
## sequences; but that information is not present within the database, but we can
## obtain it from (mostly) uniprot or genbank. That is a bit of a pain, but
## we can write a couple of perl scripts that we can use to process the resulting tables
## to harvest species from there. Hopefully these species will be other species within
## Ensembl; otherwise we have to do something with a taxonomic tree, to find some reasonable
## information.

## export the gene as well as the transcript information
get.transcript.support.feature <- function(db, cred=db.cred){
    ## limit oursevles to canonical transcripts as that is what we have been looking at
    query1 <- paste("select a.stable_id as gene, b.stable_id as transcript, d.hit_name, d.evalue, d.perc_ident, e.db_name",
                   "from gene a inner join transcript b on a.canonical_transcript_id=b.transcript_id",
                   "inner join transcript_supporting_feature c on b.transcript_id=c.transcript_id and c.feature_type='dna_align_feature'",
                   "inner join dna_align_feature d on c.feature_id=d.dna_align_feature_id",
                   "inner join external_db e on d.external_db_id=e.external_db_id;")
    query2 <- gsub("dna_align", "protein_align", query1);
    db <- ensembl.connect(cred, db)
    tbl1 <- dbGetQuery(db, query1)
    tbl2 <- dbGetQuery(db, query2)
    dbDisconnect(db)
    list("dna"=tbl1, "protein"=tbl2)
}

transcript.feat <- lapply(core.dbs, get.transcript.support.feature, cred=db.cred)

## we can merge support derived from nucleic acids (which we expect will come from the same
## species (though this will need to be checked) with support from RNA-seq data.

all( names(transcript.cred) == names(transcript.feat) )
## this is TRUE

transcript.nuc <- mapply( function(x, y){
    union( x$stable_id, y$dna$transcript )},
    transcript.cred, transcript.feat )

## then we can make a plot as before, but here using the union of RNA-seq and
## dna_align_features
##

orth.tr.nuc <- t(sapply(names(sp.y), function(sp){
    tr <- unique(orth$tr[ !is.na( orth$tr[,sp] ), sp])
    m <- c('n'=length(tr),
           'seq'=sum( tr %in% transcript.nuc[[sp]] ))
    m <- c(m, 'r'=m[2] / m[1])
}))
rownames(orth.tr.nuc) <- names(sp.y)


## we can make a plot of the total number of transcripts with supporting information
## for these.
### 
par(mar=c(5.1, 8.1, 4.1, 2.1))
sp <- names(sp.y)
y <- barplot( orth.tr.nuc[,'r.seq'], las=2, names.arg="", col=sp.col[sp], xaxs='i', horiz=T )
mtext( sp, side=2, at=y, cex=0.5, las=2)

## export the data from transcript.nuc so that we can assign taxon identities to all
## identifiers
invisible(sapply(names(transcript.feat), function(x){
    write.table( transcript.feat[[x]]$dna, paste( "dna_align_features/", x, "_dna_align_feature.txt", sep=""), quote=FALSE, sep="\t" )
    write.table( transcript.feat[[x]]$protein, paste("protein_align_features/", x, "_protein_align_feature.txt", sep=""), quote=FALSE, sep="\t" )
}))


## having now obtained the taxons for the sequences supporting these features we can
## reload the data sets:

dna.sup.tax.files <- list.files("dna_align_features", pattern="tax.txt", full.names=TRUE)
prot.sup.tax.files <- list.files("protein_align_features", pattern="tax.txt", full.names=TRUE)

## read in:

dna.sup.tax <- lapply(dna.sup.tax.files, read.table, sep="\t", stringsAsFactors=FALSE, quote="", header=TRUE )
prot.sup.tax <- lapply(prot.sup.tax.files, read.table, sep="\t", stringsAsFactors=FALSE, quote="" )

names(dna.sup.tax) <- sub("[^/]+/([^.]+\\.[^_]+).+", "\\1", dna.sup.tax.files)
names(prot.sup.tax) <- sub("[^/]+/([^.]+\\.[^_]+).+", "\\1", prot.sup.tax.files)

barplot(sapply(dna.sup.tax, nrow), horiz=T, las=2)
barplot(sapply(prot.sup.tax, nrow), horiz=T, las=2)

## far more support from protein than from dna
## how about species in each one.

lapply( dna.sup.tax, function(x){ table(x$taxon) })
## that is generally only sequences from the species itself. But in any case for most species
## that just isn't that much data.

## Let us also create a new table similar to orth.tr.cred
## but which also includes this information:
sum.dna.cred <- function(sp){
    tr <- unique(orth$tr[ !is.na( orth$tr[,sp] ), sp])
    genus.species <- strsplit(sp, ".", fixed=TRUE)[[1]]
    dna.sup.b <- with(dna.sup.tax[[sp]],
                      grepl(genus.species[1], taxon, ignore.case=TRUE) &
                      grepl(genus.species[2], taxon, ignore.case=TRUE))
    rna.seq.b <- tr %in% transcript.cred[[sp]][,'stable_id']
    dna.all <- tr %in% dna.sup.tax[[sp]]$transcript
    dna.in <- tr %in% dna.sup.tax[[sp]]$transcript[ dna.sup.b ]
    dna.out <- tr %in% dna.sup.tax[[sp]]$transcript[ !dna.sup.b ]
    m <- c('n'=length(tr), 'rna.seq'=sum(rna.seq.b), 'dna.all'=sum(dna.all),
           'dna.in'=sum(dna.in), 'dna.out'=sum(dna.out), 'comb'=sum(rna.seq.b | dna.in))
}

orth.tr.cred.2 <- t(sapply(names(sp.y), sum.dna.cred))
rownames(orth.tr.cred.2) <- names(sp.y)

cairo_pdf("rna_supporting_evidence.pdf", width=10, height=20, onefile=TRUE)
par(mar=c(5.1, 8.1, 4.1, 2.1))
sp <- names(sp.y)
y <- barplot( orth.tr.cred.2[,'rna.seq'], las=2, names.arg="", col=sp.col[sp], xaxs='i', horiz=T, main='RNA seq' )
mtext( sp, side=2, at=y, cex=0.5, las=2)
##
sp <- names(sp.y)
y <- barplot( orth.tr.cred.2[,'dna.in'], las=2, names.arg="", col=sp.col[sp], xaxs='i', horiz=T, main='RNA within' )
mtext( sp, side=2, at=y, cex=0.5, las=2)
##
sp <- names(sp.y)
y <- barplot( orth.tr.cred.2[,'dna.all'], las=2, names.arg="", col=sp.col[sp], xaxs='i', horiz=T, main='RNA all' )
mtext( sp, side=2, at=y, cex=0.5, las=2)
##
sp <- names(sp.y)
y <- barplot( orth.tr.cred.2[,'comb'], las=2, names.arg="", col=sp.col[sp], xaxs='i', horiz=T, main='RNA seq or RNA within' )
mtext( sp, side=2, at=y, cex=0.5, las=2)
##
sp <- names(sp.y)
y <- barplot( orth.tr.cred.2[,'comb'] / orth.tr.cred.2[,'n'], las=2, names.arg="", col=sp.col[sp], xaxs='i', horiz=T,
             main='RNA seq or RNA within, proportions' )
mtext( sp, side=2, at=y, cex=0.5, las=2)
##
dev.off()

## note that we redefine prot.sup.tax.tbl further down after we have filled in species for
## ens databases
prot.sup.tax.tbl <- lapply( prot.sup.tax, function(x){ sort(table(x$taxon), decreasing=TRUE) })
sapply( prot.sup.tax.tbl, head )

prot.sup.db.tbl <- lapply( prot.sup.tax, function(x){ sort(table(x$db_name), decreasing=TRUE) })
sapply( prot.sup.db.tbl, head )

## We need to handle non-uniprot databases. Thesea are usually, "Ens_Xy_" where
## X_y indicate the species

sort( table(unlist( lapply( prot.sup.tax, function(x){ grep("Uniprot", x$db_name, invert=TRUE, value=TRUE) }) )))
##                Ens_Gg_gene                Ens_Hs_gene 
##                          5                         11 
##                        PRF                        PDB 
##                         19                         34 
##               IMGT/GENE_DB               IMGT/LIGM_DB 
##                         97                        218 
##         Ens_Mm_translation         Ens_Dr_translation 
##                        297                        323 
## Genoscope_pred_translation          Ens_Gg_transcript 
##                        491                        593 
##         Ens_Ac_translation         Ens_Bt_translation 
##                        682                       1118 
##         Ens_Tg_translation         Ens_Cf_translation 
##                       1483                       1971 
##         Ens_Gg_translation                 protein_id 
##                       6406                       6414 
##             RefSeq_peptide         Ens_Ga_translation 
##                       7127                      18179 
##          Ens_Dr_transcript          Ens_Mm_transcript 
##                      28173                     133947 
##         Ens_Hs_translation          Ens_Hs_transcript 
##                     204459                     308058 

## We can map the Ens db names to species
## from danio_rerio_core_98_11
##
ext.db.query <- 'select db_name, db_display_name from external_db where db_name regexp "^Ens"'
db <- ensembl.connect(db.cred, "danio_rerio_core_98_11")
ext.ens.db <- dbGetQuery( db, ext.db.query )
dbDisconnect(db)

## we need to manually insert the species names here:
ext.ens.db

##                     db_name                                    db_display_name
## 1               Ens_Hs_gene                                 Ensembl Human Gene
## 2         Ens_Hs_transcript                           Ensembl Human Transcript
## 3        Ens_Hs_translation                          Ensembl Human Translation
## 4                      ENSG                                       Ensembl gene
## 5                ENST_ident  Ensembl transcript having exact match with Havana
## 6                  ENST_CDS         Ensembl transcript sharing CDS with Havana
## 7         Ens_Mm_transcript                           Ensembl Mouse Transcript
## 8        Ens_Mm_translation                          Ensembl Mouse Translation
## 9        Ens_Cf_translation                            Ensembl Dog Translation
## 10        Ens_Dr_transcript                       Ensembl Zebrafish Transcript
## 11       Ens_Dr_translation                      Ensembl Zebrafish Translation
## 12              Ens_Gg_gene                               Ensembl Chicken Gene
## 13        Ens_Gg_transcript                         Ensembl Chicken Transcript
## 14       Ens_Gg_translation                        Ensembl Chicken Translation
## 15        Ens_Tr_transcript                            Ensembl Fugu Transcript
## 16       Ens_Ga_translation                    Ensembl Stickleback Translation
## 17        Ens_Mg_transcript                          Ensembl Turkey Transcript
## 18        Ens_Tg_transcript                      Ensembl Zebrafinch Transcript
## 19       Ens_Fc_translation                            Ensembl Cat Translation
## 20             ENS_LRG_gene                        LRG display in Ensembl gene
## 21       ENS_LRG_transcript                  LRG display in Ensembl transcript
## 22              Ens_Ga_gene                           Ensembl Stickleback Gene
## 23        Ens_Ga_transcript                     Ensembl Stickleback Transcript
## 24       Ens_Ss_translation                            Ensembl Pig Translation
## 25       Ens_Tg_translation                     Ensembl Zebrafinch Translation
## 26              Ens_Lc_gene                            Ensembl Coelacanth Gene
## 27       Ens_Ac_translation                   Ensembl Anole Lizard Translation
## 28       Ens_Bt_translation                            Ensembl Cow Translation
## 29           Ensembl_Plants                                     Ensembl Plants
## 30       Ens_Mg_translation                         Ensembl Turkey Translation
## 31               ENSP_ident Ensembl translation having exact match with Havana
## 32          Ensembl_Metazoa                                    Ensembl Metazoa
## 33      Ens_Mmu_translation                        Ensembl Macaque Translation
## 34 ensembl_internal_synonym                           Ensembl internal synonym
## 35            Ensembl_Fungi                                      Ensembl Fungi
## 36         Ensembl_Protists                                   Ensembl Protists
## 37           Ensembl_Plants                                     Ensembl Plants

ext.ens.db.sp <- c("Homo sapiens", "Homo sapiens", "Homo sapiens", NA, NA, NA,
                   rep("Mus musculus", 2), "Canis familiaris", rep("Danio rerio", 2),
                   rep("Gallus gallus", 3), "Takifugu rubripes", "Gasterosteus aculeatus",
                   "Meleagris gallopavo", "Taeniopygia guttata", "Felis catus", NA,
                   NA, rep("Gasterosteus aculeatus", 2) , "Sus scrofa", "Taeniopygia guttata",
                   "Latimeria chalumnae", "Anolis carolinensis", "Bos taurus", NA, "Meleagris gallopavo",
                   NA, NA, "Macaca mulatta", NA, NA, NA, NA)

ext.ens.db.sp[ ext.ens.db.sp %in% sp.names ]
ext.ens.db.sp[ ! ext.ens.db.sp %in% sp.names ] ## NAs

sum( is.na(ext.ens.db.sp) )  ## 12
sum( !is.na(ext.ens.db.sp) )  ## 25
sum( !is.na(ext.ens.db.sp) && ! ext.ens.db.sp %in% sp.names  )  ## 0

## so we can use the sp.names as a translation:
names(ext.ens.db.sp) <- ext.ens.db[,1]

for(i in 1:length( prot.sup.tax )){
    if(nrow(prot.sup.tax[[i]]) > 1){
        b1 <- is.na(prot.sup.tax[[i]][,'taxon']) 
        b2 <- prot.sup.tax[[i]][,'db_name'] %in% names(ext.ens.db.sp)
        print(paste(names(prot.sup.tax)[i], sum(b1), sum(b2), sum(b1 & b2)))
        prot.sup.tax[[i]][b1 & b2,'taxon'] <- ext.ens.db.sp[ prot.sup.tax[[i]][b1 & b2, 'db_name'] ]
    }
}

## and now we should remake the prot.sup.tax.tbl
prot.sup.tax.tbl <- lapply( prot.sup.tax, function(x){ sort(table(x$taxon), decreasing=TRUE) })
sapply(prot.sup.tax.tbl, head)

## we can make a reasonable table out of that:
prot.sup.tax.tbl.head <- t(sapply( prot.sup.tax.tbl, function(x){
    n <- 5
    if(length(x) == 0)
        return(rep("NA", n))
    e <- ifelse( length(x) >= n, n, length(x) );
    str <- paste(names(x)[1:e], x[1:e], sep="\n")
    if(e < n)
        str <- c(str, rep("NA", n-e))
    str
}))

## from
which(! unique( sapply( strsplit( as.character(prot.sup.tax.tbl.head), "\n"), function(x){ x[1] }) ) %in% sp.names)
## [1]  8 10 20 22 25 29 30 32 34 36 39 41

unique( sapply( strsplit( as.character(prot.sup.tax.tbl.head), "\n"), function(x){ x[1] }) )
##  [1] "Danio rerio"                "Homo sapiens"               "Gallus gallus"              "Bos taurus"                
##  [5] "Callorhinchus milii"        "Mus musculus"               "Gasterosteus aculeatus"     "NA"                        
##  [9] "Xenopus tropicalis"         "Salmo salar"                "Callithrix jacchus"         "Ictalurus punctatus"       
## [13] "Macaca mulatta"             "Neovison vison"             "Fukomys damarensis"         "Pan troglodytes"           
## [17] "Rattus norvegicus"          "Macaca fascicularis"        "Pongo abelii"               "Crotalus adamanteus"       
## [21] "Sus scrofa"                 "Xenopus laevis"             "Canis familiaris"           "Cavia porcellus"           
## [25] "Mustela putorius furo"      "Taeniopygia guttata"        "Ictidomys tridecemlineatus" "Oryctolagus cuniculus"     
## [29] "Lithobates catesbeianus"    "Canis lupus familiaris"     "Capra hircus"               "Heterocephalus glaber"     
## [33] "Mesocricetus auratus"       "Cricetulus griseus"         "Monodelphis domestica"      "Micrurus fulvius"          
## [37] "Oreochromis niloticus"      "Ailuropoda melanoleuca"     "Neotoma lepida"             "Felis catus"               
## [41] "Crotalus horridus"         

## we want to extend, names.sp and sp.col so that we can set domain colours for the table

names.sp <- c(names.sp, 'NA'='NA', 'Salmo salar'='salmo.salar', 'Crotalus adamanteus'='crotalus.adamanteus',
  'Xenopus laevis'='xenopus.laevis', 'Mustela putorius furo'='mustela.putorius',
  'Lithobates catesbeianus'='lithobates.catesbeianus', 'Canis lupus familiaris'='canis.familiaris',
  'Heterocephalus glaber'='heterocephalus.glaber', 'Cricetulus griseus'='cricetulus.griseus',
  'Micrurus fulvius'='micrurus.fulvius', 'Neotoma lepida'='neotoma.lepida', 'Crotalus horridus'='crotalus.horridus')

which(! unique( sapply( strsplit( as.character(prot.sup.tax.tbl.head), "\n"), function(x){ x[1] }) ) %in% names(names.sp) )
## 0.
## then extend the colours..

tmp.col <- c('NA'=class.col['others'], salmo.salar=class.col['teleostei'], crotalus.adamanteus=class.col['sauria'],
             xenopus.laevis=class.col['others'], mustela.putorius=class.col['eutheria'], lithobates.catesbeianus=class.col['others'],
             heterocephalus.glaber=class.col['eutheria'], cricetulus.griseus=class.col['eutheria'],
             micrurus.fulvius=class.col['sauria'], neotoma.lepida=class.col['eutheria'], crotalus.horridus=class.col['sauria'])

## because R concatenates any existing names of that which is assigned
names(tmp.col) <- sub("(.+)(\\..+$)", "\\1", names(tmp.col))

sp.col <- c(sp.col, tmp.col)

prot.sup.table <- function(xlim=c(0,100), ylim=c(0,100), x=5, y=95, c.widths=NULL, cex=NA, row.margin=1, ...){
    plot.new()
    plot.window(xlim=xlim, ylim=ylim, xaxs='i', yaxs='i')
    sp <- rev(names(sp.y))
    prot.sup.tax.tbl.head[ prot.sup.tax.tbl.head == "NA" ] <- "\n"
    tbl <- cbind( sp.names[ sp ], prot.sup.tax.tbl.head[sp,] )
    tbl.sp <- names.sp[ sapply( strsplit( tbl, "\n" ), function(x){ x[1] }) ]
    txt.col <- matrix( sp.col[tbl.sp], nrow=nrow(tbl), ncol=ncol(tbl))
    ## find a reasonable cex..
    if(is.na(cex))
        cex <- 1
    pos <- plotTable(x,y, tbl,
                     c.widths=c.widths, text.adj=c(0,1), text.col=txt.col,
                     txt.col.byrow=FALSE, cex=cex, row.margin=row.margin,
                     doPlot=FALSE, ... )
    ## take the median of the differences
    cell.h <- median( -(diff(pos$t)) )
    cell.r <- cell.h / strheight("A", cex=cex)
##
    h <- strheight("A", cex=cex) * cell.r * nrow(tbl)
    while( h > diff(ylim) ){
            cex <- cex * (diff(ylim) / (1.01 * h))
            h <- strheight("A", cex=cex) * cell.r * nrow(tbl)
    }
    pos <- plotTable(x,y, tbl,
                     c.widths=c.widths, text.adj=c(0,1), text.col=txt.col,
                     txt.col.byrow=FALSE, cex=cex, row.margin=row.margin, ... )
    invisible(c(pos, cex=cex))
}

prot.sup.table()
## consider drawing a tree on the left of this table to make the point more clear.
require(ape)

ex.align.2.k2.nj <- readRDS( "../R_172_genomes/ex_align_2_k2_nj.rds" )
## for consistency
ex.align.2.k2.nj$tip.label <- sub("_", ".", ex.align.2.k2.nj$tip.label)

## we can also read in the inferred sizes. Remember these are discretised to
## log2(size) / 10

int.s.inf.2 <- readRDS("../R_trees_distances_2/int_s_inf_2.rds")

tree.lines <- readRDS("../R_trees_distances/tree_lines.rds")
clade.nodes <- readRDS("../R_trees_distances/clade_nodes.rds")

## this depends on internal ordering of nodes, but seems to work ok.
nodes.col <- rep('black', length(clade.nodes$mammalia))

for(clade in names(clade.nodes))
    nodes.col[ clade.nodes[[clade]] ] <- class.col[clade]

plot.tree.table <- function(lwd=1, tbl.x=50, tbl.cex=NA, row.margin=1){
    row.bg=c(rgb(0.9, 0.9, 0.9), rgb(0.8, 0.8, 0.8))
    tbl.pos <- prot.sup.table(x=tbl.x, y=100, cex=tbl.cex, row.bg=row.bg, row.margin=row.margin)
##    tbl.pos <- prot.sup.table(x=tbl.x, y=100, cex=tbl.cex, row.bg=row.bg, c.widths=rep((100-tbl.x)/6, 6), row.margin=row.margin)
    t.y <- with(tbl.pos, (t + b)/2)
##    abline(h=t.y)
##    plot.new()
##    with(tree.lines, plot.window(xlim=c(0, 1.5*max(x[,2])), ylim=c(0.25, 1.01 * max(y)), xaxs='i', yaxs='i'))
    xw <- tbl.x
    xr <- with(tree.lines, c(0, max(x[,1:2])))
    tree.lines$x[,1:2] <- xw * (tree.lines$x[,1:2] - xr[1]) / diff(xr)
    tree.lines$v[,'x'] <- xw * (tree.lines$v[,'x'] - xr[1]) / diff(xr)
    tree.lines$y <- with(tree.lines, ((y - min(y)) / diff(range(y))) * diff(range(t.y)) + min(t.y))
    tree.lines$v[,2:3] <- with(tree.lines, ((v[,2:3]-min(v[,2:3])) / diff(range(v[,2:3]))) * diff(range(t.y)) + min(t.y))
    with(tree.lines, segments(x[,1], y, x[,2], y, col=nodes.col[x[,'p']], lwd=lwd) )
    with(tree.lines, segments(v[,'x'], v[,'y1'], v[,'x'], v[,'y2'], col=nodes.col[v[,'node']], lwd=lwd))
    with(tree.lines, {
        b <- !(nodes %in% ex.align.2.k2.nj$edge[,1])
        segments( x[b,2], y[b], tbl.x, y[b], col=nodes.col[ x[b,'p'] ], lty=2)
        ## offset <- strwidth("9999", cex=0.6) * 1.25
        ## text( x[b,2] + offset, y[b], sp.names[ex.align.2.k2.nj$tip.label[ nodes[b] ]],
        ##      cex=tbl.cex, adj=c(0,0.5) )
    })
    legend('bottomleft', legend=names(class.col), text.col=class.col, box.lty=0,
           inset=c(0,0.05))
    invisible(tbl.pos)
}

prot.sup.table()

pdf("protein_supporting_evidence_source.pdf", width=10, height=37)
pos <- plot.tree.table(lwd=1, tbl.x=35, row.margin=0.25, tbl.cex=1)
dev.off()


### We also have a complaint about not having enough teleost coverage..
### the only way to consider this problem is to look at the total known
### tree. For this we have to read in the tree data from NCBI..

node.lines <- readLines("~/genomes/ncbi_taxonomy_2021_09/nodes.dmp")
node.lines <- strsplit(node.lines, split="\t|\t", fixed=TRUE)

table( sapply(node.lines, length) )
##      13 
## 2359831

nodes.matrix <- t(sapply(node.lines, eval))
## the comments field is a bit screwed up. 

## we only need these to build the tree. Later we can extract
## subsets and when doing so we can calculate suitable x and y
## positions
tax.nodes <- data.frame('id'=as.numeric(nodes.matrix[,1]),
                        'p.id'=as.numeric(nodes.matrix[,2]),
                        'rank'=nodes.matrix[,3], stringsAsFactors=FALSE)

## and then we probably need to do the same here.
names.lines <- readLines("~/genomes/ncbi_taxonomy_2021_09/names.dmp")
names.lines <- strsplit(names.lines, "\t\\|\t?" )
table( sapply(names.lines, length))  ## all 4 which is correct

names.m <- t(sapply(names.lines, eval))
table( names.m[,4] )
##         acronym           authority          blast name         common name 
##            1587              579655                 228               14556 
## equivalent name     genbank acronym genbank common name     genbank synonym 
##           52130                 485               29987                1095 
##         in-part            includes     scientific name             synonym 
##             733               61926             2359831              199762 
##   type material 
##          168237 

## we are only really interested in the scientific name..
b <- names.m[,4] == "scientific name"

tax.names <- data.frame('id'=names.m[,1], 'name'=names.m[,2], 'uname'=names.m[,3], stringsAsFactors=FALSE )
tax.sciname <- tax.names[b,]

length(unique(tax.sciname$id)) == nrow(tax.sciname)
## TRUE
rownames(tax.sciname) <- tax.sciname$id

nrow(tax.nodes) - sum( tax.nodes$id %in% tax.sciname$id )
## 9 are missing.. still we should be able to do:

tax.nodes <- cbind( tax.nodes, 'name'=tax.sciname[ as.character(tax.nodes$id), 'name' ], stringsAsFactors=FALSE )

## let us make a tree from a parent..
## unfortunately has to use cbind to extend
## the tree...
## this relies on a global global.y
## and is incredibly slow.
extend.tree <- function(root, p.id, x, nodes, filter='unclassified|environmental', min.rank='genus'){
    i <- which(nodes$id == root)
    i.c <- which(nodes$p.id == root)
    i.c <- i.c[ !grepl(filter, nodes$name[i.c]) ]
    x <- x + 1
    ## print(paste(root, x))
    ## print(nodes[i,])
    if(length(i.c) == 0 || nodes[i,'rank'] == min.rank){
        y <- global.y
        global.y <<- global.y + 1
        return( c( 'node'=root, 'x'=x, 'y'=y, 'i'=i, 'p'=p.id ))
    }
    children.l <- lapply( nodes[i.c, 'id'], extend.tree, x=x, nodes=nodes, p.id=root, filter=filter, min.rank=min.rank )
    children <- do.call(rbind, children.l)
    y <- mean( children[,'y'] )
    global.y <<- global.y + 1
    rbind( children, c(root, x, y, i, p.id) )
}

make.tree <- function(root, nodes, x=0, p.id=-1, filter='unclassified|environmental', min.rank='genus'){
    tree <- extend.tree( root, p.id, x, nodes, filter=filter, min.rank=min.rank )
    tree <- cbind( tree, 'p.i'=match(tree[,'p'], tree[,'node']) )
    cbind( as.data.frame(tree), nodes[ tree[,'i'], c('name', 'rank')], stringsAsFactors=FALSE )
}

make.tree.segments <- function(tree){
    x0 <- tree$x
    x1 <- tree$x[ tree$p.i ]
##    xm <- (x0 + x1) / 2
    y0 <- tree$y
    y1 <- tree$y[ tree$p.i ]
    df <- data.frame(x0=c(x0,x1), y0=c(y0,y0), x1=c(x1,x1), y1=c(y0,y1), name=c(tree$name, tree$name), i=rep(1:nrow(tree), 2))
    na.b <- with(df, is.na(x1) | is.na(y1))
    df[ !na.b, ]
}

## use bezier curves
make.fancy.tree.segments <- function(tree, b.n=10, ctl.p=1){
    x0 <- tree$x
    x1 <- tree$x[ tree$p.i ]
##    xm <- (x0 + x1) / 2
    xm <- x0 + (x1-x0) * ctl.p
    y0 <- tree$y
    y1 <- tree$y[ tree$p.i ]
    x <- rbind( x0, xm, x1 )
    y <- rbind( y0, y0, y1 )
    input.pts <- cbind(x=as.numeric(x), y=as.numeric(y))
    t.p <- seq(0,1, length.out=b.n)
    ur <- 2:b.n
    lr <- ur - 1
    seg.l <- lapply( seq(1,nrow(input.pts), 3), function(i){
        pts <- bezier.pts(input.pts[i:(i+2),], t=t.p)
        cbind(x0=pts[lr,'x'], y0=pts[lr,'y'], x1=pts[ur,'x'], y1=pts[ur,'y'], i=(1 + i %/% 3))
    })
    as.data.frame(do.call(rbind, seg.l))
}

ascend.tree <- function(i, tree=teleost.tree){
    if(is.na(tree[i,'p.i']))
        return( tree[i, ] )
    ancestors <- ascend.tree( tree[i, 'p.i'], tree )
    rbind( tree[i, ], ancestors )
}



tmp <- make.tree( 7925, tax.nodes, min.rank='species' )
tmp.seg <- make.tree.segments( tmp )
plot.new()
with(tmp.seg, plot.window(xlim=range(c(x0, x1)), ylim=range(c(y0, y1)) ))
do.call( segments, tmp.seg )


## 32443 is Teleostei
global.y <- 0
teleost.tree <- make.tree( 32443, tax.nodes, min.rank='genus' )
teleost.tree.seg <- make.tree.segments( teleost.tree )

teleost.tree.fseg <- make.fancy.tree.segments( teleost.tree, ctl.p=1, b.n=30 )

plot.new()
with(teleost.tree.seg, plot.window(xlim=range(c(x0, x1)), ylim=range(c(y0, y1)) ))
do.call( segments, c(teleost.tree.fseg[,1:4], col='grey') )

global.y <- 0
vertebrate.tree <- make.tree( 7742, tax.nodes, min.rank='genus' )
vertebrate.tree.seg <- make.tree.segments( vertebrate.tree )
vertebrate.tree.fseg <- make.fancy.tree.segments( vertebrate.tree, ctl.p=1, b.n=30 )

## that sort of works. We can modify the tree.segments to make prettier lines,
## but that isn't that important. For now though we want to see if we can map the species
## to this tree.

sp.genus <- sapply( strsplit(sp.names, " ", fixed=TRUE), function(x){ x[1] })

b <- sp.genus %in% teleost.tree$name
sum(b) ## only 54...
sum( sp.class[,'teleostei'] ) ## 54

genus.b <- teleost.tree$name %in% sp.genus

## lets check the vertebrate tree
sum( sp.genus %in% vertebrate.tree$name ) ## 172 good.

## so we actually do have all of these. We can actually draw the lines
## associated with them, because of the way that the teleost.tree.seg are
## arranged:

b <- teleost.tree$name %in% sp.genus
## only 48, some duplication due to duplication from multiple
## species.
b <- teleost.tree$name[ teleost.tree.fseg$i ] %in% sp.genus ## though repeating might work too
with( teleost.tree.fseg, segments( x0[b], y0[b], x1[b], y1[b], col='red', lwd=2 ))
## one major and minore branch missing:

with(teleost.tree, identify(x, y, name))

taxons <- unique( do.call(rbind, lapply( b.i, ascend.tree )))
taxons.b <- teleost.tree$name[teleost.tree.fseg$i] %in% taxons$name ## 368
with( teleost.tree.fseg, segments( x0[b], y0[b], x1[b], y1[b], col='blue', lwd=2 ))

id.i <- with(teleost.tree, identify(x, y, name, pos=TRUE))

## lets add percomorphacae to the labels:
id.i$ind <- c(id.i$ind, 5133, 1804)
id.i$pos <- c(id.i$pos, 3, 1)
          
## lets make a pdf that's big enough that we can maybe see some
## details:
cairo_pdf("taxon_representation.pdf", width=10, height=30, pointsize=8)
plot.new()
with(teleost.tree.seg, plot.window(xlim=c(min(x1), max(x0)*1.1), ylim=range(c(y0, y1)) ))
do.call( segments, c(teleost.tree.fseg[,1:4], col='black', lwd=0.1) )
with( teleost.tree.fseg[taxons.b,], segments( x0, y0, x1, y1, col='blue', lwd=2 ))
with(teleost.tree, text( x[id.i$ind], y[id.i$ind], name[id.i$ind], pos=id.i$pos) )
glab.y <- teleost.tree$y[genus.b]
glab <- teleost.tree$name[genus.b]
glab.y2 <- space.text.y( glab.y, glab )
##text( teleost.tree$x[genus.b], glab.y2, glab, pos=4 )
with(teleost.tree, text( max(x[genus.b]), glab.y2, glab, pos=4 ))
with(teleost.tree, segments( x[genus.b], glab.y, max(x[genus.b]), glab.y2, lty=3 ))
##with(teleost.tree, text( x[genus.b], y[genus.b], name[genus.b], pos=4))
dev.off()

vert.taxons <- unique( do.call( rbind, lapply( which( vertebrate.tree$name %in% sp.genus ), ascend.tree, tree=vertebrate.tree )))
vert.b <- vertebrate.tree$name[ vertebrate.tree.fseg$i ] %in% vert.taxons$name
vert.gb <- vertebrate.tree$name %in% sp.genus
    
cairo_pdf("vertebrate_representation.pdf", width=10, height=30, pointsize=8)
plot.new()
with(vertebrate.tree.seg, plot.window(xlim=c(min(x1), max(x0)*1.1), ylim=range(c(y0, y1)) ))
do.call( segments, c(vertebrate.tree.fseg[,1:4], col='black', lwd=0.1) )
with( vertebrate.tree.fseg[vert.b,], segments( x0, y0, x1, y1, col='blue', lwd=2 ))
##with(vertebrate.tree, text( x[id.i$ind], y[id.i$ind], name[id.i$ind], pos=id.i$pos) )
glab.y <- vertebrate.tree$y[vert.gb]
glab <- vertebrate.tree$name[vert.gb]
glab.y2 <- space.text.y( glab.y, glab )
##text( teleost.tree$x[genus.b], glab.y2, glab, pos=4 )
with(vertebrate.tree, text( max(x[vert.gb]), glab.y2, glab, pos=4 ))
with(vertebrate.tree, segments( x[vert.gb], glab.y, max(x[vert.gb]), glab.y2, lty=3 ))
dev.off()

## how many teloest species are percomorphae..
sp.bi <- which( sp.genus %in% teleost.tree$name ) ## length of 54 correct
names(sp.bi) <- names(sp.genus[sp.bi])

tree.bi <- match( sp.genus[sp.bi], teleost.tree$name)
names(tree.bi) <- sp.names[ sp.bi ]  ## this gives species names

teleost.lineages <- lapply( tree.bi, ascend.tree )
## for a nice summary:
teleost.lineages.n <- sort(table(do.call(rbind, teleost.lineages)$name))

teleost.lineages.perc <- as.logical(sapply( teleost.lineages, function(x){ sum(x$name == "Percomorphaceae") }))
teleost.lineages.clupeo <- as.logical(sapply( teleost.lineages, function(x){ sum(x$name == "Clupeocephala") }))
teleost.lineages.otomorpha <- as.logical(sapply( teleost.lineages, function(x){ sum(x$name == "Otomorpha") }))
sum(teleost.lineages.perc)  ## 42
sum(teleost.lineages.clupeo) ## 52

i <- match( sapply( strsplit(names(teleost.lineages), " ", fixed=TRUE), function(x){ x[1] }), teleost.tree$name )
o <- order( teleost.tree$y[i] )
teleost.lineages <- teleost.lineages[ o ]

perc.sp <- sapply( names(teleost.lineages)[ teleost.lineages.perc ], function(x){ tolower(sub(" ", "_", x)) })
n.perc.sp <- sapply( names(teleost.lineages)[ !teleost.lineages.perc ], function(x){ tolower(sub(" ", "_", x)) })


## we should consider looking at the distributions of the subsets of these
## species. We can for the sake of completeness consider both the full
## set of introns and the orthologous ones:

## just read in the data again.
orth.stats <- read.table("../family_members/orthologue_transcripts/exon_intron_stats.csv", sep="\t", stringsAsFactors=FALSE)
colnames(orth.stats) <- c('sp', 'db', 'family', 'gene', 'transcript', 'chr', 'strand', 'pos', 'ex.s', 'in.s')

gene.stats <- read.table("../R_basic_stats/all_genes_exon_intron_stats.csv", sep="\t", stringsAsFactors=FALSE)
colnames(gene.stats) <- c('sp', 'db', 'gene', 'transcript', 'chr', 'strand', 'pos', 'ex.s', 'in.s')

## orth.stats has about 24% of the rows of gene.stats; as expected.
## get intron sizes
orth.int.s <- lapply(strsplit(orth.stats$in.s, ',', fixed=TRUE ), as.numeric)
gene.int.s <- lapply(strsplit(gene.stats$in.s, ',', fixed=TRUE ), as.numeric)

## we can now obtain histograms for all of the data.
hist.all <- hist( log2( unlist( gene.int.s)), breaks=75 )

gene.int.s.h <- tapply( gene.int.s, gene.stats$sp, function(x){ hist(log2(unlist(x)), breaks=hist.all$breaks, plot=FALSE ) })
orth.int.s.h <- tapply( orth.int.s, orth.stats$sp, function(x){ hist(log2(unlist(x)), breaks=hist.all$breaks, plot=FALSE ) })

sp.sel <- c('danio_rerio','paramormyrops_kingsleyae', 'scleropages_formosus', 'parambassis_ranga',
            'gadus_morhua', 'esox_lucius', 'hucho_hucho', 'clupea_harengus', 'denticeps_clupeoides',
            'ictalurus_punctatus', 'electrophorus_electricus')

sp.sel <- n.perc.sp
sp.sel <- perc.sp[1:16]
sp.sel <- perc.sp[17:32]

par(mfrow=c(4,4))
invisible( lapply(sp.sel, function(x){
    plot(gene.int.s.h[[x]], main=sub("_", " ", x))
    abline(v=8, col='red')
}))


### let us get the histograms for the Ensembl 104 extended data set:
gene.stats.104 <- read.table("ensembl_104_exon_intron_stats.csv", sep="\t", stringsAsFactors=FALSE)
colnames(gene.stats.104) <- c('sp', 'db', 'gene', 'transcript', 'chr', 'strand', 'pos', 'ex.s', 'in.s')

gene.int.104.s <- lapply(strsplit(gene.stats.104$in.s, ",", fixed=TRUE ), as.numeric )

gene.int.104.s.all.h <- hist( log2( unlist(gene.int.104.s)), breaks=300 )

gene.int.104.s.h <- tapply( gene.int.104.s, gene.stats.104$sp,
                           function(x){ hist(log2(unlist(x)), breaks=gene.int.104.s.all.h$breaks, plot=FALSE) })

par(mfrow=c(4,4))
invisible( lapply(perc.sp[17:32], function(x){
    plot(gene.int.104.s.h[[x]], main=sub("_", " ", x))
    abline(v=8, col='red')
}))

par(mfrow=c(4,4))
invisible( lapply(n.perc.sp, function(x){
    plot(gene.int.104.s.h[[x]], main=sub("_", " ", x))
    abline(v=c(8,11), col='red')
}))

## to get the lineages of all the species...
## in the 104 data set
sp.104.names <- names( gene.int.104.s.h )
sp.104.names <- sub("_", " ", sp.104.names, fixed=TRUE )
substring(sp.104.names, 1, 1) <- toupper( substring(sp.104.names, 1, 1))
names(sp.104.names) <- sub("_", ".", names(gene.int.104.s.h), fixed=TRUE)
sp.104.genus <- sapply( strsplit(sp.104.names, " ", fixed=TRUE), function(x){ x[1] })

sp.bi <- which( sp.104.genus %in% teleost.tree$name ) ## 84 entries
names(sp.bi) <- names(sp.104.genus[sp.bi])

tree.bi <- match( sp.104.genus[sp.bi], teleost.tree$name)
names(tree.bi) <- sp.104.names[ sp.bi ]  ## this gives species names

teleost.104.lineages <- lapply( tree.bi, ascend.tree )
i <- match(sp.104.genus[sp.bi], teleost.tree$name)
o <- order( teleost.tree[i, 'y'] )

teleost.104.lineages <- teleost.104.lineages[ o ]

cairo_pdf("teleost_104_intron_sizes.pdf", width=18, height=24, onefile=TRUE)
par(mfrow=c(5,4))
for(nm in names( teleost.104.lineages )){
    plot( gene.int.104.s.h[[ sub(" ", "_", tolower(nm)) ]], main=nm )
    abline(v=c(8,11), col='red')
    with(par(),
         text(usr[2], usr[4], paste( teleost.104.lineages[[nm]][,'name'], collapse="\n"),
              adj=c(1,1))
         )
}
dev.off()

## let us calculate the proportion of introns that are below 75 bp.
lt75.104 <- tapply( gene.int.104.s, gene.stats.104$sp, function(x){ s = unlist(x); c( sum(s < 75), length(s) ) })
lt75.104 <- t(sapply(lt75.104, eval))

lt75.104 <- cbind(lt75.104, lt75.104[,1] / lt75.104[,2])
colnames( lt75.104 ) <- c('q', 'n', 'r')

lt75.104.col <- rgb( lt75.104[,'r'] / max(lt75.104[,'r']), 0, 0 )
names(lt75.104.col) <- rownames(lt75.104)

## have a closer look at:
sp.sel <- c("cyprinus_carpio", "sinocyclocheilus_anshuiensis",
            "electrophorus_electricus", "denticeps_clupeoides", "clupea_harengus",
            "salmo_salar", "esox_lucius", "hippocampus_comes", "paramormyrops_kingsleyae")

source("~/R/general_functions.R")
plot.sel <- function(sel, h=gene.int.104.s.h, col=hsvScale(1:length(sel))){
    dens <- sapply( h[sel], function(x){ x$density } )
    mids <- sapply( h[sel], function(x){ x$mids } )
    plot(mids[,1], dens[,1], xlim=range(mids), ylim=range(dens), xlab='log2 length', ylab='density', type='n')
    for(i in 1:ncol(dens))
        lines(mids[,i], dens[,i], col=col[i], lwd=3)
    legend('topright', sel, fill=col)
}

par(mfrow=c(2,1))
plot.sel( c(sp.sel[1:6]) )
plot.sel( c(sp.sel[7:10], 'danio_rerio') )

par(mfrow=c(2,2))
plot.sel( c('betta_splendens', 'danio_rerio', sp.sel[1:2]))
plot.sel( c('betta_splendens', 'danio_rerio', sp.sel[3:5]))
plot.sel( c('betta_splendens', 'danio_rerio', sp.sel[6:7]))
plot.sel( c('betta_splendens', 'danio_rerio', sp.sel[8:9], 'homo_sapiens'))


### let us do a PCA with the Ensembl 104 species. The problem I have is that I do
### not know which are teleosts and so on? But I can work that out from the
### the vertebrate.tree
vert.bi <- which(sp.104.genus %in% vertebrate.tree$name)
vert.104.bi <- match( sp.104.genus[vert.bi], vertebrate.tree$name )
vert.104.lineages <- lapply( vert.104.bi, ascend.tree, tree=vertebrate.tree )
names(vert.104.lineages) <- sub("\\.", "_", names(sp.104.names[vert.bi]))

is.lineage <- function(line, lineages=vert.104.lineages){
    as.logical( sapply( lineages, function(x){ sum( line %in% x$name ) }))
}

sp.104.class <- cbind( is.lineage('Teleostei'), is.lineage('Mammalia'), is.lineage('Eutheria'),
                      is.lineage('Sauria'), is.lineage('Aves'))
rownames(sp.104.class) <- names(vert.104.lineages)
colnames(sp.104.class) <- c('teleostei', 'mammalia', 'eutheria', 'sauria', 'aves')

sp.104.col <- rep('black', nrow(sp.104.class))
for( cn in colnames( sp.104.class )){
    sp.104.col[ sp.104.class[,cn] ] <- class.col[cn]
}
names(sp.104.col) <- rownames(sp.104.class)
                  

vert.dens <- sapply( gene.int.104.s.h[vert.bi], function(x){ x$density })
vert.dens.pca <- prcomp( t(vert.dens) )

plot(vert.dens.pca) ## reasonable amount of info in 2 first

with(vert.dens.pca, plot(x[,1], x[,2], col=sp.104.col, pch=19))
with(vert.dens.pca, identify(x[,1], x[,2], rownames(x)))

with(vert.dens.pca, plot(x[,1], x[,2], col=lt75.104.col[rownames(x)], pch=19))

with(vert.dens.pca, plot(x[,2], lt75.104[rownames(x), 'r'], pch=19))
with(vert.dens.pca, plot(x[,1], lt75.104[rownames(x), 'r'], pch=19))
## we should probably use as a cutoff 0.1, or 0.05..

vert.dens.pca.2 <- prcomp( t(vert.dens[ ,lt75.104[ colnames(vert.dens), 'r'] < 0.1 ]) )
plot(vert.dens.pca.2)

with(vert.dens.pca.2, plot(x[,1], x[,2], col=sp.104.col[rownames(x)], pch=19))
with(vert.dens.pca.2, identify(x[,1], x[,2], rownames(x)))

### I want to make a figure with some example teleost next to the taxonomic tree:
## we will want to make a selection on the basis of the position of the 104 teleost tree
## for this we can use the teleost tree and the 104 genus
## sp.104.genus

teleost.tree.104.bi <- which(teleost.tree$name %in% sp.104.genus) ## 63 (remember some genus get repated)
tel.104.taxons <- unique( do.call( rbind, lapply( teleost.tree.104.bi, ascend.tree, tree=teleost.tree )))
sum( sp.104.genus %in% teleost.tree$name ) ## 84

tel.104.genus.b <- teleost.tree$name %in% sp.104.genus
tel.104.fseg.b <- teleost.tree$name[ teleost.tree.fseg$i ] %in% tel.104.taxons$name

## and then to draw the tree with everything we can simply do:
with(teleost.tree.fseg, segments(x0, y0, x1, y1, col='grey', lwd=0.5 ))

plot.new()
with(teleost.tree.seg, plot.window(xlim=c(min(x1), max(x0)*1.1), ylim=range(c(y0, y1)) ))
do.call( segments, c(teleost.tree.fseg[,1:4], col='black', lwd=0.1) )
with( teleost.tree.fseg[tel.104.fseg.b,], segments( x0, y0, x1, y1, col='blue', lwd=2 ))
glab.y <- teleost.tree$y[tel.104.genus.b]
glab <- teleost.tree$name[tel.104.genus.b]
glab.y2 <- space.text.y( glab.y, glab )
with(teleost.tree, text( max(x[tel.104.genus.b]), glab.y2, glab, pos=4 ))
with(teleost.tree, segments( x[tel.104.genus.b], glab.y, max(x[tel.104.genus.b]), glab.y2, lty=3 ))

## let us make a selection for a set of species from different parts of the tree:

## make a list. for each member make a line plot with a legend.. 
genus.sel <- list(c('Tetraodon', 'Takifugu', 'Betta', 'Gasterosteus', 'Cynoglossus'),
                  c('Carassius', 'Cyprinus', 'Sinocyclocheilus'),
                  c('Danio', 'Astyanax', 'Ictalurus', 'Electrophorus'),
                  c('Denticeps', 'Clupea'),
                  c('Oncorhynchus', 'Salmo', 'Esox'),
                  c('Myripristis', 'Sphaeramia', 'Neogobius', 'Hippocampus'),
                  c('Monopterus', 'Mastacembelus', 'Betta' ),
                  c('Haplochromis', 'Maylandia', 'Pundamilia'),
                  c('Gasterosteus', 'Cyclopterus', 'Sander', 'Cottoperca'),
                  c('Takifugu', 'Tetraodon', 'Mola'),
                  c('Gadus', 'Scleropages', 'Paramormyrops'))



plot.sel.lines <- function(sel, data=vert.dens, x=gene.int.104.s.h[[1]]$mids, lwd=2, col=NA, ...){
    if( !all( sel %in% sp.104.genus ) )
        stop(paste("Unknown genus"), sel )
    b1 <- rep(FALSE, ncol(data))
    for(genus in sel){
        b1 <- b1 | grepl(genus, colnames(data), ignore.case=TRUE)
    }
    if(is.na(col))
        col <- rainbow( sum(b1), s=0.8, v=0.8, start=0.2 )
    data <- data[,b1]
    ## get labels
    lab <- sub("_", " ", colnames(data))
    substring(lab, 1, 1) <- toupper(substring(lab, 1, 1))
    plot(x, data[,1], ylim=range(data), col=col[1], type='l', lwd=lwd, ...)
    for(i in 2:ncol(data)){
        lines(x, data[,i], col=col[i], lwd=lwd)
    }
    legend( 'topright', legend=lab, fill=col )
}

cairo_pdf("tax_tree_104_sel.pdf", width=10, height=30, pointsize=12)
plot.new()
with(teleost.tree.seg, plot.window(xlim=c(min(x1), max(x0)*1.1), ylim=range(c(y0, y1)) ))
do.call( segments, c(teleost.tree.fseg[,1:4], col='black', lwd=0.1) )
with( teleost.tree.fseg[tel.104.fseg.b,], segments( x0, y0, x1, y1, col='blue', lwd=2 ))
glab.y <- teleost.tree$y[tel.104.genus.b]
glab <- teleost.tree$name[tel.104.genus.b]
glab.y2 <- space.text.y( glab.y, glab )
with(teleost.tree, text( max(x[tel.104.genus.b]), glab.y2, glab, pos=4, col=ifelse(glab %in% unlist(genus.sel), 'red', 'black'), cex=1 ))
with(teleost.tree, segments( x[tel.104.genus.b], glab.y, max(x[tel.104.genus.b]), glab.y2, lty=3 ))
dev.off()

cairo_pdf("selected_distributions_104.pdf", width=16, height=12)
par(mfrow=c(3,4))
for(sel in rev(genus.sel)){
    sel <- rev(sel)
    plot.sel.lines( sel, main=paste(sel, collapse=", "), xlab='log2 intron length', ylab='density', xlim=c(5,17) )
    abline(v=8, lty=2)
}
dev.off()

## I want to get the genome sizes:
require("RMySQL")
db.104 <- readLines("../core_db_104.txt")
genome.stats.104 <- lapply( db.104, function(x){ get.genome.stats( db.cred, x ) })

genome.sizes.104 <- sapply(genome.stats.104, function(x){ x[ x$statistic == 'ref_length', 'value' ] })
names(genome.sizes.104) <- sub("(.+)?_core_.+", "\\1", db.104)


### Plot pairwise data for version 98. Demonstrate that introns have decreased in
### size in Danio rerio as well as A mexicanus and P. natterei.

with(orth, plot(log2(l[,'lepisosteus.oculatus']), log2(l[,'danio.rerio']), cex=0.5,
                xlim=c(5.5, 18), ylim=c(5.5, 18), col=rgb(0,0,0,0.2)))

with(orth, plot(log2(l[,'lepisosteus.oculatus']), log2(l[,'astyanax.mexicanus']), cex=0.5,
                xlim=c(5.5, 18), ylim=c(5.5, 18), col=rgb(0,0,0,0.2)))

with(orth, plot(log2(l[,'lepisosteus.oculatus']), log2(l[,'pygocentrus.nattereri']), cex=0.5,
                xlim=c(5.5, 18), ylim=c(5.5, 18), col=rgb(0,0,0,0.2)))

with(orth, plot(log2(l[,'danio.rerio']), log2(l[,'astyanax.mexicanus']), cex=0.5,
                xlim=c(5.5, 18), ylim=c(5.5, 18), col=rgb(0,0,0,0.2)))

with(orth, plot(log2(l[,'danio.rerio']), log2(l[,'pygocentrus.nattereri']), cex=0.5,
                xlim=c(5.5, 18), ylim=c(5.5, 18), col=rgb(0,0,0,0.2)))

with(orth, plot(log2(l[,'astyanax.mexicanus']), log2(l[,'pygocentrus.nattereri']), cex=0.5,
                xlim=c(5.5, 18), ylim=c(5.5, 18), col=rgb(0,0,0,0.2)))

## to look at ancestral states..
## 77 is A. mexicanus, 76 is P. nattereri, 304 is their common ancestor
plot(int.s.inf.2[[304]]$state[,1], int.s.inf.2[[77]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)

plot(int.s.inf.2[[304]]$state[,1], int.s.inf.2[[76]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19, border=NA)

plot(int.s.inf.2[[77]]$state[,1], int.s.inf.2[[76]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19, border=NA)

## 78 is danio rerio, and 272 is the common ancestor between (A. mex, A. natt.) and d. rerio
plot(int.s.inf.2[[272]]$state[,1], int.s.inf.2[[78]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)

plot(int.s.inf.2[[272]]$state[,1], int.s.inf.2[[77]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)

plot(int.s.inf.2[[272]]$state[,1], int.s.inf.2[[76]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)

plot(int.s.inf.2[[78]]$state[,1], int.s.inf.2[[77]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)

plot(int.s.inf.2[[78]]$state[,1], int.s.inf.2[[76]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)

plot(int.s.inf.2[[77]]$state[,1], int.s.inf.2[[76]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)

## 215 is the eutherian common ancestor.
## But let us consider the common ancestor of mouse and humans..
## human is 155, mouse is 127
## common ancestor is: 202

plot(int.s.inf.2[[155]]$state[,1], int.s.inf.2[[127]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)

plot(int.s.inf.2[[202]]$state[,1], int.s.inf.2[[155]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)

plot(int.s.inf.2[[202]]$state[,1], int.s.inf.2[[127]]$state[,1], cex=0.5, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)

## Is there any commonality in which introns have grown in size?

int.growth.common <- function(anc, d1, d2, inf=int.s.inf.2, thr=80, min.s=55, invert=FALSE){
    b <- inf[[d1]]$state[,1] > min.s & inf[[d2]]$state[,1] > min.s & inf[[anc]]$state[,1] > min.s
    sm <- cbind( inf[[anc]]$state[b,1], inf[[d1]]$state[b,1], inf[[d2]]$state[b,1] )
    a.short <- sm[,1] < thr
    if(!invert){
        d1.grown <- a.short & sm[,2] > thr
        d2.grown <- a.short & sm[,3] > thr
    }else{
        d1.grown <- a.short & sm[,2] < thr
        d2.grown <- a.short & sm[,3] < thr
    }
    q <- sum(d1.grown & d2.grown)
    m <- sum(d1.grown)
    n <- sum(a.short) - sum(d1.grown)
    k <- sum(d2.grown)
    p <- phyper( q, m, n, k, lower.tail=FALSE )
    exp <- k * m / sum(a.short)
    c(q=q, m=m, n=n, k=k, p=p, exp=exp)
}

## this takes the output of simulate.is.evol and does the same thing as above..
int.growth.common.2 <- function(anc, d1, d2, tree, thr=80, min.s=55){
    sm <- sapply( c(anc, d1, d2), function(i){ tree[[i]]$size })
    sm <- sm[ apply(sm, 1, min) >= min.s, ]
    a.short <- sm[,1] < thr
    d1.grown <- a.short & sm[,2] > thr
    d2.grown <- a.short & sm[,3] > thr
    q <- sum(d1.grown & d2.grown)
    m <- sum(d1.grown)
    n <- sum(a.short) - sum(d1.grown)
    k <- sum(d2.grown)
    p <- phyper( q, m, n, k, lower.tail=FALSE )
    exp <- k * m / sum(a.short)
    c(q=q, m=m, n=n, k=k, p=p, exp=exp, r=q/exp)
}

length.retained <- function(anc, d1, d2, inf, thr=80, min.s=55){
    sm <- sapply( c(anc, d1, d2), function(i){ inf[[i]]$state[,1] })
    sm <- sm[ apply(sm, 1, min) >= min.s, ]
    a.long <- sm[,1] > thr
    d1.long <- a.long & sm[,2] > thr
    d2.long <- a.long & sm[,3] > thr
    q <- sum(d1.long & d2.long)
    m <- sum(d1.long)
    n <- sum(a.long) - sum(d1.long)
    k <- sum(d2.long)
    p <- phyper( q, m, n, k, lower.tail=FALSE )
    exp <- k * m / sum(a.long)
    c(q=q, m=m, n=n, k=k, p=p, exp=exp, r=q/exp)
}

length.retained.2 <- function(anc, d1, d2, tree, thr=80, min.s=55){
    sm <- sapply( c(anc, d1, d2), function(i){ tree[[i]]$size })
    sm <- sm[ apply(sm, 1, min) >= min.s, ]
    a.long <- sm[,1] > thr
    d1.long <- a.long & sm[,2] > thr
    d2.long <- a.long & sm[,3] > thr
    q <- sum(d1.long & d2.long)
    m <- sum(d1.long)
    n <- sum(a.long) - sum(d1.long)
    k <- sum(d2.long)
    p <- phyper( q, m, n, k, lower.tail=FALSE )
    exp <- k * m / sum(a.long)
    c(q=q, m=m, n=n, k=k, p=p, exp=exp, r=q/exp)
}
    

## D. rerio: 78
## A. mexicanus: 77
## P. nattereri: 76
## I. punctatus: 30
## E. electrophorus: 9
## their common ancestor: 272

int.growth.common( 272, 78, 77, inf=int.s.inf.2, thr=80, min.s=55 )
##        q         m         n         k         p       exp 
## 1377.000  3355.000 13083.000  7562.000     1.000  1543.406 

int.growth.common( 272, 78, 76, inf=int.s.inf.2, thr=80, min.s=55 )
##        q         m         n         k         p       exp 
## 1223.000  3470.000 13580.000  7138.000     1.000  1452.719

int.growth.common( 272, 78, 30, inf=int.s.inf.2, thr=80, min.s=55 )
##       q          m          n          k          p        exp 
## 46.0000  3439.0000 13353.0000  3674.0000     1.0000   752.4348 

int.growth.common( 272, 78, 9, inf=int.s.inf.2, thr=80, min.s=55 )
##        q          m          n          k          p        exp 
## 204.0000  3293.0000 12445.0000  2920.0000     1.0000   610.9773 

## so those that have got larger are not the same set of introns

## Some of these numbers are a little unexpected,
## but reasonable if we look at the following plots.
par(mfrow=c(1,3))
plot(int.s.inf.2[[272]]$state[,1], int.s.inf.2[[78]]$state[,1], cex=0.75, xlim=c(55, 180), ylim=c(55, 180),
     col=rgb(0,0,0,0.05), pch=19)
abline(0,1, col='red'); abline(h=80, v=80, col='red')
plot(int.s.inf.2[[272]]$state[,1], int.s.inf.2[[77]]$state[,1], cex=0.75, xlim=c(55, 180), ylim=c(55, 180),
      col=rgb(0,0,0,0.05), pch=19)
abline(0,1, col='red'); abline(h=80, v=80, col='red')
plot(int.s.inf.2[[272]]$state[,1], int.s.inf.2[[76]]$state[,1], cex=0.75, xlim=c(55, 180), ylim=c(55, 180),
      col=rgb(0,0,0,0.05), pch=19)
abline(0,1, col='red'); abline(h=80, v=80, col='red')


## we have to do some simulation to make sure that these p-values
## actually mean something. Simulate the random situation from an
## ancestor and some descendants. Then reinfer the tree using
## the same meethod that we have here and then redo.. 

## 1. determine a model for how size has changed..
plot.is.d <- function(a, d, inf=int.s.inf.2, ylim=c(-100, 100)){
    d <- inf[[d]]$state[,1] - inf[[a]]$state[,1]
    plot(inf[[a]]$state[,1], d, cex=0.75, col=rgb(0,0,0,0.05), pch=19, ylim=ylim)
}

par(mfrow=c(1,3))
x <- -150:150
plot.is.d( 272, 78, inf=int.s.inf.2 )
lines(x, 10*log2( 2^(x/10) + 2000 ) - x, col='red')
abline(v=80, col='red')
plot.is.d( 272, 77, inf=int.s.inf.2 )
lines(x, 10*log2( 2^(x/10) + 1000 ) - x, col='red')
abline(v=80, col='red')
plot.is.d( 272, 76, inf=int.s.inf.2 )
lines(x, 10*log2( 2^(x/10) + 1000 ) - x, col='red')
abline(v=80, col='red')

## also consider the common ancestor of P. natter. and A. mexic.
## #304

plot.is.d( 304, 76, inf=int.s.inf.2 )
lines(x, 10*log2( 2^(x/10) + 1000 ) - x, col='red')
##
plot.is.d( 304, 77, inf=int.s.inf.2 )
lines(x, 10*log2( 2^(x/10) + 1000 ) - x, col='red')
##
plot.is.d( 272, 77, inf=int.s.inf.2 )
lines(x, 10*log2( 2^(x/10) + 1000 ) - x, col='red')
##
## probability of change increases with time, but the
## magnitude does not seem to change much..

require(entropy) ## for discretize2d; but turns out better to make a different one
## has a blur function, but seems to be otherwise a massive package
require(spatstat)  ## but it is not good for this

## instead use my own convolve function
dyn.load("~/R/c/convolve_ks.so")

diff.2d <- function(anc, desc, min.s=10*log2(76),
                    inf=int.s.inf.2, blur.sd=6, blur.w=10,
                    abline.args=NULL, do.plot=TRUE
                    ){
    int.s <- cbind(inf[[anc]]$state[,1], inf[[desc]]$state[,1])
    b <- apply(int.s, 1, min) > min.s
    int.s <- int.s[b,]
    int.d <- int.s[,2] - int.s[,1]
    int.s.r <- range(int.s)
    int.d.r <- range(int.d)
    dummy <- suppressWarnings( cbind(int.s.r[1]:int.s.r[2], int.d.r[1]:int.d.r[2] ) )
    dummy.tbl <- table(dummy[,1], dummy[,2])
    s.d <- table( c(dummy[,1], int.s[,1]), c(dummy[,2], int.d) )
    s.d <- s.d - dummy.tbl
    s1.d <- table( c(int.d.r[1]:int.d.r[2], int.d) ) - 1
    s1.d <- s1.d / sum(s1.d)
    ## s.d <- discretize2d( int.s[,1], int.d, numBins1=(max(int.s[,1])-min(int.s[,1])),
    ##                     numBins2=(max(int.d) - min(int.d)))
    kernel <- dnorm(0:blur.w, mean=0, sd=blur.sd)
    v.sd <- apply(s.d, 2, function(x){ .Call('convolve_ks', as.double(x), kernel) })
    s.d.blurred <- t(apply(v.sd, 1, function(x){ .Call('convolve_ks', as.double(x), kernel) }))
##    s.d.blurred <- v.sd + h.sd
    s.d.blurred <- s.d.blurred / rowSums(s.d.blurred)
    par(mfrow=c(2,3))
    x <- int.s.r[1]:int.s.r[2]
    y <- int.d.r[1]:int.d.r[2]
    if(do.plot){
        plot(int.s[,1], int.d, cex=0.5, col=rgb(0,0,0,0.1))
        if(!is.null(abline.args)) do.call( abline, abline.args )
        image(x=x, y=y, s.d )
        if(!is.null(abline.args)) do.call( abline, abline.args )
        image(x=x, y=y, s.d.blurred )
        if(!is.null(abline.args)) do.call( abline, abline.args )
        ## image(x=x, y=y, log(s.d.blurred) )
        ## if(!is.null(abline.args)) do.call( abline, abline.args )
        image(x=x, y=y, t(scale(t(s.d), center=FALSE)))
        ##    image(x=x, y=y, log( s.d / rowSums(s.d) ))
        if(!is.null(abline.args)) do.call( abline, abline.args )
        image(x=x, y=y, scale(s.d))
        if(!is.null(abline.args)) do.call( abline, abline.args )
        image(x=x, y=y, t(scale(t(s.d.blurred))))
        if(!is.null(abline.args)) do.call( abline, abline.args )
    }
    invisible(list(x=x, y=y, sh=s.d, bl=s.d.blurred, pt=cbind(as=int.s[,1], delta=int.d)))
}

diff.triplet <- function(anc, d1, d2, min.s=10*log2(76),
                         inf=int.s.inf.2, blur.sd=6, blur.w=10,
                         abline.args=NULL){
    int.s <- sapply( c(anc, d1, d2), function(x){ inf[[x]]$state[,1] })
    int.s <- int.s[ apply(int.s, 1, min) > min.s, ]
    int.d1 <- int.s[,2] - int.s[,1]
    int.d2 <- int.s[,3] - int.s[,1]
    d1.r <- range(int.d1)
    d2.r <- range(int.d2)
    dummy <- suppressWarnings( cbind(d1.r[1]:d1.r[2], d2.r[1]:d2.r[2] ) )
    dummy.tbl <- table(dummy[,1], dummy[,2])
    s.d <- table( c(dummy[,1], int.d1), c(dummy[,2], int.d2) )
    s.d <- s.d - dummy.tbl
    kernel <- dnorm(0:blur.w, mean=0, sd=blur.sd)
    v.sd <- apply(s.d, 2, function(x){ .Call('convolve_ks', as.double(x), kernel) })
    s.d.blurred <- t(apply(v.sd, 1, function(x){ .Call('convolve_ks', as.double(x), kernel) }))
##    s.d.blurred <- s.d.blurred / rowSums(s.d.blurred)
    x <- d1.r[1]:d1.r[2]
    y <- d2.r[1]:d2.r[2]
    par(mfrow=c(1,3))
    plot( int.d1, int.d2, cex=0.5, col=rgb(0,0,0,0.1))
    if(!is.null(abline.args)) do.call(abline, abline.args)
    image( x=x, y=y, s.d )
    if(!is.null(abline.args)) do.call(abline, abline.args)
    image( x=x, y=y, s.d.blurred )
    if(!is.null(abline.args)) do.call(abline, abline.args)
}

## make a simulate function on the basis of something like this:
simulate.is.evol <- function(anc, desc, gen.n, min.s=10*log2(76),
                             inf=int.s.inf.2, blur.sd=6, blur.w=10,
                             size.depend=TRUE, norm.sd=NA, use.blurred=TRUE){
    int.s <- cbind(inf[[anc]]$state[,1], inf[[desc]]$state[,1])
    b <- apply(int.s, 1, min) > min.s
    int.s <- int.s[b,]
    int.d <- int.s[,2] - int.s[,1]
    int.s.r <- range(int.s)
    int.d.r <- range(int.d)
    dummy <- suppressWarnings( cbind(int.s.r[1]:int.s.r[2], int.d.r[1]:int.d.r[2] ) )
    dummy.tbl <- table(dummy[,1], dummy[,2])
    s.d <- table( c(dummy[,1], int.s[,1]), c(dummy[,2], int.d) )
    s.d <- s.d - dummy.tbl
    s1.d <- table( c(int.d.r[1]:int.d.r[2], int.d) ) - 1
    s1.d <- s1.d / sum(s1.d)
    ## s.d <- discretize2d( int.s[,1], int.d, numBins1=(max(int.s[,1])-min(int.s[,1])),
    ##                     numBins2=(max(int.d) - min(int.d)))
    kernel <- dnorm(0:blur.w, mean=0, sd=blur.sd)
    v.sd <- apply(s.d, 2, function(x){ .Call('convolve_ks', as.double(x), kernel) })
    s.d.blurred <- t(apply(v.sd, 1, function(x){ .Call('convolve_ks', as.double(x), kernel) }))
##    h.sd <- t(apply(s.d, 1, function(x){ .Call('convolve_ks', as.double(x), kernel) }))
##    s.d.blurred <- v.sd + h.sd
    s.d.blurred <- s.d.blurred / rowSums(s.d.blurred)
    s.d[ rowSums(s.d) == 0, ] <- 1
    s.d <- s.d / rowSums(s.d)
    par(mfrow=c(2,3))
    x <- int.s.r[1]:int.s.r[2]
    y <- int.d.r[1]:int.d.r[2]
    ## image(x=x, y=y, s.d )
    ## image(x=x, y=y, s.d.blurred )
    ## image(x=x, y=y, t(scale(t(s.d))))
    ## image(x=x, y=y, t(scale(t(s.d.blurred))))
    ## we now consider how to create a new set of values...
    ## dd.set <- c(x1[1], x2)
    ## evolve by sampling from the distributions:
    make.descendants <- function(anc.i, gen, max.i){
        if(gen > gen.n)
            return(max.i)
        a <- nodes[[anc.i]]$size
        a[ a < min(x) ] <- min(x)
        a[ a > max(x) ] <- max(x)
        if(size.depend){
            if(use.blurred){
                d.d <- sapply( a, function(l){
                    sample( y, size=2, replace=TRUE, prob=s.d.blurred[ 1 + l - x[1], ] )})
            }else{
                d.d <- sapply( a, function(l){
                    sample( y, size=2, replace=TRUE, prob=s.d[ 1 + l - x[1], ] )})
            }                
        }else{
            if(is.na(norm.sd)){
                d.d <- matrix( sample(y, size=2*length(a), replace=TRUE, prob=s1.d), nrow=2 )
            }else{
                d.d <- matrix( rnorm(n=2*length(a), sd=norm.sd), nrow=2 )
            }
        }
        d1 <- a + d.d[1,]
        d2 <- a + d.d[2,]
        n1 <- max.i + 1
        n2 <- max.i + 2
        print(paste(max.i, anc.i, n1, n2, length(nodes)))
        nodes[[n1]] <<- list(parent=anc.i, node=n1, size=d1, gen=gen+1)
        nodes[[n2]] <<- list(parent=anc.i, node=n2, size=d2, gen=gen+1)
        max.i <- make.descendants( n1, gen+1, n2 )
        make.descendants( n2, gen+1, max.i )
    }
    nodes <- vector(mode='list', length=sum( 2^(0:gen.n) ))
    nodes[[1]] <- list(parent=1, node=1, size=int.s[,1], gen=1)
    make.descendants( 1, 1, 1 )
    nodes
}

descendants <- simulate.is.evol(anc=272, desc=78, gen.n=3, blur.sd=12, blur.w=20)
desc.2 <- simulate.is.evol(anc=272, desc=78, gen.n=3, blur.sd=12, blur.w=20, size.depend=FALSE)
desc.3 <- simulate.is.evol(anc=272, desc=78, gen.n=3, blur.sd=12, blur.w=20, size.depend=FALSE, norm.sd=10)

### 288 is the common ancestor of B. splendens (3) and Takifugu (2).

desc.4 <- simulate.is.evol(anc=288, desc=2, gen.n=3, blur.sd=12, blur.w=20)
desc.5 <- simulate.is.evol(anc=288, desc=2, gen.n=3, blur.sd=12, blur.w=20, size.depend=FALSE)
desc.6 <- simulate.is.evol(anc=288, desc=2, gen.n=3, blur.sd=12, blur.w=20, size.depend=FALSE, norm.sd=20)

desc.7 <- simulate.is.evol(anc=266, desc=2, gen.n=3, blur.sd=12, blur.w=20)
desc.8 <- simulate.is.evol(anc=266, desc=3, gen.n=3, blur.sd=12, blur.w=20)
desc.9 <- simulate.is.evol(anc=266, desc=4, gen.n=3, blur.sd=12, blur.w=20)
desc.10 <- simulate.is.evol(anc=266, desc=9, gen.n=3, blur.sd=12, blur.w=20)
desc.11 <- simulate.is.evol(anc=266, desc=10, gen.n=3, blur.sd=12, blur.w=20)

## blurred
desc.bl.266 <- lapply( c(276, 267, 2, 3, 4, 9, 10, 78), function(x){
    simulate.is.evol(anc=266, desc=x, gen.n=3, blur.sd=12, blur.w=20)
})

## not blurred
desc.266 <- lapply( c(276, 267, 2, 3, 4, 9, 10, 78), function(x){
    simulate.is.evol(anc=266, desc=x, gen.n=3, blur.sd=12, blur.w=20, use.blurred=FALSE)
})


int.growth.common.2( 1, 6, 12, descendants, thr=80, min.s=55 )
##            q            m            n            k            p          exp 
## 4.861000e+03 8.514000e+03 6.515000e+03 8.434000e+03 2.781389e-03 4.777901e+03 
##            r 
## 1.017392e+00 

int.growth.common.2( 1, 6, 12, desc.2, thr=80, min.s=55 )
##            q            m            n            k            p          exp 
## 8.096000e+03 1.136100e+04 4.702000e+03 1.137200e+04 2.113596e-02 8.043161e+03 
##            r 
## 1.006569e+00 

int.growth.common.2( 1, 6, 12, desc.3, thr=80, min.s=55 )
##            q            m            n            k            p          exp 
## 2.113000e+03 5.475000e+03 8.645000e+03 5.300000e+03 1.860824e-02 2.055064e+03 
##            r 
## 1.028192e+00 

## These give me significant p-values, but I don't know why.
## The expected / observed are however very low.. 

length.retained.2(1, 6, 12, descendants, thr=80, min.s=55 )
##            q            m            n            k            p          exp 
## 1.165600e+04 1.288200e+04 1.455000e+03 1.286200e+04 1.154716e-17 1.155669e+04 
##            r 
## 1.008593e+00 

length.retained.2(1, 6, 12, desc.2, thr=80, min.s=55 )
##            q            m            n            k            p          exp 
## 1.236600e+04 1.332000e+04 1.017000e+03 1.326600e+04 6.893880e-07 1.232497e+04 
##            r 
## 1.003329e+00

length.retained.2(1, 6, 12, desc.3, thr=80, min.s=55 )
##            q            m            n            k            p          exp 
## 9.466000e+03 1.150100e+04 2.836000e+03 1.137100e+04 5.693745e-66 9.121704e+03 
##            r 
## 1.037745e+00 

## to plot the sizes going from ancestor to descendant..
par(mfrow=c(2,4))
invisible( sapply( c(1,2,4,6), function(i){ hist(descendants[[i]]$size, breaks=100) }))
invisible( sapply( c(1,2,4,6), function(i){ hist(desc.2[[i]]$size, breaks=100) }))

## and now use the leaf nodes to infer the intron sizes of ancestral lineages
## For this we need:
dyn.load("~/R/R_max_parsimony/src/max_parsimony.so")
source("~/R/R_max_parsimony/R/functions.R")


## make a wrapper function to call sankoff on the simulated tree

infer.ancestors <- function(desc){
    ## we need an edge table, this should have:
    ## parent, node
    edge <- t(sapply(desc, function(x){ c(x$parent, x$node) }))
    ## remove the first circular root node..
    ## alternatively modify the simulate.is.evol function
    edge <- edge[ edge[,1] != edge[,2], ]
    ## leaf nodes need to be numbered 1..leaf.n
    ## in order that they can correspond to a table of
    ## of states
    leaf.i <- which( !(edge[,2] %in% edge[,1]) )
    ## the rearranged edge edge.ra
    edge.ra <- edge + length(leaf.i)
    edge.ra[ leaf.i, 2 ] <- 1:length(leaf.i)
    edge.ra
    ##  I also need to make sure that the index of the nodes does not
    ##  exceed the number of nodes..
    total.nodes <- length(unique(as.integer(edge.ra)))
    node.remap <- rep(NA, max(edge.ra))
    node.remap[ sort(unique(as.integer(edge.ra))) ] <- 1:total.nodes
    edge.ra <- cbind( node.remap[ edge.ra[,1] ], node.remap[ edge.ra[,2] ])
    ## then obtain the leaf states and prepare the alphabet.
    states <- sapply( desc[leaf.i + 1], function(x){ x$size })
    al.size <- as.integer( 1 + max(states))
##    al.size <- as.integer( 1 + diff( range(states) ))
## I should rewrite the max parsimony function so that it can handle
## non 0 starts.
    al.offset <- 1L
    sub.matrix <- make.sub.matrix( al.size )
    missing.penalty <- as.integer( max(sub.matrix) + 1 )
    sub.matrix[1,] <- missing.penalty
    sub.matrix[,1] <- missing.penalty
    states.enc <- encode.dist( states, offset=al.offset )
    leaf.n <- sum( !(edge[,2] %in% edge[,1]) )
    total.nodes <- length(unique(as.integer(edge)))
    inf <- sankoff( edge.ra, c(total.nodes, leaf.n), sub.matrix,
                   c(al.offset, al.size), states.enc, check.states=FALSE)
    list(inf=inf, leaves=states)
}

d1.inf <- infer.ancestors(descendants)
## d2.inf <- infer.ancestors(desc.2) ## this seems to have too many distinct states
d3.inf <- infer.ancestors(desc.3)
## we can then ask questions of the descendants.

int.growth.common( 9, 1, 5, inf=d1.inf$inf, thr=80, min.s=55 )
##        q         m         n         k         p       exp 
## 495.0000 2019.0000 4338.0000 2012.0000    1.0000  639.0165 

int.growth.common( 9, 1, 5, inf=d3.inf$inf, thr=80, min.s=55 )
##       q         m         n         k         p       exp 
## 874.000  4234.000 10673.000  4245.000     1.000  1205.697 

length.retained( 9, 1, 5, inf=d1.inf$inf, thr=80, min.s=55 )
##            q            m            n            k            p          exp 
## 39842.000000 45371.000000  6145.000000 45395.000000     1.000000 39980.133260 
##            r 
##     0.996545 

length.retained( 9, 1, 5, inf=d3.inf$inf, thr=80, min.s=55 )
##            q            m            n            k            p          exp 
## 2.858000e+04 3.420600e+04 6.471000e+03 3.406800e+04 9.939571e-01 2.864838e+04 
##            r 
## 9.976133e-01 

int.growth.common( 9, 1, 6, inf=tmp$inf, thr=80, min.s=55 )
##        q         m         n         k         p       exp 
## 496.0000 2041.0000 4275.0000 2021.0000    1.0000  653.0812

## but this one is significant as it should be.
int.growth.common( 9, 1, 2, inf=tmp$inf, thr=80, min.s=55 )
##            q             m             n             k             p 
## 1.259000e+03  2.172000e+03  4.154000e+03  2.164000e+03 2.328655e-180 
##          exp 
## 7.429984e+02 


### how much of a problem is it that introns that are long
### in an ancestor are less likely to become short (< 256bp)
### in a descendant when we are considering the commonality
### of retention. It ought to be a problem (as suggested
### by the simulated data, but since the forbidden
### set is the same in each case this may not be that
### much of an issue (as it increases the implicit sample
### size.

selection.bias <- function(loss.max=0.75, thr=0.5, retain.prob=0.5, scale=1000){
    ## considering, N as the total number of points considered
    ## as being lost; i.e. those that are larger in the ancestor
    N <- 1 - thr
    ## m and k are simply retain.prob * (loss.max - thr)
    ## that is the number of introns that remain long in the
    ## descendants.
    m <- retain.prob * (loss.max - thr) 
    k <- retain.prob * (loss.max - thr)
    M <- m + (1-loss.max)
    K <- k + (1-loss.max)
    ## n is then N - m
    n <- N - M
    ## q is the proportion that we expect to be retained by
    ## both, plus the set
    q <- retain.prob * m + (1-loss.max)
    ## These numbers can then be scaled up to give
    p <- phyper( q * scale - 1, M * scale, n * scale, K * scale, lower.tail=FALSE )
    exp <- K * M / N
    c(q=q, m=M, n=n, k=K, p=p, exp=exp, r=q/exp)
}

## it does look like we may have an issue. To deal with this, the simplest
## thing is to divide up into separate sections..

common.retention <- function(anc, d1, d2, inf=int.s.inf.2, thr=80, breaks=seq(80, 200, 10), min.l=55){
    int.l <- sapply( inf[c(anc, d1, d2)], function(x){ x$state[,1] })
    int.l <- int.l[ apply(int.l, 1, min) >= min.l, ]
    sapply( 2:length(breaks), function(i){
        b <- int.l[,1] >= breaks[i-1] & int.l[,1] < breaks[i]
        q <- sum( b & int.l[,2] >= thr & int.l[,3] >= thr )
        m <- sum( b & int.l[,2] >= thr )
        n <- sum( b ) - m
        k <- sum(b & int.l[,3] >= thr )
        exp <- k * m / sum(b)
        r <- q/exp
        p <- phyper( q-1, m, n, k, lower.tail=FALSE )
        c(q=q, m=m, n=n, k=k, exp=exp, r=r, p=p)
    })
}

## for simulated data
common.retention.s <- function(anc, d1, d2, inf=descendants, thr=80, breaks=seq(80, 200, 10), min.l=55){
    int.l <- sapply( inf[c(anc, d1, d2)], function(x){ x$size })
    int.l <- int.l[ apply(int.l, 1, min) >= min.l, ]
    sapply( 2:length(breaks), function(i){
        b <- int.l[,1] >= breaks[i-1] & int.l[,1] < breaks[i]
        q <- sum( b & int.l[,2] >= thr & int.l[,3] >= thr )
        m <- sum( b & int.l[,2] >= thr )
        n <- sum( b ) - m
        k <- sum(b & int.l[,3] >= thr )
        exp <- k * m / sum(b)
        r <- q/exp
        p <- phyper( q-1, m, n, k, lower.tail=FALSE )
        c(q=q, m=m, n=n, k=k, exp=exp, r=r, p=p)
    })
}


## 266 us the common ancestor of (d.rerio, i.punctatus, e.electro, a.mex, ...)
## and most of the others (b. splendens, p. ranga, tetraodon, takifugu, ...)
## 1: tetraodon, 2, takifugu, 3, b. splendens,
## 78: D. rerio, 9. E. electricus, 30. I.punctatus
common.retention( 266, 78, 1, breaks=seq(80, 200, 10) )

## I. punctatus & Takifugu
common.retention( 266, 30, 2, breaks=seq(80, 200, 10) )

## I. punctatus & B. splendens
common.retention( 266, 30, 3, breaks=seq(80, 200, 10) )

## I. punctatus & Tetraodon
common.retention( 266, 30, 1, breaks=seq(80, 200, 10) )

## I. punctatus & Gasterosteus
common.retention( 266, 30, 4, breaks=seq(80, 200, 10) )

## E. electricus & Tetraodon
common.retention( 266, 9, 1, breaks=seq(80, 200, 10) )

## E. electricus & B. splendens
common.retention( 266, 9, 3, breaks=seq(80, 200, 10) )

## D. clupeoides & B. splendens
common.retention( 266, 10, 3, breaks=seq(80, 200, 10) )


common.retention.s( 1, 6, 12, inf=descendants, breaks=seq(80, 200, 20) )
common.retention.s( 1, 6, 12, inf=desc.2, breaks=seq(80, 200, 20) )
common.retention.s( 1, 6, 12, inf=desc.3, breaks=seq(80, 200, 20) )
common.retention.s( 1, 6, 12, inf=desc.4, breaks=seq(80, 200, 20) )
common.retention.s( 1, 6, 12, inf=desc.5, breaks=seq(80, 200, 20) )
common.retention.s( 1, 6, 12, inf=desc.6, breaks=seq(80, 200, 20) )

## the distributions of differences are:
diff.2d( 266, 276 )
diff.2d( 266, 267 )
diff.2d( 263, 266 )
diff.2d( 266, 282 )
diff.2d( 266, 2, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 3, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 4, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 26, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 9, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 10, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 24, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 28, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 30, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 31, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 38, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 266, 78, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )

sd.266 <- lapply( c(276, 267, 282, 2, 3, 4, 5, 6, 9, 10, 24, 28, 30, 31, 38, 78),
                 function(x){ diff.2d(266, x, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' )) })

## the following are based on 266 -> 276, 267, 2,3,4,9,10,78 (desc.7-11)
## and are probably more reasonable.
lapply( desc.bl.266, function(x){ common.retention.s( 1, 6, 12, inf=x, breaks=seq(80, 200, 10)) })
lapply( desc.bl.266, function(x){ common.retention.s( 1, 2, 3, inf=x, breaks=seq(80, 200, 10)) })
lapply( desc.266, function(x){ common.retention.s( 1, 6, 12, inf=x, breaks=seq(80, 200, 10)) })
lapply( desc.266, function(x){ common.retention.s( 1, 2, 3, inf=x, breaks=seq(80, 200, 10)) })

## low p-values for the 1-> 2,3 pairing or when using the 266-> 276,267 distribution

lapply( c(9, 10, 24, 28, 30, 31, 38, 78), function(x){
    common.retention(266, 2, x, breaks=seq(80, 200, 10))
})


### lets all teleost ancestor -> leaf relationships
teleost.int.evol <- lapply( which(sp.class[,'teleostei']), function(i){
    diff.2d( 266, i, do.plot=FALSE )
})
names(teleost.int.evol) <- ex.align.2.k2.nj$tip.label[ which(sp.class[,'teleostei']) ]

## ideally we would put these in taxonomic order, but for the moment I think that they are
## in order of genome size.
mi2col <- function(v, cm){
    hsvScale(v, val=(cm + v)/(cm + max(v, na.rm=TRUE)))
}

par(mfrow=c(3, 5))
invisible( sapply( names(teleost.int.evol)[1:15], function(nm){
    with(teleost.int.evol[[nm]], image(x, y, t(scale(t(sh))), main=nm,
                                       col=mi2col(1:1000, 300),
                                       xlab='Ancestral size',
                                       ylab='delta'))
    abline(log2(76)*10, -1, col='white', lty=3)
}))

## 259 is the root of all jawed vetebrates
diff.2d( 259, 26, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 259, 78, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 259, 49, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
## zebra finch:
diff.2d( 259, 75, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
## mus musculus
diff.2d( 259, 127, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
## homo sapiens
diff.2d( 259, 155, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )

## the common ancestor of birds and crocodiles: 252 and birds:
diff.2d( 252, 61, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )

## gallus gallus: 61, against L. oculatus.. 49
diff.2d( 49, 61,  abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )

## 263 is the common ancestor of L oculatus (49) and the teleosts
## 26 is Medaka
diff.2d( 263, 26, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 263, 49, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 263, 78, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 263, 77, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 263, 76, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
## 105 hucho hucho
diff.2d( 263, 105, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 263, 44, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 263, 10, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )
diff.2d( 263, 9, abline.args=list(log2(76)*10, -1, v=c(80, 100), col='blue' ) )

## try some other stuff.. 
diff.triplet( 266, 2, 9 )
diff.triplet( 266, 3, 9 )

diff.triplet( 266, 3, 78, abline.args=list(0, 1, v=0, h=0, col='blue', lty=3) )
diff.triplet( 266, 26, 78, abline.args=list(0, 1, v=0, h=0, col='blue', lty=3) )
diff.triplet( 266, 26, 10, abline.args=list(0, 1, v=0, h=0, col='blue', lty=3) )
