require("RMySQL")
source("~/R/general_functions.R")
source("../R_common_functions/functions.R")

## get the set of installed databases from
tmp <- read.table( "../core_dbs/ensembl_98_all_dbs", stringsAsFactors=FALSE)
db.list <- tmp[,9]
rm(tmp)

## I only installed core databases; and not for subspecies

db.list <- db.list[ grep("^[^_]+_[^_]+_core", db.list) ]
## that gives us a 177 databaes

length(unique( sub("_core_\\d+_\\d+$", "", db.list) ))
## from 177 different species!

## we can extract the species names
db.name <- sub("_core_\\d+_\\d+$", "", db.list)
db.sp <- sub("_core_\\d+_\\d+$", "", db.list)
db.sp <- sub("_", " ", db.sp)
substr(db.sp, 1, 1) <- toupper( substr(db.sp, 1, 1))
db.genus <- sub(" [^ ]+$", "", db.sp )

length(unique(db.genus))
## 156
## so I have a fair number of closely related species here as well.

## It seems that it would be a good idea to map these identifiers to some sort of phylogenetic tree
## We can use some of the functions that I developed in the previous analysis.

ens.tree <- get.ensembl.taxonomy()
sum( db.sp %in% ens.tree[,'node_name'] )
## 170, something matching..

db.sp[ which( !(db.sp %in% ens.tree[,'node_name'] ) ) ]
## [1] "Caenorhabditis elegans"   "Canis familiaris"        
## [3] "Cebus capucinus"          "Gorilla gorilla"         
## [5] "Mus caroli"               "Mus pahari"              
## [7] "Mus spretus"              "Saccharomyces cerevisiae"

ens.tree[ grep("Saccharomyces", ens.tree[,'node_name']), ]
ens.tree[ grep("Mus caroli", ens.tree[,'node_name']), ]
ens.tree[ grep("Gorilla gorilla", ens.tree[,'node_name']), ]

## necessary for other reasons
db.sp[ db.sp == 'Canis familiaris'] <- 'Canis lupus'

## the aboe are not equal, so instead we want to do:
db.tree.i <- sapply( db.sp, function(x){ grep(x, ens.tree[,'node_name'])[1] })
db.tree.sp <- ens.tree[ db.tree.i, 'node_name'] ## I need these for the lineage query..

db.lineages <- lapply( db.tree.sp, function(x){ get.ensembl.lineage( ens.tree, x )})

## what sort of ranks can we use?
sort(table( unlist( lapply( db.lineages, function(x){ x$rank }))))
## so we have some databases that do not even have a phylum associated
## with their taxonomic tree. Which ones?
 ## subspecies    subgenus superfamily       tribe   parvorder  infraorder 
 ##          2           9          10          17          20          25 
 ## infraclass       genus      cohort    subclass   subfamily      family 
 ##         32          48          52          55          68          71 
 ##   suborder       class       order  superclass  superorder      phylum 
 ##         82         102         131         169         172         174 
 ##    species     no rank 
 ##        175        1855 


which(sapply( db.lineages, function(x){ sum("phylum" %in% x$rank) == 0} ))
## 18, 49, 150

db.lineages[[18]]  ## C. elegans
db.lineages[[49]]  ## D. melanogaster
db.lineages[[150]] ## S. cerevisiae

## We can leave those out in any case..

db.sel <- sapply( db.lineages, function(x){ sum("phylum" %in% x$rank) > 0} )
## 174 databases

## extract the phylum from these
db.sel.phylum <- sapply( db.lineages[db.sel], function(x){ x[ x$rank == 'phylum', 'node_name' ] })
## these are all chordates. I think we can go with that and then simply have a look at what happens
## afterwards. At the moment I anyway have hagfish in my set of data, so it makes sense to keep this
## as it is.

## let us get some genome statistics, as the may be useful for checking things with:
genome.stats <- lapply( db.list, get.genome.stats )
genome.sizes <-  sapply( genome.stats, function(x){ x[ x$statistic == 'ref_length', 'value']  })


### Let us get the family members from all species
family.members <- lapply( db.name[db.sel], function(x){
    print(paste(x))
    get.family.members(x, db="ensembl_compara_98", db.is.name=TRUE )
})

db <- ensembl.connect("ensembl_compara_98")
family <- dbGetQuery( db, "select * from family" )
dbDisconnect(db)

## Get genes that map to non-variant scaffolds or chromosomes. This in order
## to avoid overcounting family members
ref.genes <- lapply( db.list[db.sel], get.all.ref.genes, max.coord.rank=2 )
names(ref.genes) <- db.list[db.sel]

## note that the lengths of these are very variable, with human having loads and loads
## of entries; However, the number of protein coding genes which we are using (protein families)
## are only around 20,000, so not so bad.

sapply(1:length(family.members), function(i){
    c(length(unique(family.members[[i]]$gene.id)), sum( unique(family.members[[i]]$gene.id) %in% ref.genes[[i]]$stable_id) ) })
## this gives us a big table that doesn't tell me that much.. 

## let us remove alternate genes

for(i in 1:length(family.members)){
     family.members[[i]] <- family.members[[i]][ family.members[[i]][,'gene.id'] %in% ref.genes[[i]][,'stable_id'], ]
}

sapply( family.members, nrow ) ## this gives rather large numbers of nrows; but not so many genes involved

## let us get a table of counts for each species..
family.counts <- sapply( family.members, function(x){
    table(c(family$stable_id, x$family.id[!duplicated( x$gene.id )])) - 1 })

dim(family.counts)
## [1] 32295     6

## we are more interested where there is only a single member of
## each family in each species. That makes life much easier.
family.counts.h <- hist( rowSums(family.counts == 1) )
family.counts.q <- sort( rowSums(family.counts == 1) )


## if we only look at vertebrate species
vert.b <- genome.sizes[db.sel] > 3e8
## this removes two ciona species from the analysis

family.counts.v.h <- hist( rowSums(family.counts[,vert.b] == 1), breaks=-1:sum(vert.b) + 0.5 )
## better perhaps to do:

fam.single.v <- table( rowSums(family.counts[,vert.b] > 0 & family.counts[,vert.b] == 1) )

##     0     1     2     3     4     5     6     7     8     9    10    11    12 
## 12105   997   473   306   251   250   198   207   176   173   158   156   146 
##    13    14    15    16    17    18    19    20    21    22    23    24    25 
##   140   145   140   118   113   124   123    95   101    94   101    93    96 
##    26    27    28    29    30    31    32    33    34    35    36    37    38 
##    81    81    79    92    88    88    89    84    90    74    68    82    83 
##    39    40    41    42    43    44    45    46    47    48    49    50    51 
##    65    71    88    92   101    78    85    72    92   107   110   115   105 
##    52    53    54    55    56    57    58    59    60    61    62    63    64 
##   117   106    98    78    79    69    68    76    61    70    62    77    69 
##    65    66    67    68    69    70    71    72    73    74    75    76    77 
##    70    79    87    83    77    98    82    79    82    80    87    74    81 
##    78    79    80    81    82    83    84    85    86    87    88    89    90 
##    61    76    68    45    66    79    56    71    77    74    68    68    54 
##    91    92    93    94    95    96    97    98    99   100   101   102   103 
##    69    69    67    86    77    59    78    70    85    84    93   111   112 
##   104   105   106   107   108   109   110   111   112   113   114   115   116 
##    97   128    99   130   131   112   158   149   155   132   146   147   124 
##   117   118   119   120   121   122   123   124   125   126   127   128   129 
##   113   120    90    86   102    93    77    84    57    67    81    85    67 
##   130   131   132   133   134   135   136   137   138   139   140   141   142 
##    82    83    77    80    88   105    97    90   127    99   107    84   104 
##   143   144   145   146   147   148   149   150   151   152   153   154   155 
##   104   112   114   105   137   139   160   160   180   185   188   229   226 
##   156   157   158   159   160   161   162   163   164   165   166   167   168 
##   247   256   298   242   274   276   293   233   245   167   138    84    59 
##   169   170   171   172 
##    21    14     4     1 

cumsum(rev(fam.single.v))
## 130 members out of 172 gives us 6114 families; we can probably use something
## like that for the analyses.
## 130/172 = 0.756
## 15/20   = 0.75
## so the proportion is similar to what I did for the smaller (15/20 species having
## a member
##
## We can start with this ..

fam.1.130.members <- lapply( rownames(family.counts)[ rowSums(family.counts[,vert.b] == 1) >= 130 ], function(fam){
    lapply( family.members[ vert.b ], function(x){
        unique( x$gene.id[ x$family.id == fam ] )
    })
})

names(fam.1.130.members) <- rownames(family.counts)[ rowSums(family.counts[,vert.b] == 1) >= 130 ]

## that actually allows for some gene duplications within a species; but most will be single genes. Let us output this list
## as a table
fam.1.130.tbl <- t(sapply( fam.1.130.members, function(x){
    sapply(x, paste, collapse=',') }))

colnames(fam.1.130.tbl) <- (db.list[ db.sel ])[ vert.b ]

write.table( fam.1.130.tbl, "vertebrate_family_members_1_130.txt", quote=FALSE, sep="\t" )
