# Evidence supporting transcripts and intron positions

The question as to what extent the similarity of intron lengths in closely related species is an artefact
of the annotation process has been raised by one of the peer reviewers. Although this is unlikely to be
an issue since transcript prediction is heavily exon biassed and intron size changes are likely to be primarily
caused by mutations affecting the intron sequence rather than those affecting the exons (as implied by the
exponential nature of intron length distributions), it would be desirable to:

a. establish what proportions of transcripts are supported by within species RNA sequences
b. establish the nature of the distributions of orthologous introns across species.

In this directory we will look at question q) using the Ensembl databases.

## Uniprot -> species
We have mappings of Uniprot ids to NCBI taxonomy ids in:

/home/lmj/genomes/uniprot/idmapping_selected.tab.gz

This has columns:

	1. UniProtKB-AC
	2. UniProtKB-ID
	3. GeneID (EntrezGene)
	4. RefSeq
	5. GI
	6. PDB
	7. GO
	8. UniRef100
	9. UniRef90
	10. UniRef50
	11. UniParc
	12. PIR
	13. NCBI-taxon
	14. MIM
	15. UniGene
	16. PubMed
	17. EMBL
	18. EMBL-CDS
	19. Ensembl
	20. Ensembl_TRS
	21. Ensembl_PRO
	22. Additional PubMed

Which means we want to map column 1 to column 13 and then get the node name from:

/home/lmj/genomes/ncbi_taxonomy_2021_09/names.dmp

Get the relevant ids from the files in dna_align_features

use:
protein_taxon.pl

## EMBL -> taxon id

For this we need to look up identifiers in:

/home/lmj/genomes/ncbi_taxonomy_2021_09/nucl_gb.accession2taxid.gz

Which and then map to the species ids.

use
gb_taxon.pl
