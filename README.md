# Source and data files used to analyse vertebrate intron size evolution (ws2)

This repository contains a set of source code and data files used to
analyse how intron sizes have evolved during vertebrate evolution.

The source and data files are present throughout the directories, but are
described in two separate lists below. Many of the data files used as input
and listed below were created within a separate directory structure on a different
computer (ws1).

## Source files

### Project specific source code

- `R/`
  - `functions.R`  
  General plotting and data transformation functions. Somewhat redundant with
  those in `R_common_funnctions/` (see below).
- `R_172_genomes/`
  - `R_172_genomes.R`  
  Initial analysis of intron size changes across the vertebrates. Was used to
  create the intron orthology and Kimura two factor distances used for
  phylogenetic data. All processed transcript alignments were called from this
  file.

  The analyses were based on sequence and exon coordinates defined in analyses
  on ws1 in `R_i72_genomes.R` and `family_members/`
- `R_basic_stats/`
  - `basic_stats.R`  
  Statistics and plots based on all introns in all species without
  consideration of orthology.
- `R_common_functions/`
  - `functions.R`  
  Primarily database access functions, but also some general plotting functions.
  - `src/`  
  Compiled functions for extracting sequence masks from tabular blast
  data.
    - `bl_cov.c`
    - `test.R`
- `R_dr_repeats/`
  - `dr_repeats.R`  
  An analysis of repeat sequence influenece on intron size in *D. rerio*. 
- `R_gene_ontology/`
  - `gene_ontology.R`  
  Analyses used to investigate enrichment and depletion statistics for gene
  ontology group membership.
- `R_genome_pos/`
  - `genome_pos.R`  
  Relationship between genome location and intron size. Incomplete.
- `R_intron_alignments/`
  - `intron_alignments.R`  
  First attempts to assess the relationship between intron length and sequence conservation;
  the code was refactored in later analyses due to being too messy.
- `R_intron_alignments_2/`  
  Primarily a cleanup of the code in `R_intron_alignments`, but also includes
  blast of aligned sequences to a range of genomes (see subdirectories). Note that
  these analyses were carried out on (slightly) incorrectly defined sets of introns.
  - `functions.R`  
  Functions specific to intron alignments, sourced in several other R sessions.
  - `intron_alignments_2.R`
  - `extracted_introns/`  
    Directories and files containing extracted intron sequences and scripts for
    running blast on these against genome sequence databases.
    - `extract_introns.pl`
    - `extract_introns.sh`
    - `blast/`
      - `find_dbs.sh`
      - `run_blast.sh`
      - `run_single_blast.sh`
  - `extracted_seqs/`
    - `blast/`
      - `blast_alignments.R`
      - `find_dbs.sh`
      - `run_blast.sh`
      - `run_single_blast.sh`
- `R_intron_alignments_3/`
  - `intron_alignments_3.R`  
  Similar analyses to those in `R_intron_alignments_2`, but using
  *G. aculeatus* as the seed species.
- `R_intron_alignments_4/`
  - `intron_alignments_4.R`  
  Alignments as in prior analyses but using reverse complemented *D. rerio* sequences
  as control sequences.
- `R_intron_alignments_5/`  
  Inclusion of additional intron sets.
  - `intron_alignments_5.R`  
  - `extracted_seqs/`  
    Directories and files containing extracted intron sequences and scripts for
    running blast on these against genome sequence databases.
    - `blast/`
      - `find_dbs.sh`
      - `run_blast.sh`
      - `run_blast_2.sh`
      - `run_single_blast.sh`
- `R_intron_alignments_6/`  
  All alignments redone due to a bug found in the code selecting specific sets
  of introns which meant that the definitions of the different intron sets was
  not completely correct. Note that the error was minor, and that no
  difference in overall behaviour was observed.
  - `functions.R`
  - `intron_alignments_6.R`
- `R_intron_alignments_7/`  
  All orthologous introns from two closely related species aligned. 
  Incomplete analysis.
  - `functions.R`
  - `intron_alignments_7.R`
- `R_intron_alignments_8/`  
  Alignments of introns with a teleost minimal length between 256 and 90 bp long.
  - `intron_alignemnts_8.R`
- `R_intron_alignments_summary/`  
  A combined analyses of alignments in `R_intron_alignments_2` and
  `R_intron_alignments_5`. Note that these used the incorrectly defined intron sets.
  - `functions.R`
  - `intron_alignments_summary.R`
  - `blast/`  
    Directories and files containing extracted intron sequences and scripts for
    running blast on these against genome sequence databases.
    - `find_dbs.sh`
    - `run_blast.sh`
    - `run_blast_ctl.sh`
    - `run_blast_ctl_2.sh`
    - `run_blast_ctl_3.sh`
    - `run_blast_med.sh`
    - `run_single_blast.sh`
- `R_intron_alignments_summary_2/`  
  A combined analyses of alignments in `R_intron_alignments_6` and
  `R_intron_alignments_8` (i.e. the corrected intron sets).
  - `functions.R`  
    Local functions.
  - `intron_alignments_summary_2.R`
  - `blast/`  
    Directories and files containing extracted intron sequences and scripts for
    running blast on these against genome sequence databases.
    - `find_dbs.sh`
    - `run_blast.sh`
    - `run_ctl_blast.sh`
    - `run_ctl_blast_2.sh`
    - `run_short.2_blast.sh`
    - `run_short_blast.sh`
    - `run_single_blast.sh`
    - `run_single_blast_2.sh`
- `R_intron_orthology/`  
  Definition of a quality score for intron orthologue predictions.
  - `intron_orthology.R`
- `R_rand_alignments/`  
  Alignment of introns to random sequences in order to evaluate the alignment
  statistics used.
  - `rand_alignments.R`
- `R_trees_distances/`  
  Inference of ancestral intron sizes.
  - `trees_distances.R`
- `R_trees_distances_2/`  
  Inference of ancestral intron sizes but with finer grained intron size discretization.
  - `trees_distances_2.R`
- `danio_rerio_orth_introns/`  
  An aborted attempt to define repeats within *D. rerio* introns by
  self-blast. This turns out to be to too expensive, so I used Ensembl repeat
  annotation instead (see `R_dr_repeats/`).
  - `extract_dr_introns.pl`
  - `run_dr_blast.sh`
- `family_members/`  
  Files synchronised from ws1 (see descriptions above). Included in the
  directory structure since the resulting outputs were used as input data for
  downstream analyses.
  - `compile_seqs.pl`
  - `compile_seqs.sh`
  - `count_nucleotides.pl`
  - `family_members.R`

### External dependancies

A mixture of compiled (C) and R code used by the project specific
code. Note that the absolute locations of these files was different
when used for the described analyses, and that it is necessary to
amend all `source()` and `dyn.load()` statements in the above
described code in order to reflect the new location.

Both `R_max_parsimony` and `exon_aligneR` have their own github 
repositories; the ones given here reflect the code used in the
described analyses.

- `external/R_max_parsimony/`  
  Source for wrappers and compiled functions implementing Sankoff maximimum
  parsimony as anextension to `R`.
  - `README.md`
  - `R/`  
    `R` wrapper functions.
    - `functions.R`
  - `src/`
    - `max_parsimony.c`
	- `tree.c`
	- `tree.h`
- `external/exon_aligneR/`  
  Source for wrappers and compiled functions implementing a range of different
  pairwise alignment methods.
  - `README.md`
  - `functions.R`
  - `src/`
    - `cigar.h`
    - `exon_aligneR.c`
    - `needleman_wunsch.c`
    - `needleman_wunsch.h`
    - `util.c`
    - `util.h`
- `external/general_functions.R`
  A collection of small utility functions used by some of the analyses
  (primarily the `hsvScale()` colour generating function).


## Data files

Note that some of the files listed below have are not included in the repository
due to their excessive size. Nevertheless it should be possible to create most
of these files from the scripts included. In some cases these may require the
use of external databases such as genome blast databases.

### Compute node 2 (ws2)

- `R_172_genomes/`  
  Files resulting from initial analysis of orthologous introns from 172
  species, and data files obtained from `ws1`.
  - `assemblies.rds`  
    Copy of file from ws1. Not used due to variable representation.
  - `dr_intron_fam.txt`
  - `dr_intron_orthology_i.txt`
  - `dr_intron_orthology_id.txt`
  - `dr_intron_orthology_l.txt`
  - `dr_intron_orthology_tr.txt`  
    Files with names starting with `dr_intron_orthology` defined the intron
	orthology used here. Each one is a table with one column for each species
	and one row for each *D. rerio* intron found in gene members of the selected
	protein family set. The files contain:  `_i` intron rank; `_id` gene
	id, `_tr` transcript id, `_l` length. The table `dr_intron_fam.txt`
	contains a single column containing the family identifier for each row of
	these tables.
  - `dr_mammal_var.txt`
  - `dr_sauria_var.txt`
  - `dr_teleost_var.txt`  
    Tables with names matching the `dr_[a-z]+_var.txt` pattern contain intron
    length variance data for the respective species groups.
  - `ex_align_2.rds`  
    Alignment statistics derived from alignments of transcripts from 472 protein
	families used to define phylogenetic distances. A list based tree (protein family ->
	species 1 -> species 2 -> alignment statistics).
  - `ex_align_2_jc.rds`
  - `ex_align_2_jcg.rds`
  - `ex_align_2_id.txt`  
	Protein families and gene identities for genes aligned all-against-all in order
	to define phylogenetic distances.
  - `ex_align_2_k2.rds`  
    Files with names matching `ex_align_2_(jc|jcg|k2).rds` contain pairwise
    phylogenetic distances: `_jc` Jukes Cantor, `_jcg` Jukes Cantor gapped,
    `_k2` Kimura two factor.
  - `ex_align_2_k2_nj.rds`  
    A neighbour joining tree derived from the distances in `ex_align_2_k2.rds`
    using the `nj()` function from the `ape` package (CRAN).
  - `ex_align_id.rds`  
    Sequence alignments and statistics used to define the intron orthology. A
	list based tree structure (protein family -> seed species -> alignment information).  
	Note that this file has been split into five smaller files with suffixes `_aa` to `_ae`.
	The original file can be reformed by concatenation.
  - `genome_sizes.txt`
  - `int_s_inf.rds`  
    Inferred lengths (log~2~ transformed) of introns ancestral species.
  - `int_s_l.rds`  
    Intron lengths (log~2~ transformed).
  - `int_sp_mi_dr_m.rds`  
    Matrix of species mutual information statistics based on orthologous
    intron lengths.
  - `intron_files.txt`  
    A mapping of protein family identifiers to the files containing the intron sequences.
  - `leaf_int_d.rds`  
    Cumulative sums of mean intron size and phylogenetic distance from the
    common ancestor to each individual leaf node.
  - `mammal_b.txt`
  - `sauria_b.txt`
  - `teleost_b.txt`  
    Logical values; is the species a mammal, sauria, or teleost.
  - `ncbi_lineages.rds`  
    Taxonomic information obtained from the Ensembl compara database (ws1).
  - `seq_region.rds`  
    Contig and scaffold information for each species (ws1).
  - `species_lineages.rds`  
    The taxonomic lineage of individual species (i.e. the nodes traversed from
    the common root to each leaf).
- `R_basic_stats/`
  - `all_genes_exon_intron_stats.csv`  
    Copy of `family_members/all_genes/exon_intron_stats.csv` created on ws1.
- `R_dr_repeats/`
  - `danio_rerio_introns_mask.txt`  
    Repeat masking of *D. rerio* introns. 
	Incomplete analysis.
- `R_intron_alignments/`  
  Alignments of selected sets of *D. rerio* introns to orthologues in all
  other species.
  - `al_ctl_500.rds`
  - `al_ctl_500_ctl.rds`
  - `al_top_500.rds`
  - `al_top_500_ctl.rds`
- `R_intron_alignments_2/`  
  A set of files complementing the analysis in `R_intron_alignments`. Not used
  for the final analysis.
  - `aligns_top_al_cl.rds`
  - `extracted_introns/`
    - `long.ctl.txt`
    - `long.fasta`
    - `long.txt`
    - `sampled.ctl.txt`
    - `sampled.fasta`
    - `sampled.txt`  
	Intron sequences and coordinates from aligned sequences. Created by shell and
	perl scripts in the same directory. Not used.
    - `blast/`  
	  A directory containing series of blast output files created by the shell scripts in this directory.
	  Outputfiles not included due to size, and since they were not used.
  - `extracted_seqs/`  
	Sequence files containing the aligned parts of the *D. rerio* sequence which aligned with
	the highest score to sequences from the indicated clades (Mammalia, Sauria and Teleostei). These
	files were extracted by `intron_alignments_2.R`.
    - `mam_long.ctl.fasta`
    - `mam_long.fasta`
    - `mam_sampled.ctl.fasta`
    - `mam_sampled.fasta`
    - `sau_long.ctl.fasta`
    - `sau_long.fasta`
    - `sau_sampled.ctl.fasta`
    - `sau_sampled.fasta`
    - `tel_long.ctl.fasta`
    - `tel_long.fasta`
    - `tel_sampled.ctl.fasta`
    - `tel_sampled.fasta`
    - `blast/`  
	  Blast outpu files. Not included due to excessive size.
- `R_intron_alignments_5/`
  Similar to `R_intron_alignments_2` with a similar directory structure, but for an extended set of introns. Also
  not used in final analysis due to incorrect intron set selection.
  - `ctl_aligns.rds`
  - `ctl_aligns_top_cl.rds`
  - `extracted_seqs/`
    - `mam_ctl.fasta`
    - `mam_long_2.fasta`
    - `sau_ctl.fasta`
    - `sau_long_2.fasta`
    - `tel_ctl.fasta`
    - `tel_long_2.fasta`
    - `blast/`
  - `long_2_aligns.rds`
  - `long_3_aligns.rds`
- `R_intron_alignments_6/`  
  Repeat of alignments performed in `R_intron_alignments_2` and
  `R_intron_alignments_5`, but with the intron set memberships corrected. The `.rds` files
  here were used in `R_intron_alignments_summary_2`.
  - `ctl_aligns.rds`
  - `ctl_ctl_aligns.rds`
  - `ctl_i.rds`
  - `long_aligns.rds`
  - `long_ctl_aligns.rds`
  - `long_i.rds`
  - `med_aligns.rds`
  - `med_ctl_aligns.rds`
  - `med_i.rds`
  - `tel_var.rds`
- `R_intron_alignments_7/`  
  Data from alignments of all orthologous introns for two pairs of closely
  related species. Incomplete analysis.
  - `sp1_al.rds`
  - `sp2_al.rds`
- `R_intron_alignments_8/`  
  An extended set of alignments. Used in `R_intron_alignments_summary_2` below.
  - `ctl_2_aligns.rds`
  - `ctl_2_i.rds`
  - `ctl_3_aligns.rds`
- `R_intron_alignments_summary/`  
  A combined analysis of alignments performed in `R_intron_alignments` and `R_intron_alignments_5`.
  Note that these are for slightly incorrectly defined intron sets.
  - `aligns_ctl_top.rds`
  - `aligns_top.rds`
  - `bl_cov.rds`
  - `bl_ctl_cov.rds`
  - `blast/`  
    Blast of aligned sequences from the indicated sets to a set of complete genomes.
	Run using the shell scripts within the directory.
  - `extracted_seqs/`  
    Aligned intronic sequences.
    - `ctl.ctl_mammalia.fasta`
    - `ctl.ctl_teleostei.fasta`
    - `ctl_mammalia.fasta`
    - `ctl_teleostei.fasta`
    - `long.ctl_mammalia.fasta`
    - `long.ctl_teleostei.fasta`
    - `long_mammalia.fasta`
    - `long_teleostei.fasta`
    - `med_mammalia.fasta`
    - `med_teleostei.fasta`
- `R_intron_alignments_summary_2/`  
  A combined analysis of the alignments peformed in `R_intron_alignments_6`
  and `R_intron_alignments_8`.
  - `aligns_ctl_top.rds`
  - `aligns_short_2_top.rds`
  - `aligns_short_top.rds`
  - `aligns_top.rds`  
    Files with names matching the pattern `aligns_.+?_top.rds` contain
    alignment information from the top scoring local alignment for each species.
  - `bl_cov.rds`
  - `bl_cov_2.rds`
  - `bl_cov_3.rds`
  - `bl_ctl_cov.rds`  
    `bl_.+?_cov_.?.rds` files contain blast coverage data (i.e. how many
    sequences from *D. rerio* aligned to the *D. rerio* intron sequences.
  - `bl_data_3.rds`
  - `bl_files_3.rds`
  - `blast/`  
    Blast of aligned sequences from the indicated sets to a set of complete
    genomes. Data not included as excessively large.
  - `extracted_seqs/`  
    The sequences used for the above blast (extracted from the highest scoring alignments per intron
	and indicated clade.
    - `ctl.ctl_mammalia.fasta`
    - `ctl.ctl_teleostei.fasta`
    - `ctl_mammalia.fasta`
    - `ctl_teleostei.fasta`
    - `long.ctl_mammalia.fasta`
    - `long.ctl_teleostei.fasta`
    - `long_mammalia.fasta`
    - `long_teleostei.fasta`
    - `med.ctl_mammalia.fasta`
    - `med.ctl_teleostei.fasta`
    - `med_mammalia.fasta`
    - `med_teleostei.fasta`
    - `short.2_mammalia.fasta`
    - `short.2_teleostei.fasta`
    - `short_mammalia.fasta`
    - `short_teleostei.fasta`
  - `mam_var.rds`
  - `tel_var.rds`  
    Mammailan and teleost intron length variances recalculated to ensure
    internal correspondence between row and column orders of datastructures
    used in the analysis.
  - `orth.rds`
  - `orth_sc_dr.rds`
- `R_intron_orthology/`
  - `intron_orth_lscore.rds`  
    Contains information about the localised alignment score at the inferred
    intron orthologue posisitions.
- `R_trees_distances/`  
  Analysis of inferred ancestral intron sizes.
  - `all_intron_s_n.rds`  
    Numbers of introns for each species that are annotated as shorter than
    2^5^ or 2^8^ bases.
  - `class_col_3.rds`
  - `int_s_b.rds`
  - `sp_class_3.rds`
  - `sp_col_3.rds`
  - `sp_y.rds`  
    Data files giving species classification, colours and a logical taxonomic
    ordering (`sp_y.rds`).
- `R_trees_distances_2/`
  - `sp_class_3.rds`  
    More finely grained taxonomic information.
- `family_members/`  
  Data copied from `ws1` in order to faciliate extended analyses.
  - `introns_nuc_counts.txt`
  - `orthologue_transcripts/`
    - `exon_intron_stats.csv`
  - `vertebrate_family_members_1_130.txt`
