# Determinatio of taxon sizes

Taxon sizes were determined by descending taxonomy trees downloaded
from the catalogue of life project
(<https://www.catalogueoflife.org/data/metadata>) using the
`genomes/catalogue_of_life/perl_scripts/species_counts.pl`
script. Since the catalogue of life does not recognise the taxa
Vertebrata or Teleostia, the sizes of the immediate children of the
Chordata and Actinopterygii classes were determined using the above
script.  The classes were then defined as being members of Vertebrata
or Teleostia by inspecting the NCBI taxonomy tree which includes these
divisions.

All scripts were called from shell scripts giving the command arguments
and the output file names.
