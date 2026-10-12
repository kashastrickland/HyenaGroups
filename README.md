All code and data required to reproduce results from the manuscript:
"Linking social structure to fine-scale genetic structure and genetic variance in fitness in spotted hyenas". 

Note that individual hyena names were given new ID's for all data that are not the same as the original name in the field database for the project, but these ID's are the same between datasets in this repository.

CODE

"Figure 1.R" reproduces the map in Figure 1.

'geneFlowDispersalHyenas.R' is used for all population genetic analyses of genetic differentiation, and the relationship between dispersal and distance. Also to make Figures 1 - 2, and Table 1 (and supplementary figures and tables).

'geneticgroup_animalmodels.R' is used for running all quantitative genetics models, including making Figures 3 and 4, and Tables 2 and 3 (and supplementary figures and tables).

in path /data

'clan_dispersal.csv' ; rates of dispersal between all clans

'idDist_deID.csv' ; geographic distance between genotyped individuals (in metres)

'life-history_ids_deID.csv' ; fitness and meta data for individuals included in quantitative genetic analyses

'meta_geneticdata_deID.csv' ; meta data for genotyped individuals

'pedH_deID.csv' ; full population pedigree

'udDists440_2.csv' ; geographic distance between clans (in metres)

in path /shapefiles

all shape files required to reproduce Figure 1. 

Please note the processed SNP genotype file "hyenas_deID.vcf" can be found at 10.5281/zenodo.17609304. This file is required to reproduce genetic analyses, and is the genotype file post-bioinformatics and QC and filtering as described in the manuscript.
