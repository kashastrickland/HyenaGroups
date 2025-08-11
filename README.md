All code and data required to reproduce results from tha manuscript "Microevolutionary consequences of social structure in wild spotted hyenas"

CODE

'preparedata.R' is used to create the data tables required for all analyses from the back-end tables of the hyena database. It requires permission to access the hyenaR R package and database, but here to show how main data were generated. 

'geneFlowDispersalHyenas_git.R' is used for all population genetic analyses of genetic differentiation, and the relationship between dispersal and distance. Also to make Figures 1 - 2, and Table 1 (and supplementary figures and tables).

'geneticgroup_animalmodels_git.R' is used for running all quantitative genetics models, including making Figures 3 and 4, and Tables 2 and 3 (and supplementary figures and tables).

DATA

'clan_dispersal.csv' ; rates of dispersal between all clans

'idDist.csv' ; geographic distance between genotyped individuals (in metres)

'life-history_ids.csv' ; fitness and meta data for individuals included in quantitative genetic analyses

'meta_geneticIDs.csv' ; meta data for genotyped individuals

'pedigree_clans.csv' ; full population pedigree

'udDists440_2.csv' ; geographic distance between clans (in metres)
