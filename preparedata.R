
##requires permission
drat::addRepo("hyenaproject")
install.packages("hyenaR")

####Generate meta data for sequenced individuals, dispersal information between clans, full pedigree and life history data for genetic group animal models

####Kasha Strickland
####13th March 2025

library(tidyverse)
library(hyenaR)
library(dartR)
library(data.table)
library(adehabitatHR)
library(sf)

##load database - requires permission 
load_package_database.full(db.path="~/Documents/hyenaR_files/hyena140325.sqlite")

########
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Hyenas/group_geneticdifferences_ms")

##read in genetic data
gl<-dartR::gl.read.vcf("data/Hyenas_2ndRun_mincov10MaxCov110menderr01maf1miss70_dupsremovalidswap.recode.vcf")

##generate metadata for genotyped individuals
meta<-data.frame(ID=indNames(gl))
starting_table <- create_id_starting.table(from = find_pop_date.observation.first(),
                                         to = find_pop_date.observation.last())
starting_table |>
  mutate(sex = fetch_id_sex(ID),
         clan_birth = fetch_id_clan.birth(ID),
         birth_date = fetch_id_date.birth(ID),
         birth_year = year(birth_date),
         ) -> all_ids

meta_geneticdata<-left_join(meta,all_ids)

##rename all non-crater clans as "X"
meta_geneticdata[meta_geneticdata$clan_birth%in%c("C","X","R","B","U","Y"),]$clan_birth<-"X"

#########################################################
#####calculate pairwise dispersal between clans##########
#########################################################

##code from Alex Courtiol - 13/05/2024
##extract all dispersal events
create_id_life.transition.table() |> 
  mutate(origin = recode(origin,"C"="X","R"="X","B"="X","U"="X")) |>
  mutate(destination = recode(destination,"C"="X","R"="X","B"="X","U"="X")) |>
  filter(origin %in% c("X",find_clan_name.all()),
         destination %in% c("X",find_clan_name.all())) |> 
  mutate(year = lubridate::year(date),
         sex = fetch_id_sex(ID)) -> dispersingIDs

dispersingIDs %>% 
  group_by(ID) %>% 
  slice(tail(row_number(),1)) -> dispersingIDs_last

##check it's worked - find IDs that have more than observation and confirm last is taken
table(dispersingIDs$ID)
dispersingIDs[dispersingIDs$ID=="N-144",]
dispersingIDs_last[dispersingIDs_last$ID=="N-144",]
##looks all good

##post-dispersing clan membership
dispersingIDs_dest<-dispersingIDs_last[c("ID","destination")]

##merge with genetic data 
meta_geneticdata2<-left_join(meta_geneticdata,dispersingIDs_dest)

##replace with clan of birth if individual did not disperse
meta_geneticdata2[is.na(meta_geneticdata2$destination),]$destination<-meta_geneticdata2[is.na(meta_geneticdata2$destination),]$clan_birth
colnames(meta_geneticdata2)[6]<-"clan_postdisp"

##calculate number of dispersal events per year from one clan (origin) to another (destination)
create_id_life.transition.table() |> 
  mutate(origin = recode(origin,"C"="X","R"="X","B"="X","U"="X")) |>
  mutate(destination = recode(destination,"C"="X","R"="X","B"="X","U"="X")) |>
  filter(origin %in% c("X",find_clan_name.all()),
         destination %in% c("X",find_clan_name.all())) |> 
  mutate(year = lubridate::year(date),
         sex = fetch_id_sex(ID)) |> 
  filter(sex == "male") |> 
  summarise(events = n(), .by = c("origin", "destination", "year")) -> dispProbs

## total number of dispersal events per year per clan
dispProbs |> 
  group_by(year,origin) |>
  summarise(total = sum(events)) -> sumDisp

dispProbs2<-merge(dispProbs,sumDisp,by=c("origin","year"))
dispProbs2$probDisp<-dispProbs2$events/dispProbs2$total

##now use the above to summarise across whole time period - i.e., sum up across years
dispProbs2 |>
  group_by(origin,destination) |>
  summarise(totalN = sum(events)) -> dispProbsTotal

dispProbsTotal |> 
  group_by(origin) |>
  summarise(total = sum(totalN)) -> allDisp

allDisp2<-merge(dispProbsTotal,allDisp,by=c("origin"))
allDisp2$probDisp<-allDisp2$totalN/allDisp2$total

#####################################################################
#######create clan home ranges to calculate territory centroids######
#####################################################################

##using same time frame for clan territories as used in quantitative genetics analyses as below
clans_tracktable<-fetch_clan_track.table(find_clan_name.all(main.clans = TRUE),
                                         from = find_pop_date.observation.first(),
                                         to = find_pop_date.observation.last() - years(10))

clans_tracktable<-do.call("rbind",clans_tracktable)
str(clans_tracktable)

##keep main crater clans only
clans_tracktable<-clans_tracktable[clans_tracktable$sighting_clan%in%unique(clans_tracktable$clan),]
unique(clans_tracktable$sighting_clan)

##format spatial data for calculation of utilisation distributions 
hrxydata<-SpatialPointsDataFrame(clans_tracktable[,c(6,5)],clans_tracktable[,"sighting_clan"])
proj4string(hrxydata) <- CRS("+proj=longlat +datum=WGS84")  ## for example
hrxydataUTM <- spTransform(hrxydata, CRS("+proj=utm +zone=36 ellps=WGS84"))
x <- seq(min(hrxydataUTM@coords[,1])-1000,max(hrxydataUTM@coords[,1])+1000,by=250) # where resolution is the pixel size you desire
y <- seq(min(hrxydataUTM@coords[,2])-1000,max(hrxydataUTM@coords[,2])+1000,by=250)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)
plot(xy)

unique(clans_tracktable$sighting_clan)
hrxydata<-SpatialPointsDataFrame(clans_tracktable[,c("longitude", "latitude")],clans_tracktable[,"sighting_clan"])
proj4string(hrxydata) <- CRS("+proj=longlat +datum=WGS84")  ## for example
hrxydata_UTM <- spTransform(hrxydata, CRS("+proj=utm +zone=36 ellps=WGS84"))
uds440<-kernelUD(hrxydata_UTM[,1] , h = 440, grid=xy)
##h = 440 identified previously via use of h = href, and then standardizing to max h value

##calculate 88% polygon of UDs, centroid of that polygon and then geographic distance between centroids
clan_polys440<-adehabitatHR::getverticeshr(uds440,percent = 88)
poly440_sf<-st_as_sf(clan_polys440)
cent_440<-sf::st_centroid(poly440_sf)
udDist440<-sf::st_distance(cent_440)

udDist440_2<-matrix(udDist440[1:8,],8,8)
rownames(udDist440_2)<-colnames(udDist440_2)<-cent_440$id

###estimate individuals' centroids

##extract sightings tables to calculate number of sightings per indivudual
##keep IDs that have at least 5 sightigs
starting_table <- create_id_starting.table(from = find_pop_date.observation.first(),        # set new margins based on previous insights
                                           to = find_pop_date.observation.last(),
                                           lifestage = "cub",
                                           lifestage.overlap = "start")

genIDswithinStudy<-starting_table[starting_table$ID%in%indNames(gl),]

sighID<-fetch_id_sighting.table(genIDswithinStudy$ID,
                                from = find_pop_date.observation.first(),
                                to = find_pop_date.observation.last())

NsighID2<-lapply(sighID,nrow)
sighID2<-NsighID2[NsighID2>4]

##draw individual home ranges as MCPs using hyenaR functions
idHR<-fetch_id_homerange(names(sighID2))
idHR
idHRcent<-sf::st_centroid(idHR)
idDist<-sf::st_distance(idHRcent)

colnames(idDist)<-rownames(idDist)<-names(sighID2)
dim(idDist)

##save datasets for downstream analysis
write.csv(meta_geneticdata2,file="data/meta_geneticIDs.csv")
write.csv(allDisp2,file="data/clan_dispersal.csv")
write.csv(udDist440_2,file="data/udDist440_2.csv")
write.csv(idDist,file="data/idDist.csv")

#################################################
#### extract life-history data for QG models ####
#################################################

##find the full pedigree
pedH<-create_pop_pedigree()
pedH<-data.frame(pedH)
pedH$birth_year<-year(fetch_id_date.birth(pedH$ID))

##append clan membership to individuals for use in gg models
pedH$birth_clan<-fetch_id_clan.birth(pedH$ID)
pedH<-left_join(pedH,dispersingIDs_dest)
pedH[is.na(pedH$destination),]$destination<-pedH[is.na(pedH$destination),]$birth_clan
colnames(pedH)[ncol(pedH)]<-"clan_postdisp"

##replace and group all non-crater clan names to be X
pedH[pedH$birth_clan%in%c("B","C","R","U"),]$birth_clan<-"X"
pedH[pedH$clan_postdisp%in%c("B","C","R","U"),]$clan_postdisp<-"X"
table(pedH$birth_clan)
table(pedH$clan_postdisp)

###########
###extract life history and phenotypic data for individuals that will be included in analyses
starting_table <- create_id_starting.table(from = find_pop_date.observation.first(),
                                           to = find_pop_date.observation.last() - years(10),
                                           lifestage = "cub",
                                           lifestage.overlap = "start")

starting_table |>
  mutate(mother = fetch_id_id.mother.genetic(ID),
         sex = fetch_id_sex(ID),
         birth_clan = fetch_id_clan.birth(ID),
         birth_date = fetch_id_date.birth(ID),
         birth_year = year(birth_date),
         birth_rank = fetch_id_rank.birth(ID),
         alive_at_end = fetch_id_is.alive(ID, at = find_pop_date.observation.last()),
         lifespan = fetch_id_duration.lifespan(ID),
         LRS = fetch_id_number.offspring(ID)) -> life_history_ids

##add on clan post dispersal
life_history_ids<-left_join(life_history_ids,dispersingIDs_dest)
life_history_ids[is.na(life_history_ids$destination),]$destination<-life_history_ids[is.na(life_history_ids$destination),]$birth_clan
colnames(life_history_ids)[ncol(life_history_ids)]<-"clan_postdisp"

##calculate clan size at birth
life_history_ids |>
  mutate(clan_size = fetch_clan_number(clan = birth_clan, at = birth_date)) -> life_history_ids

##rename all non-crater clans to be "X"
life_history_ids[life_history_ids$birth_clan%in%c("B","C","R","U"),]$birth_clan<-"X"
life_history_ids[life_history_ids$clan_postdisp%in%c("B","C","R","U"),]$clan_postdisp<-"X"

table(life_history_ids$birth_clan)
table(life_history_ids$clan_postdisp)

##remove individuals still alive at last observation date or with missing key info
life_history_ids<-life_history_ids[life_history_ids$alive_at_end==FALSE&
                                     !is.na(life_history_ids$mother)&
                                     !is.na(life_history_ids$birth_rank)&
                                     !is.na(life_history_ids$birth_year),]

table(life_history_ids$birth_clan)
table(life_history_ids$clan_postdisp)

##write out data for downstream analyses
write.csv(life_history_ids,file="data/life_history_ids.csv",row.names = TRUE)
write.csv(pedH,file="data/pedigree_clans.csv",row.names = TRUE)
