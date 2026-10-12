
##map n'gorongoro crater
## requires access to hyena database - not provided
##Kasha Strickland November 2025

library(adehabitatHR)
library(raster)
library(sf)
library(ggplot2)
library(ggspatial)
library(hyenaR)
library(lubridate)

load_package_database.full("")##requires access to hyenaR database - not provided
setwd("")

#### read in shape files
list.files("Crater arc view")
list.files("Crater shapefiles 1")
list.files("Marta shapefiles")
list.files("fixedshpfiles")
poly_crater <- st_read("Crater shapefiles 1/Crater boundary.shp")
water <- st_read("fixedshpfiles/fixedLakes.SHP")
woodland <- st_read("fixedshpfiles/fixedwoodland_notswampcropped.shp")
rivers <- st_read("fixedshpfiles/rivers_edited.shp")
swamps <- st_read("fixedshpfiles/swamps_fixed.shp")

##convert projections of everything above
st_crs(water)<-sp::CRS("+proj=longlat +datum=WGS84")
#st_crs(woodland)<-sp::CRS("+proj=longlat +datum=WGS84")
st_crs(rivers)<-sp::CRS("+proj=utm +zone=36 +south datum=WGS84")
st_crs(swamps)<-sp::CRS("+proj=longlat +datum=WGS84")

water2<-st_transform(water,crs=st_crs(poly_crater))
woodland2<-st_transform(woodland,crs=st_crs(poly_crater))
rivers2<-st_transform(rivers,crs=st_crs(poly_crater))
swamps2<-st_transform(swamps,crs=st_crs(poly_crater))

rivers2<-st_intersection(rivers2,poly_crater)
plot(rivers2$geometry)

woodland_cropped<-st_intersection(woodland2, poly_crater)
plot(woodland$geometry)
plot(woodland_cropped$geometry)

##crop from swamps
woodland_swamp_intersection<-st_intersection(woodland_cropped,swamps2)
woodland_swampcropped<-st_difference(woodland_cropped,st_union(st_geometry(
  woodland_swamp_intersection)))

plot(woodland_swampcropped$geometry)
plot(swamps2$geometry)

##create clan home ranges
clans_tracktable<-fetch_clan_track.table(find_clan_name.all(main.clans = TRUE),
                                         from = find_pop_date.observation.first(),        # set new margins based on previous insights
                                         to = find_pop_date.observation.last() - years(10))

clans_tracktable<-do.call("rbind",clans_tracktable)
str(clans_tracktable)
clans_tracktable<-clans_tracktable[clans_tracktable$sighting_clan%in%unique(clans_tracktable$clan),]
unique(clans_tracktable$sighting_clan)
hrxydata<-SpatialPointsDataFrame(clans_tracktable[,c("longitude","latitude")],clans_tracktable[,"sighting_clan"])
proj4string(hrxydata) <- CRS("+proj=longlat +datum=WGS84")  ## for example
hrxydataUTM <- spTransform(hrxydata, CRS("+proj=utm +zone=36 +south datum=WGS84"))
x <- seq(min(hrxydataUTM@coords[,1])-1000,max(hrxydataUTM@coords[,1])+1000,by=250) # where resolution is the pixel size you desire
y <- seq(min(hrxydataUTM@coords[,2])-1000,max(hrxydataUTM@coords[,2])+1000,by=250)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)
plot(xy)

udshref<-kernelUD(hrxydataUTM, h = "href", grid = xy)
udshref$T@h

clan_polys<-list()

uds440<-kernelUD(hrxydataUTM, h = 440, grid = xy)
tmpPoly<-adehabitatHR::getverticeshr(uds440,percent = 85)

clan_polys[[1]]<-adehabitatHR::getverticeshr(uds440[[1]],percent = 85)
plot(clan_polys[[1]])
clan_polys[[2]]<-adehabitatHR::getverticeshr(uds440[[2]],percent = 80)
plot(clan_polys[[2]])
clan_polys[[3]]<-adehabitatHR::getverticeshr(uds440[[3]],percent = 85)
plot(clan_polys[[3]])
clan_polys[[4]]<-adehabitatHR::getverticeshr(uds440[[4]],percent = 85)
plot(clan_polys[[4]])
clan_polys[[5]]<-adehabitatHR::getverticeshr(uds440[[5]],percent = 80)
plot(clan_polys[[5]])
clan_polys[[6]]<-adehabitatHR::getverticeshr(uds440[[6]],percent = 85)
plot(clan_polys[[6]])
clan_polys[[7]]<-adehabitatHR::getverticeshr(uds440[[7]],percent = 85)
plot(clan_polys[[7]])
clan_polys[[8]]<-adehabitatHR::getverticeshr(uds440[[8]],percent = 80)
plot(clan_polys[[8]])

all_polys<-do.call(bind,clan_polys)
all_polys$id<-names(uds440)
plot(all_polys)  

##plot to check
clan_polys_sf<-st_as_sf(all_polys)
centroids<-st_centroid(clan_polys_sf)
plot(centroids$geometry)

ggplot(poly_crater) +
  geom_sf()+
  #geom_sf(data=polyN_sf)+
  geom_sf(data=clan_polys_sf)+
  geom_sf(data=centroids)
  #geom_sf(data=polyT_sf)

#######
###plot map
##colour scheme for clan territories
oCols<-c("royalblue","gold","darkcyan","tomato2","chartreuse3","grey70")

pdf("tables_figures/Fig1A.pdf",width=12,height=8)
ggplot(poly_crater) +
  geom_sf(aes(fill = "Savannah"),alpha=0.3) +
  geom_sf(data=water2,aes(fill = "Water"),alpha=0.4)+
  geom_sf(data=woodland_swampcropped,aes(fill="Forest"),alpha=0.3)+
  geom_sf(data=swamps2,aes(fill="Swamp"),alpha=0.3)+
  geom_sf(data=rivers2,color="steelblue2",alpha=0.4)+
  #scale_color_manual(name="",guide="none",values=c("River"="steelblue2"))+
  scale_fill_manual(name="Habitat",
                    values=c("Savannah"="khaki","Water"="steelblue2","Swamp"="#CFA127","Forest"="olivedrab4"))+
  ggnewscale::new_scale_fill()+
  ###plot the territories
  geom_sf(data=clan_polys_sf,aes(col = c("Airstrip","Engitati","Forest","Lemala","Munge","Ngoitokitok","Shamba","Triangle"),
                              fill=c("Airstrip","Engitati","Forest","Lemala","Munge","Ngoitokitok","Shamba","Triangle")),alpha=0.15,linewidth=1)+
  ggspatial::annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
    style = ggspatial::north_arrow_nautical(fill = c("grey40", "white"),
                                            line_col = "grey20")
  ) +
  ggspatial::annotation_scale(
    location = "tl",
    bar_cols = c("grey60", "white")
  ) +
  scale_color_manual(values=c("Ngoitokitok"="gray20","Airstrip"=oCols[[1]],"Engitati"=oCols[[2]],"Forest"=oCols[[3]],
                              "Lemala"=oCols[[4]],"Munge"=oCols[[5]],"Shamba"=oCols[[6]],"Triangle"="mediumorchid3"),guide="none")+
  scale_fill_manual(name="Clan territory",values=c("Ngoitokitok"="gray20","Airstrip"=oCols[[1]],"Engitati"=oCols[[2]],"Forest"=oCols[[3]],
                             "Lemala"=oCols[[4]],"Munge"=oCols[[5]],"Shamba"=oCols[[6]],"Triangle"="mediumorchid3"))+
  theme_classic()+ 
  theme(legend.position = c(1.1,0.75),
        legend.title.align=0.5,
        axis.text.x = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 12, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"))+
  guides(fill = guide_legend(override.aes = list(col=c(oCols[1:5],"gray20",oCols[6],"mediumorchid3"),
                                                 fill=c(oCols[1:5],"gray20",oCols[6],"mediumorchid3"),
                                                 alpha=0.2,
                                                 linewidth=0.5)))

dev.off()
