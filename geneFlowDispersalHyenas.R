
##hyena gene flow analyses
##Kasha Strickland Novemeber 2025

library(colorspace)
library(adegenet)
library(dartR)
library(gdata)
library(dplyr)
library(data.table)
library(cowplot)

setwd("")

meta_geneticdata<-read.csv("meta_geneticdata_deID.csv")
allDisp2<-read.csv("clan_dispersal.csv")
idDist<-read.csv("idDist_deID.csv",row.names = 1)
udDist<-read.csv("udDist440_2.csv",row.names = 1)

if(!dir.exists("tables_figures")){
  dir.create("tables_figures")  
}

##https://github.com/flystar233/pyscript-bio/blob/master/glPcaFast.R
source("glPcaFast.R")

#############################
###### genomic PCA ##########
#############################

#read in filtered SNP dataset
gl<-dartR::gl.read.vcf("hyenas_deIDd.vcf")
gl2<-gl.compliance.check(gl)

popHet<-gl.report.heterozygosity(gl2)

##assigning population as individual ID's for managing the genlight object later
pop(gl2)<-indNames(gl2)
  
##removing any individuals from genetic data that does not have some metadata associated - only one individual
meta_geneticdata2<-meta_geneticdata[!is.na(meta_geneticdata$clan_birth),]
gl2<-gl.drop.ind(gl2,setdiff(indNames(gl2), meta_geneticdata2$ID))

##run the PCA - keep all axes
hPCA<-glPcaFast(gl2)

##calculate proportion of variance explained by each axis
var_frac <- hPCA$eig/sum(hPCA$eig)

##plot results of gpca and colour and cluster individuals based on their clan
oCols<-c("royalblue","gold","darkcyan","tomato2","chartreuse3","gray15",
         "grey60","mediumorchid3","darkorange")

pdf("tables_figures/Fig1C.pdf")
plot(hPCA$scores[,1],hPCA$scores[,2],type="n",cex.lab=1.4,cex.axis=1.2,
            xlab=paste0("First PC (", round(100 * var_frac[1], 1), "%)"),
            ylab=paste0("Second PC (", round(100 * var_frac[2], 1), "%)"))
s.class(hPCA$scores,fac=as.factor(meta_geneticdata2$clan_birth),add.plot=T,col=alpha(oCols,0.5))
dev.off()

pdf("tables_figures/FigS1.pdf")
plot(hPCA$scores[,3],hPCA$scores[,4],type="n",cex.lab=1.2,cex.axis=1.2,
            xlab=paste0("Third PC (", round(100 * var_frac[3], 1), "%)"),
            ylab=paste0("Fourth PC (", round(100 * var_frac[4], 1), "%)"))
s.class(hPCA$scores[,3:4],fac=as.factor(meta_geneticdata2$clan_birth),add.plot=T,col=alpha(oCols,0.5))
dev.off()

#########################
#####fsts between clans##
#########################

glC<-gl2
##reformat to code population within genlight object as clan of birth
poprecode<-data.frame(pop=pop(glC),pop2=meta_geneticdata2$clan_birth)
write.table(poprecode,"poprecode.csv",row.names=FALSE,col.names=FALSE, sep=",")
glC<-gl.recode.pop(glC,pop.recode = "poprecode.csv")
table(pop(glC))

##heterozygosity per clan
clHet<-gl.report.heterozygosity(glC)
##make table
Table1<-clHet[c("pop","Ho","HoSD","uHe","uHeSD","FIS","FISSD")]
Table1$pop<-c("Airstrip","Engitati","Forest","Lemala","Munge","Ngoitokitok","Shamba","Triangle","Non-crater")
Table1<-data.frame(Clan=Table1$pop,
                   Heterozygosity = paste0(round(Table1$Ho,3)," (",round(Table1$HoSD,3),")"),
                   Expected_heterozygosity = paste0(round(Table1$uHe,3)," (",round(Table1$uHeSD,3),")"),
                   Fis = paste0(round(Table1$FIS,3)," (",round(Table1$FISSD,3),")"))
write.csv(Table1,"tables_figures/Table1.csv",row.names = FALSE)

##calculate pairwise fsts between clans
fst<-StAMPP::stamppFst(glC)
##and format the data for later usage
fst$Fsts[upper.tri(fst$Fsts)]<-t(fst$Fsts)[upper.tri(fst$Fsts)]
fstV<-na.omit(as.data.frame(unmatrix(fst$Fsts)))
fstV$origin<-sapply(strsplit(rownames(fstV),":"), `[`, 1)
fstV$destination<-sapply(strsplit(rownames(fstV),":"), `[`, 2)

##plot fst heatmap
#reorder the matrix cols and rows
clan_order<-c("A","E","F","L","M","N","S","T","X")
fst_a<-fst$Fsts[clan_order,clan_order]
fst_a[lower.tri(fst_a)]<-NA
diag(fst_a)<-0
fst_m<-reshape2::melt(fst_a)
fst_m<-na.omit(fst_m)
fst_m$value<-as.numeric(fst_m$value)
str(fst_m)
fst_m$value2<-ifelse(fst_m$Var1==fst_m$Var2,NA,fst_m$value)

fst_m$Var2<-factor(fst_m$Var2,levels=c("X","T","S","N","M","L","F","E","A"))
fst_m$Var1<-factor(fst_m$Var1,levels=c("A","E","F","L","M","N","S","T","X"))

##plot

pdf("tables_figures/Fig1D.pdf")
ggplot(data = fst_m, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(name="Fst\n",low = "#2166AC", high = "#D6604D", mid = "#F7F7F7", 
                       midpoint = 0.02, limit = c(0,0.05))+
  geom_text(data = fst_m,aes(label = round(as.numeric(value2),3),size=5),show.legend = FALSE)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=16),legend.title.align=0.5)+
  theme(
    legend.position = c(.98, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )+
  xlab("\nClan")+
  ylab("Clan\n")+
  coord_fixed()
dev.off()

#########################################
####fst clan post-dispersal##############
#########################################

glCd<-gl2
##reformat to code population within genlight object as clan of birth
poprecodeCd<-data.frame(pop=pop(glCd),pop2=meta_geneticdata2$clan_postdisp)
write.table(poprecodeCd,"poprecodeCd.csv",row.names=FALSE,col.names=FALSE, sep=",")
glCd<-gl.recode.pop(glCd,pop.recode = "poprecodeCd.csv")
table(pop(glCd))

##calculate pairwise fsts between clans post dispersal
fstCd<-StAMPP::stamppFst(glCd)
##and format the data for later usage
fstCd$Fsts[upper.tri(fstCd$Fsts)]<-t(fstCd$Fsts)[upper.tri(fstCd$Fsts)]
fstCdV<-na.omit(as.data.frame(unmatrix(fstCd$Fsts)))
fstCdV$origin<-sapply(strsplit(rownames(fstCdV),":"), `[`, 1)
fstCdV$destination<-sapply(strsplit(rownames(fstCdV),":"), `[`, 2)

fstCd_a<-fstCd$Fsts[clan_order,clan_order]

##plot fst heatmap
##keep only upper or lower diagonal
fstCd_a[lower.tri(fstCd_a)]<-NA
diag(fstCd_a)<-0
fstCd_m<-reshape2::melt(fstCd_a)
fstCd_m<-na.omit(fstCd_m)
fstCd_m$value<-as.numeric(fstCd_m$value)
str(fstCd_m)
fstCd_m$value2<-ifelse(fstCd_m$Var1==fstCd_m$Var2,NA,fstCd_m$value)
##plot heatmap
fstCd_m$Var2<-factor(fstCd_m$Var2,levels=c("X","T","S","N","M","L","F","E","A"))
fstCd_m$Var1<-factor(fstCd_m$Var1,levels=c("A","E","F","L","M","N","S","T","X"))

##plot
pdf("tables_figures/FigS2.pdf")
ggplot(data = fstCd_m, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(name="Fst\n",low = "#2166AC", high = "#D6604D", mid = "#F7F7F7", 
                       midpoint = 0.02, limit = c(0,0.05))+
  geom_text(data = fstCd_m,aes(label = round(as.numeric(value2),3),size=5),show.legend = FALSE)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=16),legend.title.align=0.5)+
  theme(
    legend.position = c(.98, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )+
  xlab("\nClan")+
  ylab("Clan\n")+
  coord_fixed()
dev.off()

#################################################
####test relationship between spatial distance, dispersal and fst####
#################################################

####plot pairwise dispersal as heatmap
##need to reorganise to keep zeros and have the rate of philopatry on the diagonal
allDisp2$destination<-factor(allDisp2$destination,levels=c("X","T","S","N","M","L","F","E","A"))
allDisp2$origin<-factor(allDisp2$origin,levels=c("A","E","F","L","M","N","S","T","X"))

pdf("tables_figures/Figure1B.pdf",width=8,height=8)
ggplot(data = allDisp2, aes(x=origin, y=destination, fill=probDisp)) + 
  geom_tile()+
  scale_fill_gradient2(low = "#2166AC", high = "#D6604D", 
                       midpoint = 0.12, limit = c(0,0.6))+
  geom_text(data = allDisp2,aes(label = round(as.numeric(probDisp),3),size=5),show.legend = FALSE)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=18),
        legend.title=element_blank())+
  xlab("\nOrigin")+
  ylab("Destination\n")+
coord_fixed()
dev.off()

##merge with pairwise Fst 
fst_both<-merge(fstV,fstCdV,by=c("origin","destination"),all.x=TRUE)
colnames(fst_both)[3:4]<-c("observedfst_birthClan","observedfst_dispClan")
disp_fst<-merge(fst_both,allDisp2,by=c("origin","destination"),all.x=TRUE)

disp_fst[is.na(disp_fst)]<-0

##test relationship between spatial distance, genetic distance and dispersal probability

str(udDist)##check the matrix of spatial distances

##make matrices of fst and probabilty of dispersal

##blank data frame to capture all and organise matrices 
##removing X clan - which is non-crater individuals
tdf<-data.frame(origin=unique(disp_fst$origin),
                destination=unique(disp_fst$origin),
                observedfst_birthClan=NA,
                observedfst_dispClan=NA,
                totalN=NA,
                total=NA,
                probDisp=NA)

disp_fstM<-rbind(disp_fst,tdf)
disp_fstM<-disp_fstM[order(disp_fstM$origin),]
disp_fstM<-disp_fstM[order(disp_fstM$destination),]
dispM<-matrix(disp_fstM$probDisp,ncol=9,nrow=9)
colnames(dispM)<-rownames(dispM)<-unique(disp_fstM$origin)
dispM<-dispM[rownames(dispM)!="X",colnames(dispM)!="X"]

fst_M<-matrix(disp_fstM$observedfst_birthClan,ncol=9,nrow=9)
colnames(fst_M)<-rownames(fst_M)<-unique(disp_fstM$origin)
fst_M<-fst_M[rownames(fst_M)!="X",colnames(fst_M)!="X"]

fstCd_M<-matrix(disp_fstM$observedfst_dispClan,ncol=9,nrow=9)
colnames(fstCd_M)<-rownames(fstCd_M)<-unique(disp_fstM$origin)
fstCd_M<-fstCd_M[rownames(fstCd_M)!="X",colnames(fstCd_M)!="X"]

udDist<-as.matrix(udDist)
diag(udDist)<-NA

asnipe::mrqap.dsp(dispM~udDist, test.statistic = "beta")
asnipe::mrqap.dsp(fst_M~dispM+udDist, test.statistic = "beta")
asnipe::mrqap.dsp(fst_M~udDist, test.statistic = "beta")
asnipe::mrqap.dsp(fstCd_M~dispM+udDist, test.statistic = "beta")

df_fst_plots<-data.frame(disp=as.vector(dispM),fst=as.vector(fst_M),dist=as.vector(udDist))

pdf("tables_figures/FigureS4.pdf")
ggplot(df_fst_plots,aes(x=dist,y=disp))+
  geom_point()+
  theme_classic()+
  xlab("\nGeographic distance between clans")+
  ylab("Probability of dispersal between clans\n")+
  theme(text = element_text(size=14),
        plot.margin = unit(c(1,2,1,2), "lines"))
dev.off()

fig2B<-ggplot(df_fst_plots,aes(x=disp,y=fst))+
  geom_point()+
  geom_smooth(method = "lm", se = F,colour="black",alpha=1/2) +
  geom_line(aes(y = stage(start = fst, after_stat = ymin)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  geom_line(aes(y = stage(fst, ymax)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  #geom_line(stat="smooth",method="glm",se=FALSE)+
  theme_classic()+
  xlab("\nProbability of dispersal between clans")+
  ylab("Fst between clans\n")+
  theme(text = element_text(size=14),
        plot.margin = unit(c(1,2,1,2), "lines"))+
  ylim(0.01,0.045)

fig2C<-ggplot(df_fst_plots,aes(x=dist,y=fst))+
  geom_point()+
  geom_smooth(method = "lm", se = F,colour="black",alpha=1/2) +
  geom_line(aes(y = stage(start = fst, after_stat = ymin)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  geom_line(aes(y = stage(fst, ymax)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  #geom_line(stat="smooth",method="glm",se=FALSE)+
  theme_classic()+
  xlab("\nGeographic distance between clans (m)")+
  ylab("Fst between clans\n")+
  theme(text = element_text(size=14),
        plot.margin = unit(c(1,2,1,2), "lines"))+
  ylim(0.01,0.045)

df_fstCd_plots<-data.frame(disp=as.vector(dispM),fst=as.vector(fstCd_M),dist=as.vector(udDist))

pdf("tables_figures/FigureS3.pdf")
ggplot(df_fstCd_plots,aes(x=disp,y=fst))+
  geom_point()+
  geom_smooth(method = "lm", se = F,colour="black",alpha=1/2) +
  geom_line(aes(y = stage(start = fst, after_stat = ymin)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  geom_line(aes(y = stage(fst, ymax)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  #geom_line(stat="smooth",method="glm",se=FALSE)+
  theme_classic()+
  xlab("\nProbability of dispersal between clans")+
  ylab("Fst between clans\n")+
  theme(text = element_text(size=14),
        plot.margin = unit(c(1,2,1,2), "lines"))+
  ylim(0.01,0.045)
dev.off()

###############################################################################
##null model to test the probability of observed differentiation given random##
##dispersal of males - females retain natal clan membership####################

glNd<-glC
fst_randD<-list()
for (j in 1:1000){
  glTmp2<-glNd
  popRd<-meta_geneticdata2[c("ID","sex","clan_postdisp")]
  popRd$sex[is.na(popRd$sex)]<-"unknown"
  clans <- unique(popRd$clan_postdisp)
  popRd$clan_postdisp[popRd$sex=="male"]<-sample(clans, sum(popRd$sex=="male"), replace=TRUE)
  glTmp2$pop<-as.factor(popRd$clan_postdisp)
  fstTmp2<-StAMPP::stamppFst(glTmp2)
  fstTmp2$Fsts[upper.tri(fstTmp2$Fsts)]<-t(fstTmp2$Fsts)[upper.tri(fstTmp2$Fsts)]
  fstTmp2V<-na.omit(as.data.frame(unmatrix(fstTmp2$Fsts)))
  fstTmp2V$origin<-sapply(strsplit(rownames(fstTmp2V),":"), `[`, 1)
  fstTmp2V$destination<-sapply(strsplit(rownames(fstTmp2V),":"), `[`, 2)
  fst_randD[[j]]<-fstTmp2V
  cat(j, "\r")
}

##reordering dataframes generated in above code
fst_randDO<-list()
or<-unique(fstV$destination)
for (w in 1:length(fst_randD)){
  df<-fst_randD[[w]]
  df1 <- df %>% 
    arrange(factor(origin, levels = or))
  df2 <- df1 %>% 
    arrange(factor(destination, levels = or))
  fst_randDO[[w]]<-df2[,1]
  cat(w)
}

fst_randDOdf<-data.frame(df2[,2:3],do.call("cbind",fst_randDO))
names(fstV)[1]<-"observed"
fst_randDOdf<-merge(fstV,fst_randDOdf,by = c("origin","destination"))

nullD_qs<-data.frame(fst_randDOdf[,1:3],
                     as.data.frame(t(apply(fst_randDOdf[,4:ncol(fst_randDOdf)],1,
                                           quantile,probs=c(0.05,0.95)))))

names(nullD_qs)
str(nullD_qs)
nullDiagdf<-data.frame(origin=unique(nullD_qs$origin),
                destination=unique(nullD_qs$origin),
                observed=NA,
                X5.=NA,
                X95.=NA)

colnames(nullDiagdf)<-colnames(nullD_qs)
nullD_qsM<-rbind(nullD_qs,nullDiagdf)
nullD_qsM<-nullD_qsM[order(nullD_qsM$origin),]
nullD_qsM<-nullD_qsM[order(nullD_qsM$destination),]

mat95<-matrix(nullD_qsM$X95.,
  length(unique(c(nullD_qsM$origin,nullD_qsM$destination))),
  length(unique(c(nullD_qsM$origin,nullD_qsM$destination))))

rownames(mat95)<-colnames(mat95)<-unique(unique(c(nullD_qsM$origin,nullD_qsM$destination)))

fig2A<-ggplot(nullD_qs,aes(x=X95.,y=observed))+
     geom_point()+
     xlab("\nExpected Fst (upper 95% quantile)")+
     ylab("Fst between clans\n")+
     theme_classic()+
     theme(text = element_text(size=14),
           plot.margin = unit(c(1,2,1,2), "lines"))+
     xlim(0.005,0.045)+
     ylim(0.01,0.045)

###############################
######IBD analysis#############
###############################

##genetic distance matrix [nei distance] -- takes a while to run
gen<-StAMPP::stamppNeisD(gl2,FALSE)## takes a while to run

# Ensure symmetric distance matrices by matching colnames to rownames
colnames(gen)<-rownames(gen)
idDist<-as.matrix(idDist)
colnames(idDist)<-rownames(idDist)
gen2<-gen[rownames(gen)%in%rownames(idDist),colnames(gen)%in%colnames(idDist)]
idDist2<-idDist[rownames(gen2),colnames(gen2)]

##check order
identical(rownames(gen2), rownames(idDist2))

##run isolation-by-distance analyses
diag(gen2)<-NA
diag(idDist2)<-NA
vegan::mantel(gen2,idDist2)

##plot isolation-by-distance
pedH<-read.csv(file="pedH_deID.csv")
Amat <- as.matrix(nadiv::makeA(pedH[, 1:3]))
AmtV<-data.frame(unmatrix(Amat))
AmtV$pair<-rownames(AmtV)

ibd_df<-data.frame(dist=as.vector(idDist2),
                   gen=as.vector(gen2),
                   pair=rownames(data.frame(unmatrix(gen2))))

ibd_dfA<-left_join(AmtV[,c(2,1)],ibd_df)

fig2D<-ggplot(ibd_df,aes(x=dist,y=gen))+
  geom_point(col="darkgrey",alpha=1/7)+
  geom_smooth(method = "lm", se = F,colour="black",alpha=1/2) +
  geom_line(aes(y = stage(start = gen, after_stat = ymin)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  geom_line(aes(y = stage(gen, ymax)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  theme_classic()+
  xlab("\nGeographic distance between individuals (m)")+
  ylab("Genetic distance between individuals (Nei's D)\n")+
  theme(text = element_text(size=14),
        plot.margin = unit(c(1,2,1,2), "lines"))
  
pdf("tables_figures/Figure2.pdf",height=10,width=11)
plot_grid(fig2A,fig2B,fig2C,fig2D,ncol=2,labels = c("A","B","C","D"))
dev.off()

save.image("")
