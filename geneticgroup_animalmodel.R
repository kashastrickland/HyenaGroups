

##Hyena genetic group animal models
##Kasha Strickland November 2025

library(RColorBrewer)
require(hyenaR)
require(dplyr)
require(tidyr)
library(ggplot2)
require(MCMCglmm)
require(brms)
require(nadiv)
library(pedigreemm) 
library(data.table)
library(tidybayes)
library(wesanderson)
library(pedigree)
library(reshape2)

setwd("")

if(!dir.exists("tables_figures")){
  dir.create("tables_figures")  
}

if(!dir.exists("modelResults_Robjects")){
  dir.create("modelResults_Robjects")  
}

##read data
life_history_ids<-read.csv(file="life_history_ids_deID.csv")
pedH<-read.csv(file="pedH_deID.csv")

##unknown sex - randomly assign a sex (removing selectively removes those that die young)
life_history_ids[is.na(life_history_ids$sex),]$sex<-sample(c("male","female"),nrow(life_history_ids[is.na(life_history_ids$sex),]),
                                                           replace = TRUE)

str(life_history_ids)

######################################
##original inverse A matrix
ainvOut <- inverseA(pedH[, 1:3])
Ainv <- ainvOut$Ainv
pedH$f <- ainvOut$inbreeding
Ainv  <- as(Ainv, "dgCMatrix") 
Ainv@Dimnames[[1L]] <- as.character(pedH$ID)

##some formatting to replace NA's in parentage with clan of birth
naSire <- which(is.na(pedH[, 2]))
naDam  <- which(is.na(pedH[, 3]))

pedH$damOriginal <- pedH$dam
pedH$sireOriginal <- pedH$sire

pedHc<-pedH

str(life_history_ids)
str(pedHc)

####check sample sizes per clan

table(life_history_ids$birth_clan)

##remove IDs not born into a crater clan from life history data
##sampling means estimates of fitness will be poor

life_history_ids_crater<-life_history_ids[!life_history_ids$birth_clan=="X",]
table(life_history_ids_crater$birth_clan)
str(life_history_ids_crater)

#### keeping "immigrant" ancestry in pedigree for all the following calculations
#### to account for individuals true ancestry 

##replace NA in parentage with clan
pedHc[naSire, 2] <- pedHc[naSire, "birth_clan"]
pedHc[naDam, 3] <- pedHc[naDam, "birth_clan"]

##
Qc <- ggcontrib(pedHc[, 1:3], ggroups = unique(pedHc$birth_clan))
pedHc <- cbind(pedHc, Qc)

##
pedc <- pedigree(sire=pedHc[,"sireOriginal"],dam=pedHc[,"damOriginal"],label=pedHc[,"ID"])
TInvc <- as(pedc, "sparseMatrix") 
DDc <- Dmat(pedc)
DInvc <- Diagonal(x=1/DDc)
invA_testc<- t(TInvc) %*% DInvc %*% (TInvc)
sum(invA_testc - Ainv)

#########################################
##now create scaled A matrices per clan##
#########################################

##X - immigrant/non-crater
scalingXc <- ifelse(pedHc$X>0,1/(pedHc$X)^2,1e12)
DD0cX <- 1- pedHc$X *(1-DDc)
DInvTildecX <- Diagonal(x=1/(DD0cX) * scalingXc) 
AinvcX <- crossprod(sqrt(DInvTildecX) %*% TInvc)

##Airstrip clan - A
scalingAc <- ifelse(pedHc$A>0,1/(pedHc$A)^2,1e12)
DD0cA <- 1- pedHc$A *(1-DDc)
DInvTildecA <- Diagonal(x=1/(DD0cA) * scalingAc) 
AinvcA <- crossprod(sqrt(DInvTildecA) %*% TInvc)

##Munge clan - M
scalingMc <- ifelse(pedHc$M>0,1/(pedHc$M)^2,1e12)
DD0cM <- 1- pedHc$M *(1-DDc)
DInvTildecM <- Diagonal(x=1/(DD0cM) * scalingMc) 
AinvcM <- crossprod(sqrt(DInvTildecM) %*% TInvc)

##Lemala clan - L
scalingLc <- ifelse(pedHc$L>0,1/(pedHc$L)^2,1e12)
DD0cL <- 1- pedHc$L *(1-DDc)
DInvTildecL <- Diagonal(x=1/(DD0cL) * scalingLc) 
AinvcL <- crossprod(sqrt(DInvTildecL) %*% TInvc)

##Engitati clan - E
scalingEc <- ifelse(pedHc$E>0,1/(pedHc$E)^2,1e12)
DD0cE <- 1- pedHc$E *(1-DDc)
DInvTildecE <- Diagonal(x=1/(DD0cE) * scalingEc) 
AinvcE <- crossprod(sqrt(DInvTildecE) %*% TInvc)

##Ngoitokitok clan - N
scalingNc <- ifelse(pedHc$N>0,1/(pedHc$N)^2,1e12)
DD0cN <- 1- pedHc$N *(1-DDc)
DInvTildecN <- Diagonal(x=1/(DD0cN) * scalingNc) 
AinvcN <- crossprod(sqrt(DInvTildecN) %*% TInvc)

##Forest clan - F
scalingFc <- ifelse(pedHc$F>0,1/(pedHc$F)^2,1e12)
DD0cF <- 1- pedHc$F *(1-DDc)
DInvTildecF <- Diagonal(x=1/(DD0cF) * scalingFc) 
AinvcF <- crossprod(sqrt(DInvTildecF) %*% TInvc)

## Triangle clan - T
scalingTc <- ifelse(pedHc$T>0,1/(pedHc$T)^2,1e12)
DD0cT <- 1- pedHc$T *(1-DDc)
DInvTildecT <- Diagonal(x=1/(DD0cT) * scalingTc) 
AinvcT <- crossprod(sqrt(DInvTildecT) %*% TInvc)

## Shamba clan - S
scalingSc <- ifelse(pedHc$S>0,1/(pedHc$S)^2,1e12)
DD0cS <- 1- pedHc$S *(1-DDc)
DInvTildecS <- Diagonal(x=1/(DD0cS) * scalingSc) 
AinvcS <- crossprod(sqrt(DInvTildecS) %*% TInvc)

##format matrices correctly
AinvcX <- as(AinvcX, "dgCMatrix") 
AinvcX@Dimnames[[1L]] <- as.character(pedHc$ID)

AinvcA <- as(AinvcA, "dgCMatrix") 
AinvcA@Dimnames[[1L]] <- as.character(pedHc$ID)

AinvcL <- as(AinvcL, "dgCMatrix") 
AinvcL@Dimnames[[1L]] <- as.character(pedHc$ID)

AinvcM <- as(AinvcM, "dgCMatrix") 
AinvcM@Dimnames[[1L]] <- as.character(pedHc$ID)

AinvcE <- as(AinvcE, "dgCMatrix") 
AinvcE@Dimnames[[1L]] <- as.character(pedHc$ID)

AinvcN <- as(AinvcN, "dgCMatrix") 
AinvcN@Dimnames[[1L]] <- as.character(pedHc$ID)

AinvcF <- as(AinvcF, "dgCMatrix") 
AinvcF@Dimnames[[1L]] <- as.character(pedHc$ID)

AinvcT <- as(AinvcT, "dgCMatrix") 
AinvcT@Dimnames[[1L]] <- as.character(pedHc$ID)

AinvcS <- as(AinvcS, "dgCMatrix") 
AinvcS@Dimnames[[1L]] <- as.character(pedHc$ID)

##organise and format the life history data for models
##create the clan indices for the model below

pedHc$IndexX <- pedHc$IndexA <- 
  pedHc$IndexL <- 
  pedHc$IndexM <- 
  pedHc$IndexE <- 
  pedHc$IndexN <- 
  pedHc$IndexF <- 
  pedHc$IndexT <- pedHc$IndexS <- pedHc$ID

pedHc$IndexX <- ifelse(pedHc$X>0, pedHc$IndexX, NA)
pedHc$IndexA <- ifelse(pedHc$A>0, pedHc$IndexA, NA)
pedHc$IndexL <- ifelse(pedHc$L>0, pedHc$IndexL, NA)
pedHc$IndexM <- ifelse(pedHc$M>0, pedHc$IndexM, NA)
pedHc$IndexE <- ifelse(pedHc$E>0, pedHc$IndexE, NA)
pedHc$IndexN <- ifelse(pedHc$N>0, pedHc$IndexN, NA)
pedHc$IndexF <- ifelse(pedHc$F>0, pedHc$IndexF, NA)
pedHc$IndexT <- ifelse(pedHc$T>0, pedHc$IndexT, NA)
pedHc$IndexS <- ifelse(pedHc$S>0, pedHc$IndexS, NA)

lh_ids_crater<-left_join(life_history_ids_crater,pedHc[c(1,7,10:27)])
str(lh_ids_crater)

#####################################################
###relative fitness with a gaussian distribution#####
#####################################################

lh_ids_crater$lrs_r<-lh_ids_crater$LRS/mean(lh_ids_crater$LRS)
lh_ids_crater$lifespan_r<-lh_ids_crater$lifespan/mean(lh_ids_crater$lifespan)

########################################################
##############run models with homogeneous variances#####
########################################################

priors_c <- list(R = list(V = diag(1), nu = 0.001),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000),
                        G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000),
                        G3 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000),
                        G4 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000)))

GGrelLRS_c <- MCMCglmm(lrs_r ~ sex*birth_rank +  f + clan_size + A + E + F + L + M + N + T + S, 
                       random = ~ ID + mother + birth_year + birth_clan,
                       ginverse = list(ID=Ainv), 
                       data = lh_ids_crater,
                       pr=TRUE,
                       prior = priors_c,
                       nitt = 21000, burn = 1000, thin = 10)
summary(GGrelLRS_c)
plot(GGrelLRS_c)

GGrelLong_c <- MCMCglmm(lifespan_r ~ sex*birth_rank +  f + clan_size + A + E + F + L + M + N + T + S, 
                       random = ~ ID + mother + birth_year + birth_clan,
                       ginverse = list(ID=Ainv), 
                       data = lh_ids_crater,
                       pr=TRUE,
                       prior = priors_c,
                       nitt = 21000, burn = 1000, thin = 10)
summary(GGrelLong_c)
plot(GGrelLong_c)

#############################
####make results tables######
#############################

# Extracting fixed effects
fixed_effects_lrs <- round(summary(GGrelLRS_c)$solutions,3)
fixed_effects_lrs_df <- as.data.frame(fixed_effects_lrs)
##reorder
fixed_effects_lrs_df<-fixed_effects_lrs_df [c(rownames(fixed_effects_lrs_df[1:11,]),"S","T","sexmale:birth_rank"),]

# Extracting random effects
random_effects_lrs<-as.data.frame(GGrelLRS_c$VCV)
random_effects_lrs_df <- data.frame(mean=round(colMeans(random_effects_lrs),3),
                                        lwr5=apply(random_effects_lrs[1:5],2,quantile,0.05),
                                        upr95=apply(random_effects_lrs[1:5],2,quantile,0.95))

colnames(random_effects_lrs_df)<-colnames(fixed_effects_lrs_df[1:3])
random_effects_lrs_df[,2:3]<-round(random_effects_lrs_df[,2:3],3)

##combine fixed and random effects coefficients into one dataframe
rs_lrs<-rbind(fixed_effects_lrs_df[1:3],random_effects_lrs_df)

###lifespan
# Extracting fixed effects
fixed_effects_long <- round(summary(GGrelLong_c)$solutions,3)
fixed_effects_long_df <- as.data.frame(fixed_effects_long)
##reorder
fixed_effects_long_df<-fixed_effects_long_df [c(rownames(fixed_effects_long_df[1:11,]),"S","T","sexmale:birth_rank"),]
# Extracting random effects
random_effects_long<-as.data.frame(GGrelLong_c$VCV)
random_effects_long_df <- data.frame(mean=round(colMeans(random_effects_long),3),
                                         lwr5=apply(random_effects_long[1:5],2,quantile,0.05),
                                         upr95=apply(random_effects_long[1:5],2,quantile,0.95))

colnames(random_effects_long_df)<-colnames(fixed_effects_long_df[1:3])
random_effects_long_df[,2:3]<-round(random_effects_long_df[,2:3],3)

##combine fixed and random effects coefficients into one dataframe
rs_long<-rbind(fixed_effects_long_df[1:3],random_effects_long_df)

##combine into one dataframe - Table 1
parameter_names<-c("Intercept",
                   "SexM",
                   "Social rank",
                   "F",
                   "Size of clan",
                   "GGAirstrip",
                   "GGEngitati",
                   "GGForest",
                   "GGLemala",
                   "GGMunge",
                   "GGNgoitokitok",
                   "GGShamba",
                   "GGTriangle",
                   "SexM:Social rank",
                   "VA",
                   "VMother",
                   "VCohort",
                   "VClan",
                   "VR")

tab2<-data.frame(Parameter=parameter_names,
                   LRS=paste0(rs_lrs$post.mean," (",rs_lrs$`l-95% CI`," - ",rs_lrs$`u-95% CI`,")"),
                   Lifespan=paste0(rs_long$post.mean," (",rs_long$`l-95% CI`," - ",rs_long$`u-95% CI`,")"))

#################################
##with heterogeneous variances###
#################################

priorsGG_c <- list(R = list(V = diag(8), nu = 0.002),
                 G = list(G1 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G2 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G3 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G4 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G5 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G6 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G7 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G8 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G9 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G10 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G11 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G12 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*1000)))

GGrelLRS_hetvar_c <- MCMCglmm(lrs_r ~ sex*birth_rank + f + clan_size + A + E + F + L + M + N + T + S,
                         random = ~ IndexX + IndexA + IndexF + IndexE + IndexN + IndexT + IndexS +
                          IndexM + IndexL + mother + birth_year + birth_clan,
                          ginverse = list(IndexX = AinvcX, IndexA = AinvcA, IndexF = AinvcF, 
                                           IndexE = AinvcE, IndexN = AinvcN, IndexT = AinvcT,
                                           IndexS = AinvcS, IndexM = AinvcM, IndexL = AinvcL),
                         rcov=~idh(birth_clan):units,
                         data = lh_ids_crater,
                         family = "gaussian",
                         prior = priorsGG_c,
                         pr=TRUE,
                         nitt = 18000, burn = 3000, thin = 10)

plot(GGrelLRS_hetvar_c)
summary(GGrelLRS_hetvar_c)

############lifespan
GGrelLong_hetvar_c <- MCMCglmm(lifespan_r ~ sex*birth_rank + f + clan_size + A + E + F + L + M + N + T + S,
                               random = ~ IndexX + IndexA + IndexF + IndexE + IndexN + IndexT + IndexS +
                               IndexM + IndexL + mother + birth_year + birth_clan,
                               ginverse = list(IndexX = AinvcX, IndexA = AinvcA, IndexF = AinvcF, 
                                               IndexE = AinvcE, IndexN = AinvcN, IndexT = AinvcT,
                                               IndexS = AinvcS, IndexM = AinvcM, IndexL = AinvcL),
                              rcov=~idh(birth_clan):units,
                              data = lh_ids_crater,
                              family = "gaussian",
                              prior = priorsGG_c,
                              pr=TRUE,
                              nitt = 18000, burn = 3000, thin = 10)

plot(GGrelLong_hetvar_c)
summary(GGrelLong_hetvar_c)

####save model objects - as workspace at this point
save.image(file="modelResults_Robjects/GGmods_ranefsaved.RData")
####reload model objects if coming back later
load(file="modelResults_Robjects/GGmods.RData")


#############################
####make results tables######
#############################

# Extracting fixed effects
fixed_effects_lrs_het <- round(summary(GGrelLRS_hetvar_c)$solutions,3)
fixed_effects_lrs_het_df <- as.data.frame(fixed_effects_lrs_het)

# Extracting random effects
random_effects_lrs_het<-as.data.frame(GGrelLRS_hetvar_c$VCV)
random_effects_lrs_het<-random_effects_lrs_het[c("IndexA","IndexE","IndexF","IndexL","IndexM","IndexN","IndexS","IndexT","IndexX",
                                                   colnames(random_effects_lrs_het[10:20]))]

random_effects_lrs_het_df <- data.frame(mean=round(colMeans(random_effects_lrs_het),3),
                                        lwr5=apply(random_effects_lrs_het[1:20],2,quantile,0.05),
                                        upr95=apply(random_effects_lrs_het[1:20],2,quantile,0.95))

colnames(random_effects_lrs_het_df)<-colnames(fixed_effects_lrs_het_df[1:3])
random_effects_lrs_het_df[,2:3]<-round(random_effects_lrs_het_df[,2:3],3)

###lifespan
# Extracting fixed effects
fixed_effects_long_het <- round(summary(GGrelLong_hetvar_c)$solutions,3)
fixed_effects_long_het_df <- as.data.frame(fixed_effects_long_het)

# Extracting random effects
random_effects_long_het<-as.data.frame(GGrelLong_hetvar_c$VCV)
random_effects_long_het<-random_effects_long_het[c("IndexA","IndexE","IndexF","IndexL","IndexM","IndexN","IndexS","IndexT","IndexX",
                                                   colnames(random_effects_long_het[10:20]))]

random_effects_long_het_df <- data.frame(mean=round(colMeans(random_effects_long_het),3),
                                        lwr5=apply(random_effects_long_het[1:20],2,quantile,0.05),
                                        upr95=apply(random_effects_long_het[1:20],2,quantile,0.95))

colnames(random_effects_long_het_df)<-colnames(fixed_effects_long_het_df[1:3])
random_effects_long_het_df[,2:3]<-round(random_effects_long_het_df[,2:3],3)

#####now compile all the results from models into a single dataframe
parameter_names2<-c("VAAirstrip",
  "VAEngitati",
  "VAForest",
  "VALemala",
  "VAMunge",
  "VANgoitokitok",
  "VAShamba",
  "VATriangle",
  "VANon-crater",
  "VMother",
  "VCohort",
  "VClan",
  "VRAirstrip",
  "VREngitati",
  "VRForest",
  "VRLemala",
  "VRMunge",
  "VRNgoitokitok",
  "VRShamba",
  "VRTriangle")

tab3<-data.frame(Parameter = parameter_names2,
                 LRS = paste0(random_effects_lrs_het_df$post.mean," (",random_effects_lrs_het_df$`l-95% CI`," - ",random_effects_lrs_het_df$`u-95% CI`,")"),
                 lifespan = paste0(random_effects_long_het_df$post.mean," (",random_effects_long_het_df$`l-95% CI`," - ",random_effects_long_het_df$`u-95% CI`,")"))



##Table S1:
fixed_effects_lrs_het_df<-fixed_effects_lrs_het_df[c(rownames(fixed_effects_lrs_het_df[1:11,]),"S","T","sexmale:birth_rank"),]
fixed_effects_long_het_df<-fixed_effects_long_het_df[c(rownames(fixed_effects_long_het_df[1:11,]),"S","T","sexmale:birth_rank"),]

tabS1<-data.frame(Parameter=parameter_names[1:14],
                   LRS=paste0(fixed_effects_lrs_het_df$post.mean," (",fixed_effects_lrs_het_df$`l-95% CI`," - ",fixed_effects_lrs_het_df$`u-95% CI`,")"),
                   Lifespan = paste0(fixed_effects_long_het_df$post.mean," (",fixed_effects_long_het_df$`l-95% CI`," - ",fixed_effects_long_het_df$`u-95% CI`,")"))

write.csv(tab2,file="tables_figures/Table2.csv",row.names = FALSE)
write.csv(tab3,file="tables_figures/Table3.csv",row.names = FALSE)
write.csv(tabS1,file="tables_figures/TableS1.csv",row.names = FALSE)

####################################################
## plot all effects ###############################
####################################################

##plot effects of immigrant ancestry on fitness
##predicting across original dataframe

dummyped<-lh_ids_crater

dummyped_lrspreds <- as.data.frame(predict(GGrelLRS_c, newdata=dummyped,
                                            marginal=GGrelLRS_c$Random$formula,
                                            type="response",
                                            interval="confidence",
                                            level=0.95))


dummyped_longpreds <- as.data.frame(predict(GGrelLong_c, newdata=dummyped,
                                            marginal=GGrelLong_c$Random$formula,
                                            type="response",
                                            interval="confidence",
                                            level=0.95))
dummyped$fitLRS<-dummyped_lrspreds$fit
dummyped$lwrLRS<-dummyped_lrspreds$lwr
dummyped$uprLRS<-dummyped_lrspreds$upr

dummyped$fitL<-dummyped_longpreds$fit
dummyped$lwrL<-dummyped_longpreds$lwr
dummyped$uprL<-dummyped_longpreds$upr

pdf("tables_figures/Figure3A.pdf")
ggplot(dummyped,aes(x=X,y=lrs_r))+
  geom_jitter(alpha=1/4)+
  geom_smooth(method = "lm", se = F,colour="black",alpha=1/2) +
  geom_line(aes(y = stage(start = lrs_r, after_stat = ymin)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  geom_line(aes(y = stage(lrs_r, ymax)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  xlab("\nProportion of immigrant ancestry")+
  ylab("LRS (relative to mean)\n\n")+
  theme_classic()+
  theme(text = element_text(size=14))
dev.off()

pdf("tables_figures/Figure3B.pdf")
ggplot(dummyped,aes(x=X,y=lifespan_r))+
  geom_jitter(alpha=1/4)+
  geom_smooth(method = "lm", se = F,colour="black",alpha=1/2) +
  geom_line(aes(y = stage(start = lifespan_r, after_stat = ymin)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  geom_line(aes(y = stage(lifespan_r, ymax)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  xlab("\nProportion of immigrant ancestry")+
  ylab("Lifespan (in years relative to mean)\n")+
  theme_classic()+
  theme(text = element_text(size=14))
dev.off()

##########################################
###plot posteriors for Va across clans####
##########################################

oCols<-c("darkorange","mediumorchid3","grey60","gray15","chartreuse3","tomato2","darkcyan","gold","royalblue")

#####homogeneous Va
#lrs
var_lrsC<-as.data.frame(GGrelLRS_c$VCV)
var_lrsC_l<-melt(var_lrsC[1])
homVlrs<-data.frame(intercepts=c((quantile(var_lrsC_l$value,0.05)),
                                 (median(var_lrsC_l$value)),
                                  (quantile(var_lrsC_l$value,0.95))),
                    names=c("lower5%","median","upper5%"))

#lifespan
var_longC<-as.data.frame(GGrelLong_c$VCV)
var_longC_l<-melt(var_longC[1])
homVlong<-data.frame(intercepts=c((quantile(var_longC_l$value,0.05)),
                                  (median(var_longC_l$value)),
                                  (quantile(var_longC_l$value,0.95))),
                     names=c("lower5%","median","upper5%"))

#####heterogeneous Va
#lrs
var_lrshetvarC<-as.data.frame(GGrelLRS_hetvar_c$VCV)
var_lrshetvarC_l<-melt(var_lrshetvarC[1:9])
var_lrshetvarC_l$variable<-factor(var_lrshetvarC_l$variable,
          levels = c("IndexX","IndexT","IndexS","IndexN","IndexM","IndexL","IndexF","IndexE","IndexA"))

#lifespan
var_longhetvarC<-as.data.frame(GGrelLong_hetvar_c$VCV)
var_longhetvarC_l<-melt(var_longhetvarC[1:9])
var_longhetvarC_l$variable<-factor(var_longhetvarC_l$variable,
           levels = c("IndexX","IndexT","IndexS","IndexN","IndexM","IndexL","IndexF","IndexE","IndexA"))

########################################
####extract residual variances##########
########################################

#####homogeneous Vr
#lrs
var_lrsC_l_r<-melt(var_lrsC[5])
homVlrs_r<-data.frame(intercepts=c((quantile(var_lrsC_l_r$value,0.05)),
                                   (median(var_lrsC_l_r$value)),
                                   (quantile(var_lrsC_l_r$value,0.95))),
                      names=c("lower5%","median","upper5%"))

#lifespan
var_longC_l_r<-melt(var_longC[5])
homVlong_r<-data.frame(intercepts=c((quantile(var_longC_l_r$value,0.05)),
                                    (median(var_longC_l_r$value)),
                                    (quantile(var_longC_l_r$value,0.95))),
                       names=c("lower5%","median","upper5%"))

#####heterogeneous Vr
#lrs
var_lrshetvarC_l_r<-melt(var_lrshetvarC[13:20])
var_lrshetvarC_l_r$variable<-factor(var_lrshetvarC_l_r$variable,
                                    levels = c("birth_clanT.units","birth_clanS.units","birth_clanN.units","birth_clanM.units","birth_clanL.units","birth_clanF.units","birth_clanE.units","birth_clanA.units"))

#lifespan
var_longhetvarC_l_r<-melt(var_longhetvarC[13:20])
var_longhetvarC_l_r$variable<-factor(var_longhetvarC_l_r$variable,
                                    levels = c("birth_clanT.units","birth_clanS.units","birth_clanN.units","birth_clanM.units","birth_clanL.units","birth_clanF.units","birth_clanE.units","birth_clanA.units"))

###########################
#######plot variances######
###########################

###Va

##lrs
pdf("tables_figures/Figure4A.pdf")
ggplot(var_lrshetvarC_l,aes(y=variable,color = variable, x = value)) +
  geom_vline(data=homVlrs,aes(xintercept=intercepts,linetype=names),alpha=1/2,col="black",key_glyph = "path")+
  scale_linetype_manual(values=c("median"="dashed","lower5%"="dotted","upper5%"="dotted"), 
                        labels=c("Posterior median","95% CI posterior"),
                        breaks=c("median","lower5%"),name = expression(V[A]^~~~LRS~~paste("N'goro"))) +
  stat_pointinterval(quantiles = 50)+
  theme_classic()+
  scale_color_manual(name="",guide="none",values=oCols)+
  ylab("")+
  xlab(expression(V[A]^~~~LRS))+
  scale_y_discrete(labels=c("Non-crater","Triangle","Shamba","Ngoitokitok","Munge","Lemala","Forest","Engitati","Airstrip"))+
  theme(text=element_text(size=14),
        legend.justification = c(1, 1), legend.position = c(0.98, 0.97),
        legend.title.align=0.5)
dev.off()

##lifespan
pdf("tables_figures/Figure4B.pdf")
ggplot(var_longhetvarC_l,aes(y=variable,color = variable, x = value)) +
  geom_vline(data=homVlong,aes(xintercept=intercepts,linetype=names),alpha=1/2,col="black",key_glyph = "path")+
  scale_linetype_manual(values=c("median"="dashed","lower5%"="dotted","upper5%"="dotted"), 
                        labels=c("Posterior median","95% CI posterior"),
                        breaks=c("median","lower5%"),name = expression(V[A]^~~~Lifespan~~paste("N'goro"))) +
  stat_pointinterval()+
  theme_classic()+
  scale_color_manual(name="",guide="none",values=oCols)+
  ylab("")+
  xlab(expression(V[A]^~~~Lifespan))+
  scale_y_discrete(labels=c("Non-crater","Triangle","Shamba","Ngoitokitok","Munge","Lemala","Forest","Engitati","Airstrip"))+
  theme(text=element_text(size=14),
        legend.justification = c(1, 1), legend.position = c(0.98, 0.97),
        legend.title.align=0.5)
dev.off()

###Vr
oCols2<-c("mediumorchid3","grey60","gray15","chartreuse3","tomato2","darkcyan","gold","royalblue")

##lrs
pdf("tables_figures/FigureS5A.pdf")
ggplot(var_lrshetvarC_l_r,aes(y=variable,color = variable, x = value)) +
  geom_vline(data=homVlrs_r,aes(xintercept=intercepts,linetype=names),alpha=1/2,col="black",key_glyph = "path")+
  scale_linetype_manual(values=c("median"="dashed","lower5%"="dotted","upper5%"="dotted"), 
                        labels=c("Posterior median","95% CI posterior"),
                        breaks=c("median","lower5%"),name = expression(V[R]^~~~LRS~~paste("N'goro"))) +
  stat_pointinterval()+
  theme_classic()+
  scale_color_manual(name="",guide="none",values=oCols2)+
  ylab("")+
  xlab(expression(V[R]^~~~LRS))+
  scale_y_discrete(labels=c("Triangle","Shamba","Ngoitokitok","Munge","Lemala","Forest","Engitati","Airstrip"))+  
  theme(text=element_text(size=14),
        legend.justification = c(1, 1), legend.position = c(0.98, 0.97),
        legend.title.align=0.5)
dev.off()

##lifespan
pdf("tables_figures/FigureS5B.pdf")
ggplot(var_longhetvarC_l_r,aes(y=variable,color = variable, x = value)) +
  geom_vline(data=homVlong_r,aes(xintercept=intercepts,linetype=names),alpha=1/2,col="black",key_glyph = "path")+
  scale_linetype_manual(values=c("median"="dashed","lower5%"="dotted","upper5%"="dotted"), 
                        labels=c("Posterior median","95% CI posterior"),
                        breaks=c("median","lower5%"),name = expression(V[R]^~~~Lifespan~~paste("N'goro"))) +
  stat_pointinterval()+
  theme_classic()+
  scale_color_manual(name="",guide="none",values=oCols2)+
  ylab("")+
  xlab(expression(V[R]^~~~Lifespan))+
  scale_y_discrete(labels=c("Triangle","Shamba","Ngoitokitok","Munge","Lemala","Forest","Engitati","Airstrip"))+
  theme(text=element_text(size=14),
        legend.justification = c(1, 1), legend.position = c(0.98, 0.97),
        legend.title.align=0.5)
dev.off()


###################################################################################
## calculate posterior of differences in Va (LRS) between pairs of clans ##########
## and plot the posteriors as hisotgrams for supplement ###########################
###################################################################################

#############
##for lrs####
#############

head(var_lrshetvarC)
diffs<-list()

clanNames<-c("Non-crater","Airstrip","Forest","Engitati","Ngoitokitok",
             "Triangle","Shamba","Munge","Lemala")
  
for (i in 1:9){
  IDtmp<-clanNames[[i]]
  tmp<-list()
  for(j in 1:9){
  tmp[[j]]<-data.frame(var_lrshetvarC[i]-var_lrshetvarC[j])
  names(tmp[[j]])<-paste0(IDtmp,"_",clanNames[[j]])
  }
  diffs[[i]] <- do.call("cbind",tmp[-i])
  
  ##plot them as multipanel histogram
  pdf(paste0("tables_figures/FigS6_",IDtmp,"_LRS.pdf"))
  par(mfrow=c(4,2))
  for(k in 1:length(diffs[[i]])){
    hist(diffs[[i]][,k],breaks=30,main=names(diffs[[i]][k]),xlab="Difference in VA for LRS",col=alpha("cornflowerblue",0.2))
    abline(v=0,col="darkorange",lwd=2)
  }
  dev.off()
  }

##calculate summary statistics for differences between clans in Va
diffdDF<-do.call("cbind",diffs)
pd_lrs_mode<-posterior.mode(as.mcmc(diffdDF[1:72]))
pd_lrs_95<-t(apply(diffdDF,2,quantile,probs=c(0.05,0.95)))

pd_lrs<-cbind(pd_lrs_mode,pd_lrs_95)
write.csv(pd_lrs,file="tables_figures/pd_prs.csv")

#############
##lifespan####
#############

diffs_lifespan<-list()

for (i in 1:9){
  IDtmp<-clanNames[[i]]
  tmp<-list()
  for(j in 1:9){
    tmp[[j]]<-data.frame(var_longhetvarC[i]-var_longhetvarC[j])
    names(tmp[[j]])<-paste0(IDtmp,"_",clanNames[[j]])
  }
  diffs_lifespan[[i]] <- do.call("cbind",tmp[-i])
  
  ##plot them as multipanel histogram
  pdf(paste0("tables_figures/FigS7_",IDtmp,"_Lifespan.pdf"))
  par(mfrow=c(4,2))
  for(k in 1:length(diffs_lifespan[[i]])){
    hist(diffs_lifespan[[i]][,k],breaks=30,main=names(diffs_lifespan[[i]][k]),xlab="Difference in VA for Lifespan",col=alpha("cornflowerblue",0.2))
    abline(v=0,col="darkorange",lwd=2)
  }
  dev.off()
}

##calculate summary statistics for differences between clans in Va
diffs_lifespandDF<-do.call("cbind",diffs_lifespan)
pd_lifespan_mode<-posterior.mode(as.mcmc(diffs_lifespandDF[1:72]))
lifespan_pd_95<-t(apply(diffs_lifespandDF,2,quantile,probs=c(0.05,0.95)))

pd_lifespan<-cbind(pd_lifespan_mode,lifespan_pd_95)
write.csv(pd_lifespan,file="tables_figures/pd_lifespan.csv")

