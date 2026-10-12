
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
gl<-dartR::gl.read.vcf("data/Hyenas_2ndRun_mincov10MaxCov110menderr01maf1miss70_dupsremovalidswap.recode.vcf")
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

pdf("tables_figures/R2/Fig2A.pdf")
plot(hPCA$scores[,1],hPCA$scores[,2],type="n",cex.lab=1.2,cex.axis=1,
            xlab=paste0("First PC (", round(100 * var_frac[1], 1), "%)"),
            ylab=paste0("Second PC (", round(100 * var_frac[2], 1), "%)"))
s.class(hPCA$scores,fac=as.factor(meta_geneticdata2$clan_birth),add.plot=T,col=alpha(oCols,0.5))
dev.off()

pdf("tables_figures/R2/FigS1.pdf")
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
clHet<-gl.report.heterozygosity(glC,nboots = 1000,ncpus = 12,CI.type="perc")
popHet_test<-gl.test.heterozygosity(glC)

##make table
Table1<-clHet[c("pop","Ho","HoSD","uHe","uHeSD","FIS","FISLCI","FISHCI")]
Table1$pop<-c("Airstrip","Engitati","Forest","Lemala","Munge","Ngoitokitok","Shamba","Triangle","Non-crater")
Table1<-data.frame(Clan=Table1$pop,
                   Heterozygosity = paste0(round(Table1$Ho,3)," (",round(Table1$HoSD,3),")"),
                   Expected_heterozygosity = paste0(round(Table1$uHe,3)," (",round(Table1$uHeSD,3),")"),
                   Fis = paste0(round(Table1$FIS,3)," (",round(Table1$FISLCI,3)," - ",round(Table1$FISHCI,3),")"))
write.csv(Table1,"tables_figures/R2/Table1.csv",row.names = FALSE)

#########################################
##calculate pairwise fsts between clans##
#########################################

fst<-StAMPP::stamppFst(glC,nboots = 1000, nclusters = 5,percent=95)
##use p-values from boostrapping here to report in ms

##save confidence intervals
fst_boot<-fst$Bootstraps
# create pair names like "Pop1-Pop2"
pair_names <- paste(fst_boot[,1], fst_boot[,2], sep = "-")

# extract bootstrap columns
fst_values <- as.matrix(fst_boot[, -(1:2)])

# transpose so:
# rows = bootstrap replicates
# cols = population pairs
fst_mat <- t(fst_values)
colnames(fst_mat) <- pair_names

##calculate ∆fst between all population pairs across bootstraps to identify which are sig. different from each other
fst_pairs <- combn(colnames(fst_mat), 2, simplify = FALSE)

delta_mat <- sapply(fst_pairs, function(pair) {
  fst_mat[, pair[1]] - fst_mat[, pair[2]]
})

colnames(delta_mat) <- sapply(fst_pairs, function(pair) {
  paste(pair[1], "vs", pair[2])
})

alpha <- 0.001

ci_lower <- apply(delta_mat, 2, quantile, probs = alpha/2, na.rm = TRUE)
ci_upper <- apply(delta_mat, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)

ci_results <- data.frame(
  comparison = colnames(delta_mat),
  lower = ci_lower,
  upper = ci_upper
)

ci_results$significant <- !(ci_results$lower <= 0 & ci_results$upper >= 0)
ci_results$direction <- ifelse(ci_results$lower > 0, "first > second",
                               ifelse(ci_results$upper < 0, "first < second", "ns"))

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

# make pair IDs consistent with your matrix ("A:B" format)
pair_id <- paste(fst_boot[,1], fst_boot[,2], sep=":")

# extract bootstrap values
boot_vals <- as.matrix(fst_boot[, -(1:2)])

# compute stats
fst_ci_df <- data.frame(
  pair = pair_id,
  mean = rowMeans(boot_vals, na.rm = TRUE),
  lower = apply(boot_vals, 1, quantile, probs = 0.05, na.rm = TRUE),
  upper = apply(boot_vals, 1, quantile, probs = 0.95, na.rm = TRUE)
)

fst_m$pair <- paste(fst_m$Var1, fst_m$Var2, sep=":")

fst_ci_df$pair_sorted <- apply(
  cbind(fst_boot[,1], fst_boot[,2]),
  1,
  function(x) paste(sort(x), collapse=":")
)

fst_m$pair_sorted <- apply(
  cbind(as.character(fst_m$Var1), as.character(fst_m$Var2)),
  1,
  function(x) paste(sort(x), collapse=":")
)

# merge CI info
fst_m <- merge(fst_m, fst_ci_df[, c("pair_sorted","lower","upper")],
               by="pair_sorted", all.x=TRUE)

fst_m$label <- ifelse(
  fst_m$Var1 == fst_m$Var2,
  "",
  sprintf("%.3f\n[%.3f–%.3f]", fst_m$value, fst_m$lower, fst_m$upper)
)

##plot

pdf("tables_figures/R2/Fig2B.pdf",width=11,height=10)
ggplot(data = fst_m, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(name="Fst\n",low = "#2166AC", high = "#D6604D", mid = "#F7F7F7", 
                       midpoint = 0.02, limit = c(0,0.05))+
  geom_text(aes(label = label), size = 4, show.legend = FALSE)+
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

##calculate pairwise fsts between clans post dispersal WITH bootstrapping
fstCd <- StAMPP::stamppFst(glCd, nboots = 1000, nclusters = 5, percent = 95)
fstCdV<-na.omit(as.data.frame(unmatrix(fstCd$Fsts)))
fstCdV$origin<-sapply(strsplit(rownames(fstCdV),":"), `[`, 1)
fstCdV$destination<-sapply(strsplit(rownames(fstCdV),":"), `[`, 2)

##save confidence intervals
fstCd_boot<-fstCd$Bootstraps
# create pair names like "Pop1-Pop2"
pair_names_fstCd <- paste(fstCd_boot[,1], fstCd_boot[,2], sep = "-")

# extract bootstrap columns
fstCd_values <- as.matrix(fstCd_boot[, -(1:2)])

# transpose so:
# rows = bootstrap replicates
# cols = population pairs
fstCd_mat <- t(fstCd_values)
colnames(fstCd_mat) <- pair_names

##calculate ∆fst between all population pairs across bootstraps to identify which are sig. different from each other
fstCd_pairs <- combn(colnames(fstCd_mat), 2, simplify = FALSE)

delta_mat_fstCd <- sapply(fstCd_pairs, function(pair) {
  fstCd_mat[, pair[1]] - fstCd_mat[, pair[2]]
})

colnames(delta_mat_fstCd) <- sapply(fstCd_pairs, function(pair) {
  paste(pair[1], "vs", pair[2])
})

alpha <- 0.05

ci_lower_fstCd <- apply(delta_mat_fstCd, 2, quantile, probs = alpha/2, na.rm = TRUE)
ci_upper_fstCd <- apply(delta_mat_fstCd, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)

ci_results_fstCd <- data.frame(
  comparison = colnames(delta_mat_fstCd),
  lower = ci_lower_fstCd,
  upper = ci_upper_fstCd
)

ci_results_fstCd$significant <- !(ci_results_fstCd$lower <= 0 & ci_results_fstCd$upper >= 0)

ci_results_fstCd$direction <- ifelse(ci_results_fstCd$lower > 0, "first > second",
                               ifelse(ci_results_fstCd$upper < 0, "first < second", "ns"))

##extract bootstrap results
fstCd_boot <- fstCd$Bootstraps

# create pair IDs
pair_id_Cd <- paste(fstCd_boot[,1], fstCd_boot[,2], sep=":")

# extract bootstrap values
boot_vals_Cd <- as.matrix(fstCd_boot[, -(1:2)])

# compute CI stats
fstCd_ci_df <- data.frame(
  pair = pair_id_Cd,
  mean = rowMeans(boot_vals_Cd, na.rm = TRUE),
  lower = apply(boot_vals_Cd, 1, quantile, probs = 0.05, na.rm = TRUE),
  upper = apply(boot_vals_Cd, 1, quantile, probs = 0.95, na.rm = TRUE)
)

##make sorted pair ID for symmetric matching
fstCd_ci_df$pair_sorted <- apply(
  cbind(fstCd_boot[,1], fstCd_boot[,2]),
  1,
  function(x) paste(sort(x), collapse=":")
)

fstCd_ci_df$significant <- !(fstCd_ci_df$lower <= 0 & fstCd_ci_df$upper >= 0)

fstCd_ci_df$direction <- ifelse(fstCd_ci_df$lower > 0, "first > second",
                               ifelse(fstCd_ci_df$upper < 0, "first < second", "ns"))

##format FST matrix (same as before)
fstCd$Fsts[upper.tri(fstCd$Fsts)] <- t(fstCd$Fsts)[upper.tri(fstCd$Fsts)]
fstCd_a <- fstCd$Fsts[clan_order, clan_order]

##keep only upper triangle
fstCd_a[lower.tri(fstCd_a)] <- NA
diag(fstCd_a) <- 0

##melt for plotting
fstCd_m <- reshape2::melt(fstCd_a)
fstCd_m <- na.omit(fstCd_m)
fstCd_m$value <- as.numeric(fstCd_m$value)

##factor ordering 
fstCd_m$Var2 <- factor(fstCd_m$Var2, levels=c("X","T","S","N","M","L","F","E","A"))
fstCd_m$Var1 <- factor(fstCd_m$Var1, levels=c("A","E","F","L","M","N","S","T","X"))

##create sorted pair IDs for merge
fstCd_m$pair_sorted <- apply(
  cbind(as.character(fstCd_m$Var1), as.character(fstCd_m$Var2)),
  1,
  function(x) paste(sort(x), collapse=":")
)

##merge CI info (robust to A-B vs B-A)
fstCd_m <- merge(
  fstCd_m,
  fstCd_ci_df[, c("pair_sorted","lower","upper")],
  by = "pair_sorted",
  all.x = TRUE
)

##create label with CI
fstCd_m$label <- ifelse(
  fstCd_m$Var1 == fstCd_m$Var2,
  "",
  sprintf("%.3f\n[%.3f–%.3f]", fstCd_m$value, fstCd_m$lower, fstCd_m$upper)
)

##plot
pdf("tables_figures/R2/FigS2.pdf",width=9,height=9)
ggplot(data = fstCd_m, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  geom_text(aes(label = label), size = 3, show.legend = FALSE) +
  scale_fill_gradient2(
    name="Fst\n",
    low = "#2166AC", high = "#D6604D", mid = "#F7F7F7",
    midpoint = 0.02, limit = c(0,0.05)
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(text=element_text(size=16), legend.title.align=0.5) +
  theme(
    legend.position = c(.98, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) +
  xlab("\nClan") +
  ylab("Clan\n") +
  coord_fixed()
dev.off()

#################################################
####test relationship between spatial distance, dispersal and fst####
#################################################

####plot pairwise dispersal as heatmap
##need to reorganise to keep zeros and have the rate of philopatry on the diagonal
allDisp2$destination<-factor(allDisp2$destination,levels=c("X","T","S","N","M","L","F","E","A"))
allDisp2$origin<-factor(allDisp2$origin,levels=c("A","E","F","L","M","N","S","T","X"))

allDisp3<-as_tibble(allDisp2)
allDisp3$origin<-as.factor(allDisp3$origin)
allDisp3$destination<-as.factor(allDisp3$destination)

allDisp3 = allDisp3 %>%
  group_by(origin) %>%
  tidyr::complete(destination = unique(allDisp3$destination))

#pdf("tables_figures/R2/Figure3A.pdf",width=8,height=8)
fig3A<-ggplot(data = allDisp3, aes(x=origin, y=destination, fill=probDisp)) + 
  geom_tile()+
  scale_fill_gradient2(low = "#2166AC", high = "#D6604D", 
                       midpoint = 0.12, limit = c(0,0.6))+
  geom_text(data = allDisp2,aes(label = round(as.numeric(probDisp),3),size=5),show.legend = FALSE)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=16),
        legend.title=element_blank())+
  xlab("\nOrigin")+
  ylab("Destination\n")+
coord_fixed()
#dev.off()

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

##run mrqap across bootstraps of fst
# function to convert one bootstrap row → matrix
make_fst_matrix <- function(fst_row, pair_names, pops) {
  mat <- matrix(NA, length(pops), length(pops),
                dimnames = list(pops, pops))
  for(i in seq_along(pair_names)) {
    pair <- strsplit(pair_names[i], "-")[[1]]
    p1 <- trimws(pair[1])
    p2 <- trimws(pair[2])
    if(!(p1 %in% pops && p2 %in% pops)) next
    mat[p1, p2] <- fst_row[i]
    mat[p2, p1] <- fst_row[i]
  }
  
  diag(mat) <- 0
  mat
}

pops <- sort(setdiff(unique(disp_fstM$origin), "X"))
nboots <- nrow(fst_mat)

# store results
mrqap_results <- data.frame(
  beta_disp = numeric(nboots),
  beta_dist = numeric(nboots),
  beta_distOnly = numeric(nboots)
)

# rebuild pair names safely
pair_split <- strsplit(colnames(fst_mat), "-")
keep <- sapply(pair_split, function(x) !("X" %in% x))

fst_mat <- fst_mat[, keep, drop = FALSE]

for(i in 1:nboots){
  
  fst_boot_mat <- make_fst_matrix(
    fst_mat[i, ],
    colnames(fst_mat),
    pops
  )
  
  # ensure same ordering as other matrices
  fst_boot_mat <- fst_boot_mat[rownames(dispM), colnames(dispM)]

  fit <- asnipe::mrqap.dsp(
    fst_boot_mat ~ dispM + udDist,
    test.statistic = "beta"
  )
  
  fit_ud <- asnipe::mrqap.dsp(
    fst_boot_mat ~ udDist,
    test.statistic = "beta"
  )

  
  mrqap_results$beta_disp[i] <- fit$coefficients[2]
  mrqap_results$beta_dist[i] <- fit$coefficients[3]
  mrqap_results$beta_distOnly = fit_ud$coefficients[2]
  
}

alpha <- 0.05

summary_mrqap <- data.frame(
  predictor = c("dispM", "udDist", "udDist_only"),
  mean = c(mean(mrqap_results$beta_disp),
           mean(mrqap_results$beta_dist),
           mean(mrqap_results$beta_distOnly)),
  lower = c(quantile(mrqap_results$beta_disp, alpha/2),
            quantile(mrqap_results$beta_dist, alpha/2),
            quantile(mrqap_results$beta_distOnly, alpha/2)),
  upper = c(quantile(mrqap_results$beta_disp, 1 - alpha/2),
            quantile(mrqap_results$beta_dist, 1 - alpha/2),
            quantile(mrqap_results$beta_distOnly, 1 - alpha/2))
)

#########plot results
df_fst_plots<-data.frame(disp=as.vector(dispM),fst=as.vector(fst_M),dist=as.vector(udDist))

pdf("tables_figures/R2/FigureS4.pdf")
ggplot(df_fst_plots,aes(x=dist,y=disp))+
  geom_point()+
  theme_classic()+
  xlab("\nGeographic distance between clans")+
  ylab("Probability of dispersal between clans\n")+
  theme(text = element_text(size=14),
        plot.margin = unit(c(1,2,1,2), "lines"))
dev.off()

fig3C<-ggplot(df_fst_plots,aes(x=disp,y=fst))+
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
  theme(text = element_text(size=16),
        plot.margin = unit(c(1,2,1,2), "lines"))+
  ylim(0.01,0.045)

pdf("tables_figures/R2/Fig2D.pdf",width=9,height=9)
ggplot(df_fst_plots,aes(x=log(dist),y=(fst/(1-fst))))+
  geom_point()+
  geom_smooth(method = "lm", se = F,colour="black",alpha=1/2) +
  geom_line(aes(y = stage(start = fst/(1-fst), after_stat = ymin)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  geom_line(aes(y = stage(fst/(1-fst), ymax)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  #geom_line(stat="smooth",method="glm",se=FALSE)+
  theme_classic()+
  xlab("\nln(Geographic distance between clans (m))")+
  ylab("Fst/(1 - Fst)\n")+
  theme(text = element_text(size=16),
        plot.margin = unit(c(1,2,1,2), "lines"))+
  ylim(0.01,0.045)
dev.off()

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
  popRd<-meta_geneticdata2[c("ID","sex","clan_birth")]
  popRd$sex[is.na(popRd$sex)]<-"unknown"
  clans <- unique(popRd$clan_birth)
  popRd$clan_birth[popRd$sex=="male"]<-sample(clans, sum(popRd$sex=="male"), replace=TRUE)
  glTmp2$pop<-as.factor(popRd$clan_birth)
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

## get range of bootstrapped FSt's to calculate overlap with above null
fst_randDOdf$pair<-paste0(fst_randDOdf$origin,":",fst_randDOdf$destination)

nullD<-as.data.frame(t(apply(fst_randDOdf[,3:1002],1,
                      quantile,probs=c(0.05,0.95))))
nullD$pair<-rownames(nullD)
nullD_qs<-merge(fst_ci_df,nullD,by="pair")

nullD_qs$obs_rand<-ifelse(nullD_qs$lower>nullD_qs$`95%`,"TRUE","FALSE")

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

fig3B<-ggplot(nullD_qs,aes(x=nullD_qs$'95%',y=mean))+
     geom_point()+
     xlab("\nExpected Fst under random dispersal (upper 95% quantile)")+
     ylab("Fst between clans\n")+
     theme_classic()+
     theme(text = element_text(size=16),
           plot.margin = unit(c(1,2,1,2), "lines"))+
     xlim(0.005,0.045)+
     ylim(0.01,0.045)


###############################
######IBD analysis#############
###############################

##genetic distance matrix [nei distance] -- takes a while to run
gen<-StAMPP::stamppNeisD(gl2,FALSE)

# Ensure symmetric distance matrices by matching colnames to rownames
colnames(gen)<-rownames(gen)
idDist<-as.matrix(idDist)
colnames(idDist)<-rownames(idDist)
gen2<-gen[rownames(gen)%in%rownames(idDist),colnames(gen)%in%colnames(idDist)]
idDist2<-idDist[rownames(gen2),colnames(gen2)]

##check order
stopifnot(identical(rownames(gen2), rownames(idDist2)))

##run isolation-by-distance analyses
diag(gen2)<-NA
diag(idDist2)<-NA
vegan::mantel(gen2,idDist2)
vegan::mantel((fst_M/(1-fst_M)),log(udDist))

##plot isolation-by-distance
pedH<-read.csv(file="data/pedigree_clans.csv",row.names = 1)
Amat <- as.matrix(nadiv::makeA(pedH[, 1:3]))
AmtV<-data.frame(unmatrix(Amat))
AmtV$pair<-rownames(AmtV)

ibd_df<-data.frame(dist=as.vector(idDist2),
                   gen=as.vector(gen2),
                   pair=rownames(data.frame(unmatrix(gen2))))

ibd_dfA<-left_join(AmtV[,c(2,1)],ibd_df)

pdf("tables_figures/R2/fig2C.pdf",width=11)
ggplot(ibd_df,aes(x=dist,y=gen))+
  geom_point(col="darkgrey",alpha=1/7)+
  geom_smooth(method = "lm", se = F,colour="black",alpha=1/2) +
  geom_line(aes(y = stage(start = gen, after_stat = ymin)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  geom_line(aes(y = stage(gen, ymax)),
            stat = "smooth", method = "lm", lty = "dashed",alpha=1/2) +
  theme_classic()+
  xlab("\nGeographic distance between individuals (m)")+
  ylab("Genetic distance between individuals (Nei's D)\n")+
  theme(text = element_text(size=16),
        plot.margin = unit(c(1,2,1,2), "lines"))
dev.off()

pdf("tables_figures/R2/Figure3.pdf",height=23,width=10)
plot_grid(fig3A,fig3B,fig3C,ncol=1,labels = c("A","B","C"),
          align = "v")
dev.off()

save.image("modelResults_Robjects/popgen_R1.RData")

