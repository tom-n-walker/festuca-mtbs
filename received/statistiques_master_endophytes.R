#Script travail de master

#Endophytes#############################################################################
########################################################################################

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tmap)
library(raster)
library(lme4)
library(lmerTest)
library(ape)
library(plot.matrix)
library(ade4)
library(heatmaply)
library(reshape2)
library(vegan)
library(rgl)
library(phytools)

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Endophytes 2/Extract") 

champi <- read.table(file="Festuca_endophytes_Illumina_ITS_1285332544.all_Fungi_ASVs.tsv", sep='\t', header=TRUE)

### Phylogeny
setwd("D:/Affaires Nadline/Master thesis/Statistiques/Metabolomic")
tree <- read.tree("ITS_FestucaSwissIQTREE.contree")

### climatic niche 

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Metabolomic")
Taxon_niche = read.csv("Taxon_Predictors_Ranges.csv", sep = ",", h = T)

### On enlève les colonnes qu'on ne veut pas

champi2 <- champi[,-2]
champi2 <- champi2[,-99]
champi2 <- champi2[,-98]
champi2 <- champi2[,-97]
champi2 <- champi2[,-96]
champi2 <- champi2[,-95]
champi2 <- champi2[,-94]
champi2 <- champi2[,-93]
champi2 <- champi2[,-1]
champi22 <- champi2

#Hellinger
library(labdsv)
str(champi22)
champi22 <- apply(champi22, 2, as.numeric)
fix(champi22)
champi3 <- hellinger(champi22)
fix(champi3)

### Transposer la matrice

champi33 <- t(champi3)
class(champi33)
dim(champi33)
fix(champi33)

#Importer tableau avec les vrais ID ------
setwd("D:/Affaires Nadline/Master thesis/Statistiques/Endophytes 2") 

FestucaDNA = read.csv("Festuca_DNA.csv", sep = ";", h = T)

row.names(FestucaDNA) <- FestucaDNA$ID_DNA #On rattache la colonne "ID_Metabo" dans le tableau à la fonction
#row.names, afin de pouvoir ensuite merger par "row.names" puisque ce sera en commun

#Merger les 2 matrices --------------------

champi4 = merge(FestucaDNA,champi33, by = "row.names")
fix(champi4)

champi4<-champi4[!(champi4$PlantID=="35" | champi4$PlantID=="36"),]
champi4<-champi4[!(champi4$PlantID=="51" | champi4$PlantID=="90"),]
champi4<-champi4[!(champi4$PlantID=="135" | champi4$PlantID=="137"),]

### Faire la somme des OTUs 
champi_all <- champi4 #Rubra et JB
champi_all <- champi4[is.na(champi4$Site),]#Que JB

#Que JB
champi_all <- champi4[is.na(champi4$Site),]
#On regarde à partir de quelle colonne il y a les endophytes:
colnames(champi_all)
champi_all2 <- champi_all[,26:ncol(champi_all)] #Nom des endophytes commence à partir de colonne 29 incluse jusqu'à la fin
div_endophytes <- rowSums(champi_all2) #On fait la somme des endophytes

#Que Rubra
champi_rubra <- champi4[!is.na(champi4$Site),]
champi_rubra2 <- champi_rubra[,26:ncol(champi_all)] #Nom des molécules commence à partir de colonne 29 incluse jusqu'à la fin
#champi_rubra2 <- data.matrix(champi_rubra2, rownames.force=NA) 
div_rubra <- rowSums(champi_rubra2)

reg = lmer(div_rubra~Altitude+(1|Transect), data=champi_rubra)
anova(reg) #Significatif !

###############JB BRAY Distance #############################################################
champi_all <- champi4 # include rubra elevation gradient
champi_all2 <- champi_all[,26:ncol(champi_all)]
#champi_all2 <- data.matrix(champi_all2, rownames.force=NA)
div_endophytes <- rowSums(champi_all2)

boxplot(div_endophytes~champi_all$Altitude,las = 2)

### Rajouter la colonne div_sp (si jamais je veux sauvegarder tel quel)
champi_all2_env <- champi_all[,2:25] #on garde les colonnes 3 à 28 pour créer la matrice environnementale
champi_all2_env$div_endo <- div_endophytes #Rajouter la colonne div_endo qui contient les éléments de div_endophytes (voir au dessus)
#write.table(champi_all2_env, "matriceendo_finale.csv", sep = ";")
champi_all2_env$taxon <- gsub(" ","_",champi_all2_env$Plant.Species) #gsub(pattern = "o" , replacement = "O") donc ici
#pour remplacer l'espace entre les noms par un tiret
colnames(champi_all2_env)
#Du coup on a rajouté à la matrice environmentale deux colonnes, div_sp et taxon.

### Bray distance + add pcao metric 
#Bray curtis = structures chimiques proches, pas richesse !!
#(1) JB
champi_all2 <- champi_all[,26:ncol(champi_all)] #on garde seulement les colonnes qui
#contiennent les OTUs. 
#champi_all2 <- data.matrix(champi_all2, rownames.force=NA)
spe.bray <- vegdist(champi_all2) #Calculer la distance de Bray entre les molécules

#Multidimensional scaling (MDS) pour projeter la matrice de distance
#en plusieurs axes ------------------------------------------

spe.b.pcoa <- cmdscale(spe.bray, k=3) #k=3 car on veut une PCOA a 3 axes (dimensions)
spe.b.pcoa <- data.frame(spe.b.pcoa) #Mise en tableau des données brutes (donc en 3 axes distincts qu'on peut voir dans
#"l'environnement" (=panneau à droite de R))


champi_full_sp = merge(champi_all2_env,spe.b.pcoa, by = "row.names")
#On les merge pour pouvoir ensuite les colorer par variables environmentales.
plot(champi_full_sp$X1,champi_full_sp$X2,col=as.factor(champi_full_sp$Plant.Species),pch=16,cex=2)

PC1Species <- data.frame(champi_full_sp$X1, champi_full_sp$div_endo, champi_full_sp$PlantID)
#write.table(PC1Species, "PC1Species.csv", sep = ";")

#(2) Rubra
champi_rubra2_env <- champi_rubra[,2:25]
colnames(champi_rubra)
champi_rubra2 <- champi_rubra[,26:ncol(champi_rubra)] #Nom des endophytes commence à partir de colonne 29 incluse jusqu'à la fin
#champi_all2[champi_all2 > 0] <- 1
champi_rubra2 <- data.matrix(champi_rubra2, rownames.force=NA) #Convertir les valeurs du tableau en valeurs numériques
div_rubra <- rowSums(champi_rubra2)

#champi_rubra2_env$div_endorubra <- div_rubra

spe.bray <- vegdist(champi_rubra2)
spe.b.pcoa <- cmdscale(spe.bray, k=3)
spe.b.pcoa <- data.frame(spe.b.pcoa)
champi_full_sp2 = merge(champi_rubra2_env,spe.b.pcoa, by = "row.names")
plot(champi_full_sp2$X1,champi_full_sp2$X3,col=as.factor(champi_full_sp2$Altitude),pch=16,cex=2)
text(champi_full_sp2$X1,champi_full_sp2$X3, labels = champi_full_sp2$Altitude,cex=0.9)

############################################################################33

###Phylogeny

Plant_Species <- unique(champi_full_sp$Plant.Species)
Plant_Species <- gsub(" ", "_",Plant_Species)

##### species niche #Sélectionner les Festuca sp. dans le tableau complet de "Taxon_niche" et creation d'un data frame
matt_niche = NULL
for ( i in c (1:length(Plant_Species))) {
  
  matt_inter <- Taxon_niche[grep(Plant_Species[i],Taxon_niche$taxon),]
  matt_inter <-  matt_inter[1,]
  matt_inter$taxon <- Plant_Species[i]
  matt_niche <-  rbind(matt_niche,matt_inter)
}

champi_full_sp1 = merge(champi_full_sp,matt_niche, by = "taxon") 
#Donc désormais on possède un tableau nommé Metabo_full_sp1 qui contient la colonne DEM_Mean en plus

########################################################################################

champi_stat <- champi_full_sp1[!is.na(champi_full_sp1$Altitude),] ## delet ! to exclude rubra field data  # ! to keep rubra
#Sélectionner tous ceux qui n'ont pas de "NA" dans la colonne "Altitude" car ! signifie "différent de"
#On reprend la matrice Metabo_full_sp1 car elle contient déjà la matrice de distance.

#Calculer la richesse en endophytes par espèce pour pouvoir placer dans la phylogénie
sp_div_mean <- tapply(champi_stat$div_endo ,champi_stat$Plant.Species, mean) ## calculer la moyenne de phytodiversité par espèces
names(sp_div_mean) <- gsub(" ", "_",names(sp_div_mean))
names(sp_div_mean) <- gsub("Festuca_valesciaca", "Festuca_valesiaca",names(sp_div_mean))


#Mettre de la couleur et du texte:

myColorRamp <- function(colors, values) { 
  v <- (values - min(values))/diff(range(values)) 
  x <- colorRamp(colors)(v) 
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255) 
} 


cols <- myColorRamp(c("darkred","orange","darkgreen"),champi_stat$Altitude) ### colour en fonction elevation

plot(champi_stat$Altitude,champi_stat$div_endo,col=cols,pch=16,cex=3)

plot(champi_stat$Altitude,champi_stat$div_endo,col=cols,pch=16,cex=3,
     main = "Impact of elevation on endophytic richness",
     xlab="Altitude (m)", ylab="Endophytic richness",
     cex.axis=2, cex.lab=2, cex.main=2)

text(champi_stat$Altitude,champi_stat$div_endo,labels = champi_stat$Site)

####################################################################
#####Div_sp en fonction des altitudes pour chaque transect##########

champi_rubra_env_plot <- champi_stat[grep("Zermatt", champi_stat$Transect),] #  "Chasseron" "Chasseral" "Zermatt"   "Salgesch"  "Lavey"

plot(champi_rubra_env_plot$Altitude,champi_rubra_env_plot$div_endo, pch=16,cex=3)
text(champi_rubra_env_plot$Altitude,champi_rubra_env_plot$div_endo,labels = champi_rubra_env_plot$Site) 

#Ploter les espèces en fonction de leur optimum de niche (points noirs) vs ploter les F.rubra du gradient (points rouges)
champi_stat_full_sp <- champi_full_sp1[is.na(champi_full_sp1$Site),]
champi_stat_rubra_sp <- champi_full_sp1[!is.na(champi_full_sp1$Altitude),]
#Metabo_stat_rubra_sp <- Metabo_stat_rubra_sp[grep("Zermatt", Metabo_stat_rubra_sp$Transect),] #  "Chasseron" "Chasseral" "Zermatt"   "Salgesch"  "Lavey"

ggplot(champi_stat, aes(x= Altitude, y = div_endo))+
  geom_point()+
  stat_smooth( method="lm", formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", aes(group=1), size=1.5) 

plot(champi_stat$Altitude,champi_stat$div_endo,pch=16,col=1)
champi_stat_rubra_sp <- champi_full_sp[!is.na(champi_full_sp$Site),]
points(champi_stat_rubra_sp$Altitude,champi_stat_rubra_sp$div_endo,pch=16,col="red")

plot(champi_stat$Altitude,champi_stat$div_endo,pch=16,cex=3,
     main = "Impact of elevation on endophytic richness",
     xlab="Altitude (m)", ylab="Endophytic richness",
     cex.axis=2, cex.lab=2, cex.main=2)
points(champi_stat_rubra_sp$Altitude,champi_stat_rubra_sp$div_endo,pch=16,cex=3,col="red")

champi_stat_rubra_sp <- champi_full_sp[!is.na(champi_full_sp$Site),]
points(champi_stat_rubra_sp$Altitude,champi_stat_rubra_sp$div_endo,pch=16,col="red")

#########################################################################################################
############## phylogeny ################################################################################

#creation arbre phylo
#Parmi l'arbre phylogénétique en entier, on sélectionne seulement les noms des espèces qui apparaissent
#dans notre tableau sous "Plant_Species"
tip.labels=tree$tip.label
sub.tips<-setdiff(tip.labels, Plant_Species) #vu qu'on a remplacé précédement l'espace par _
new_tree<-drop.tip(tree, sub.tips)
tr <- root(new_tree, 1)
plot(tr)


### phylo signal chimie richness (mis ensemble en fonction de la diversité composés) // ARBRE PHYLOGENETIQUE
sp_div_mean <- tapply(champi_stat$div_sp ,champi_stat$Plant.Species, mean) ## calculer la moyenne de phytodiversité par espèces
names(sp_div_mean) <- gsub(" ", "_",names(sp_div_mean))
names(sp_div_mean) <- gsub("Festuca_valesciaca", "Festuca_valesiaca",names(sp_div_mean))

obj<-contMap(tr,sp_div_mean,plot=FALSE)
obj<-setMap(obj,invert=TRUE)

plot(obj,type="fan",outline=FALSE)
phylosig(tr,sp_div_mean,method="lambda",test=TRUE) #YEY!
#+ c'est rouge, + il y a de composés.
#Ici signal phylogénétique = structuré par rapport à phylogénie.


#--truc de distance cophenetique et de bray entre arbre phylogenetique et chimique-----------------------
#Bray Curtis regarde aux structures chimiques proches, pas la richesse. Donc savoir s'il y a des profils
#chimiques proches. Normalement, c'est l'axe X1 qui doit expliquer mieux la variance.

dist_phylo <- cophenetic(tr) 
spe.bray <- vegdist(dist_phylo)
spe.b.pcoa <- cmdscale(spe.bray, k=3)
spe.b.pcoa <- data.frame(spe.b.pcoa)
sp_div_matt <- data.frame(sp_div_mean)
matt_phylo_chimie  = merge(spe.b.pcoa,sp_div_matt, by = "row.names")

cols <- myColorRamp(c("red","orange","darkgreen"),matt_phylo_chimie$sp_div_mean) ### colour en fonction elevation

ggplot(data = matt_phylo_chimie, aes(x = X1, y = sp_div_mean))+
  geom_point(col = cols, size = 3)

cor.test(matt_phylo_chimie$X1,matt_phylo_chimie$sp_div_mean)


### phylo signal chimie structure (arbre en fonction de la structure des composés)-----------------------
champi_stat_full_sp <- champi_full_sp1[is.na(champi_full_sp1$Site),]

sp_struc <- tapply(champi_stat_full_sp$X1, champi_stat_full_sp$Plant.Species, mean) ## calculer la moyenne de structure chimique par espèces
names(sp_struc) <- gsub(" ", "_",names(sp_struc))
names(sp_struc) <- gsub("Festuca_valesciaca", "Festuca_valesiaca",names(sp_struc))

obj<-contMap(tr,sp_struc,plot=FALSE)
obj<-setMap(obj,invert=TRUE)

plot(obj,type="fan",outline=FALSE)
phylosig(tr,sp_struc,method="lambda",test=TRUE) #Oh non...

##--tangelgram between trees of distance chimique et distance phylogenetique-----------------------------
library(dendextend)


phylo.dist<-cophenetic(tr)
phylo.dist <- as.dist(phylo.dist)

dist_chimie <- champi_all[,26:ncol(champi_all)]
dist_chimie <- aggregate(dist_chimie,by=list(champi_all$Plant.Species),FUN=mean)
#fix(Metabo_species)
#colnames(Metabo_species)

row.names(dist_chimie) <- dist_chimie$Group.1
row.names(dist_chimie) <- gsub(" ", "_",row.names(dist_chimie))
row.names(dist_chimie) <- gsub("Festuca_valesciaca", "Festuca_valesiaca",row.names(dist_chimie))

dist_chimie <- dist_chimie[,2:ncol(dist_chimie)]
dist_chimie <- vegdist(dist_chimie, na.rm = T)

mantel.rtest(dist_chimie,phylo.dist ,nrepet = 999)


##Mise en dendrogramme des arbres
hc1 <- hclust(dist_chimie, method = "complete")
hc2 <- hclust(phylo.dist, method = "complete")

dend1 <- as.dendrogram(hc1)
dend2 <- as.dendrogram(hc2)

dend12 <- dendlist(dend2, dend1)
cors <- cor.dendlist(dend12)

#creation of the association matrix:
association <- cbind(dend1$tip.label, dend2$tip.label)

cophyloplot(dend1, dend2, assoc = association,
            length.line = 4, space = 28, gap = 3)

tanglegram(dend12,common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
           margin_inner=9,
           lwd=1,remove_nodePar = TRUE,rank_branches = TRUE,type = "triangle")


# dend12 %>% untangle %>% tanglegram
dend12 %>%  untangle(method = "ladderize") %>% tanglegram(common_subtrees_color_lines = TRUE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
                                                          margin_inner=9,
                                                          lwd=0.1,remove_nodePar = TRUE,rank_branches = TRUE,axes = FALSE,lab.cex = 1)


### phylo signal optimum de niche------------------------------------------------------------------------

opti_niche <- matt_niche$DEM_Mean
names(opti_niche) <- matt_niche$taxon

obj<-contMap(tr,opti_niche,plot=FALSE)
obj<-setMap(obj,invert=TRUE)


plot(obj,type="fan",outline=FALSE)
phylosig(tr,opti_niche,method="lambda",test=TRUE)

#(2)
#############################################################################################
#Statistiques OTUs
############################################################################################

rm(list=ls())

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Endophytes 2/Extract") 

champi <- read.table(file = "Festuca_endophytes_Illumina_ITS_1285332544.all_Fungi_ASVS.tsv", sep="\t", header = TRUE)
head(champi)
dim(champi)
str(champi)
plot(champi$DNA_control)
plot(champi$Festuca_acuminata_114)

fix(champi)

plot(champi$total)
mean(champi$total)
min(champi$total)
plot(champi$total, ylim=c(0,5000))

library(vegan)
library(betapart)
library(indicspecies)
library(ade4)
library(factoextra)
library(reshape)

sum(champi$DNA_control)
range(champi$total)
dim(champi)
fix(champi)
matrix <- champi[,3:93]
fix(matrix)
#Nombre d'OTUs différents en tout c'est 730 puisqu'ils s'affichent
#comme le nom des lignes et qu'il y a 730 lignes

#Intra
which( colnames(champi)=="Festuca_rubra_001" )
matrixrubra <- champi[,20:90]
#On enlève Festuca_rubra_105 car issue du JB
fix(matrixrubra)
library(dplyr)
matrixrubra2 <- matrixrubra[apply(matrixrubra[,-1], 1, function(x) !all(x==0)),]
#remove rows which have zero in all columns
fix(matrixrubra2)

#Inter
matrixsp <- champi[,c(3:19,91:93)]
fix(matrixsp)#Festuca_rubra_105 besoin car issue du JB
matrixsp1 <- matrixsp[apply(matrixsp[,-1], 1, function(x) !all(x==0)),]
fix(matrixsp1)

############################################################################################
#Catégories endophytes
############################################################################################

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tmap)
library(raster)
library(lme4)
library(lmerTest)
library(bestNormalize)

#Pour Rubra##############################################################################

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Endophytes 2")

Rubra = read.csv("matriceendo_finalerubra.csv", sep = ";", h = T)

#Part I : SLA, WC, leaf C/N, div_sp
#(1) Normalisation des données
hist(Rubra$div_endo) 
bn <- bestNormalize(Rubra$div_endo)
bn
plot(bn, leg_loc = "bottomright")
Rubra$div_endo <- predict(bn)

hist(Rubra$SLA) 
bn <- bestNormalize(Rubra$SLA)
bn
plot(bn, leg_loc = "bottomright")
Rubra$SLA <- predict(bn)

reg <- lmer(SLA~div_endo+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lm(SLA~div_endo, data=Rubra)
summary(reg)
anova(reg)
##
hist(Rubra$Water.Content) 
bn <- bestNormalize(Rubra$Water.Content)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Water.Content <- predict(bn)

reg <- lmer(Water.Content~div_endo+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lm(Water.Content~div_endo, data=Rubra)
summary(reg)
anova(reg)

Rubra$cat3 <-cut(Rubra$Water.Content, breaks=4, right=FALSE, labels=c(1:4))
ggplot(Rubra, aes(x=as.numeric(cat3), y = div_endo))+
  geom_boxplot(aes(fill = cat3))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^1))+
  labs(title = "Impact of leaf water content on endophytic richness",
       x="Leaf water content (%)",
       y="Endophytic richness")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Leaf water content (%)", 
                   limits=c("Low", "Mid", "High", "Very high"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

ggplot(Rubra, aes(x = Water.Content, y = div_endo)) +
  geom_point(col = "red") +
  geom_smooth(data=Rubra,
              aes(x=Water.Content, y=div_endo), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5) +
  geom_point(data=ClimatSpecies,
             aes(x=Leaf_C.N, y=div_sp), color="turquoise", size=1.5) +
  geom_smooth(data=ClimatSpecies,
              aes(x=Leaf_C.N, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5)


##
hist(Rubra$Leaf_C.N) 
bn <- bestNormalize(Rubra$Leaf_C.N)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Leaf_C.N <- predict(bn)

reg <- lmer(Leaf_C.N~div_endo+(1|Transect), data=Rubra)
summary(reg)
anova(reg)
##
hist(Rubra$div_sp) 
bn <- bestNormalize(Rubra$div_sp)
bn
plot(bn, leg_loc = "bottomright")

Rubra$div_sp <- predict(bn)

reg <- lmer(div_sp~div_endo+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lm(div_sp~div_endo, data=Rubra)
summary(reg)
anova(reg)



########################################################################################
###############Pour les espèces#########################################################

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tmap)
library(raster)
library(lme4)
library(lmerTest)
library(bestNormalize)

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Endophytes 2")

#Extraire la somme totale des endophytes comme il a été fait pour div_sp et ajouter cette colonne dans
#un tableau excel de Rubra

Rubra = read.csv("matriceendo_finalespecies.csv", sep = ";", h = T)

#Part I : SLA, WC, leaf C/N, div_sp
#(1) Normalisation des données

qqnorm(Rubra$div_endo)
qqline(Rubra$div_endo)
hist(Rubra$div_endo) 
bn <- bestNormalize(Rubra$div_endo)
bn
plot(bn, leg_loc = "bottomright")
Rubra$div_endo <- predict(bn)

qqnorm(Rubra$SLA)
qqline(Rubra$SLA)
hist(Rubra$SLA) 
bn <- bestNormalize(Rubra$SLA)
bn
plot(bn, leg_loc = "bottomright")
Rubra$SLA <- predict(bn)

reg <- lmer(SLA~div_endo+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lm(SLA~div_endo, data=Rubra)
summary(reg)
anova(reg)
##
qqnorm(Rubra$Water.Content)
qqline(Rubra$Water.Content)
hist(Rubra$Water.Content) 
bn <- bestNormalize(Rubra$Water.Content)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Water.Content <- predict(bn)

reg <- lmer(Water.Content~div_endo+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lm(Water.Content~div_endo, data=Rubra)
summary(reg)
anova(reg)
##
qqnorm(Rubra$Leaf_C.N)
qqline(Rubra$Leaf_C.N)
hist(Rubra$Leaf_C.N) 
bn <- bestNormalize(Rubra$Leaf_C.N)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Leaf_C.N <- predict(bn)

reg <- lm(Leaf_C.N~div_endo, data=Rubra)
summary(reg)
anova(reg)
##
qqnorm(Rubra$div_sp)
qqline(Rubra$div_sp)
hist(Rubra$div_sp) 
bn <- bestNormalize(Rubra$div_sp)
bn
plot(bn, leg_loc = "bottomright")

Rubra$div_sp <- predict(bn)

reg <- lmer(div_sp~div_endo+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lm(div_sp~div_endo, data=Rubra)
summary(reg)
anova(reg)


#Leaf water content est significatif alors on le plotte

Rubra$cat2 <-cut(Rubra$Leaf_C.N, breaks=4, right=FALSE, labels=c(1:4))
ClimatSpecies$cat3 <-cut(ClimatSpecies$Leaf_C.N, breaks=4, right=FALSE, labels=c(1:4))


ggplot() +
  geom_boxplot(data=ClimatSpecies,
               aes(x = cat3, y = div_endo, group=cat3, col="cat3"), fill = "dodgerblue",alpha=0.2,width=0.5) +
  geom_boxplot(data=Rubra,
               aes(x=cat2, y = div_endo, group = cat2, col="cat2"), fill = "gold",alpha=0.2,width=0.5) +
  
  geom_smooth(data=ClimatSpecies,
              aes(x=as.numeric(cat3), y=div_endo), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="dodgerblue", size=1.5) +
  geom_smooth(data=Rubra,
              aes(x=as.numeric(cat2), y=div_endo), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="gold", size=1.5,alpha = 0.1)+
  coord_cartesian(ylim=c(200,850)) +
  theme_classic()+
  labs(title="Impact of leaf C/N ratio on phytochemical richness",
       x="Leaf C/N ratio",
       y="Phytochemical richness (nb of molecules)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name="Leaf C/N ratio",
                   breaks=c("1", "2", "3", "4"),
                   labels=c("Small", "Medium", "Large", "Very large"))+
  labs(colour = "Plant type")+
  scale_color_manual(labels = c("F. rubra", "Festuca sp."), values = c("gold", "dodgerblue"))+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

reg <- lmer(div_sp~as.numeric(cat2)+(1|Transect.x), data = ClimatRubra)
summary(reg)
anova(reg) #Significatif

reg <- lmer(Leaf_C.N ~ div_sp + (1|Transect), data = Rubra)
summary(reg)
anova(reg) #Significatif


#############################################################################################
####################Variables edaphiques (sol)###############################################

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tmap)
library(raster)
library(lme4)
library(lmerTest)
library(bestNormalize)

#Pour Rubra##############################################################################

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Endophytes 2")

#Extraire la somme totale des endophytes comme il a été fait pour div_sp et ajouter cette colonne dans
#un tableau excel de Rubra

Rubra = read.csv("matriceendo_finalerubra.csv", sep = ";", h = T)

#Part I : SLA, WC, leaf C/N, div_sp
#(1) Normalisation des données
hist(Rubra$div_endo) 
bn <- bestNormalize(Rubra$div_endo)
bn
plot(bn, leg_loc = "bottomright")
Rubra$div_endo <- predict(bn)

hist(Rubra$CEC) 
bn <- bestNormalize(Rubra$CEC)
bn
plot(bn, leg_loc = "bottomright")
Rubra$CEC <- predict(bn)

reg <- lmer(div_endo~CEC+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lm(CEC~div_endo, data=Rubra)
summary(reg)
anova(reg)
##
hist(Rubra$PAF) 
bn <- bestNormalize(Rubra$PAF)
bn
plot(bn, leg_loc = "bottomright")
Rubra$PAF <- predict(bn)

reg <- lmer(div_endo~PAF+(1|Transect), data=Rubra)
summary(reg)
anova(reg)
##
hist(Rubra$Soil_C.N) 
bn <- bestNormalize(Rubra$Soil_C.N)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Soil_C.N <- predict(bn)

reg <- lmer(div_endo~Soil_C.N+(1|Transect), data=Rubra)
summary(reg)
anova(reg)
##
hist(Rubra$Phosphorus) 
bn <- bestNormalize(Rubra$Phosphorus)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Phosphorus <- predict(bn)

reg <- lmer(div_endo~Phosphorus+(1|Transect), data=Rubra)
summary(reg)
anova(reg)
##
hist(Rubra$Soil_Corg) 
bn <- bestNormalize(Rubra$Soil_Corg)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Soil_Corg <- predict(bn)

reg <- lmer(div_endo~Soil_Corg+(1|Transect), data=Rubra)
summary(reg)
anova(reg)
##
hist(Rubra$HR) 
bn <- bestNormalize(Rubra$HR)
bn
plot(bn, leg_loc = "bottomright")
Rubra$HR <- predict(bn)

reg <- lmer(div_endo~HR+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

Rubra$cat3 <-cut(Rubra$Phosphorus, breaks=4, right=FALSE, labels=c(1:4))
ggplot(Rubra, aes(x=as.numeric(cat3), y = div_endo))+
  geom_boxplot(aes(fill = cat3))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^1))+
  labs(title = "Impact of bioavailable phosphorus on endophytic richness",
       x="Bioavailable phosphorus (mg/g)",
       y="Endophytic richness")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Bioavailable phosphorus (mg/g)", 
                   limits=c("Low", "Mid", "High", "Very high"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

ggplot(Rubra, aes(x = Phosphorus, y = div_endo)) +
  geom_point(col = "red") +
  geom_smooth(data=Rubra,
              aes(x=Phosphorus, y=div_endo), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5)

##
hist(Rubra$div_endo) 
bn <- bestNormalize(Rubra$div_endo)
bn
plot(bn, leg_loc = "bottomright")
Rubra$div_endo <- predict(bn)

hist(Rubra$Climat) 
bn <- bestNormalize(Rubra$Climat)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Climat <- predict(bn)

reg <- lmer(div_endo~Climat+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

Rubra$cat3 <-cut(Rubra$Climat, breaks=4, right=FALSE, labels=c(1:4))
ggplot(Rubra, aes(x=as.numeric(cat3), y = div_endo))+
  geom_boxplot(aes(fill = cat3))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^1))+
  labs(title = "Impact of climate on endophytic richness",
       x="Climate (PC1)",
       y="Endophytic richness")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Climate (PC1)", 
                   limits=c("Warm+", "Warm", "Cold", "Cold+"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

ggplot(Rubra, aes(x = Climat, y = div_endo)) +
  geom_point(col = "red") +
  geom_smooth(data=Rubra,
              aes(x=Climat, y=div_endo), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5)

reg <- lmer(div_endo~as.numeric(cat3)+(1|Transect), data=Rubra)
summary(reg)
anova(reg)#Significatif

############################################################################################
###############Catégories###################################################################

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tmap)
library(raster)
library(lme4)
library(lmerTest)
library(ape)
library(plot.matrix)
library(ade4)
library(heatmaply)
library(reshape2)
library(vegan)
library(rgl)
library(phytools)
library(ggrepel)
library(bestNormalize)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Endophytes 2") 

Rubra = read.csv("matriceendo_finalerubra.csv", sep = ";", h = T)

Species = read.csv("matriceendo_finalespecies.csv", sep = ";", h = T)

#Sans catégories (points)

ggplot() +
  geom_point(data=Species,
             aes(x = Altitude, y = div_endo), col = "turquoise") +
  geom_point(data=Rubra,
             aes(x=Altitude, y = div_endo), col = "red")+
  geom_smooth(data=Species,
              aes(x=Altitude, y=div_endo), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="turquoise", size=1.5) +
  geom_smooth(data=Rubra,
              aes(x=Altitude, y=div_endo), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="red", size=1.5)+
  labs(title="Impact of elevation on endophytic richness",
       x="Elevation",
       y="Endophytic richness")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))


#Graphe au propre

Rubra$cat2 <-cut(Rubra$Altitude , breaks=4, right=FALSE, labels=c(1:4))
Species$cat3 <-cut(Species$Altitude , breaks=4, right=FALSE, labels=c(1:4))

ggplot() +
  geom_boxplot(data=Species,
               aes(x = cat3, y = div_endo, group=cat3, col = "cat3"), fill = "dodgerblue",alpha=0.2,width=0.5) +
  geom_boxplot(data=Rubra,
               aes(x=cat2, y = div_endo, group = cat2, col = "cat2"), fill = "gold",alpha=0.2,width=0.5) +
  
  geom_smooth(data=Species,
              aes(x=as.numeric(cat3), y=div_endo), method="lm",
              formula = y ~ poly(x,2), n= 40, se=FALSE, color="dodgerblue", size=1.5) +
  geom_smooth(data=Rubra,
              aes(x=as.numeric(cat2), y=div_endo), method="lm",
              formula = y ~ poly(x,1), n= 40, se=FALSE, color="gold", size=1.5,alpha = 0.1)+
  #coord_cartesian(ylim=c(200,850)) +
  theme_classic()+
  labs(title="Impact of elevation on endophytic richness",
       x="Elevation",
       y="Endophytic richness")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name="Elevation",
                   breaks=c("1", "2", "3", "4"),
                   labels=c("Low", "Mid", "High", "Very high"))+
  labs(colour = "Plant type")+
  scale_color_manual(labels = c("F. rubra", "Festuca sp."), values = c("gold", "dodgerblue"))+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

hist(Rubra$div_endo) 
bn <- bestNormalize(Rubra$HR)
bn
plot(bn, leg_loc = "bottomright")
Rubra$HR <- predict(bn)

hist(Rubra$Altitude) 
bn <- bestNormalize(Rubra$Altitude)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Altitude <- predict(bn)

reg <- lmer(div_endo~Altitude+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lmer(div_endo~as.numeric(cat2)+(1|Transect), data = Rubra)
summary(reg)
anova(reg) #Pour celui-là on ne normalise pas les données car ce sont des catégories
##
hist(Species$div_endo) 
bn <- bestNormalize(Species$div_endo)
bn
plot(bn, leg_loc = "bottomright")
Species$div_endo <- predict(bn) #Pour species on ne normalise pas car jeu de données trop petit

reg <- lm(div_endo~Altitude, data=Species)
summary(reg)
anova(reg)

reg <- lmer(div_endo~as.numeric(cat2)+(1|Transect), data = Species)
summary(reg)
anova(reg)

##############################################################################################333
###########Axe 1 PCA (PC1) structuration chimique

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tmap)
library(raster)
library(lme4)
library(lmerTest)
library(ape)
library(plot.matrix)
library(ade4)
library(heatmaply)
library(reshape2)
library(vegan)
library(rgl)
library(phytools)
library(ggrepel)
library(bestNormalize)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Endophytes 2") 

Rubra = read.csv("matriceendo_finalerubra.csv", sep = ";", h = T)

Species = read.csv("matriceendo_finalespecies.csv", sep = ";", h = T)

PC1Species$X1 <- range01(PC1Species$X1)
PC1Rubra$X1 <- range01(PC1Rubra$X1)

ggplot(Rubra, aes(x = X1, y = div_sp)) +
  geom_point(col = "red") +
  geom_smooth(data=Rubra,
              aes(x=X1, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5) +
  geom_point(data=Species,
             aes(x=X1, y=div_sp), color="turquoise", size=1.5) +
  geom_smooth(data=Species,
              aes(x=X1, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5)

#Seulement pour les espèces

Species$cat3 <-cut(Species$X1 , breaks=4, right=FALSE, labels=c(1:4))

ggplot() +
  geom_point(data=Species,
             aes(x=X1, y=div_endo), color="turquoise", size=1.5) +
  geom_smooth(data=Species,
              aes(x=X1, y=div_endo), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5)

#En boxplots

Species$cat3 <-cut(Species$X1 , breaks=4, right=FALSE, labels=c(1:4))

ggplot() +
  geom_boxplot(data=Species,
               aes(x=cat3, y=div_endo, group=cat3), color="turquoise", size=1.5) +
  geom_smooth(data=Species,
              aes(x=as.numeric(cat3), y=div_endo), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="turquoise", size=1.5)+
  theme_classic()+
  labs(title="Impact of fungal composition on endophytic richness",
       x="PC1 fungal composition",
       y="Endophytic richness")+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

reg <- lm(div_endo~X1, data=Species)
summary(reg)
anova(reg)

reg <- lm(Altitude~X1, data=Species)
summary(reg)
anova(reg)