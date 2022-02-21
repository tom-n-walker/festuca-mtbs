#Script travail de master

#Phytochimie#############################################################################
########################################################################################

#D?but du script r?alis? avec Manu

install.packages("reshape2")
install.packages("plot.matrix")
install.packages("ade4")
install.packages("phytools")

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

setwd("/Users/rasmanns/Dropbox/Lavoro/Manuscripts/2021/Nadline_festuca/analyses") 
#setwd("D:/nadline_process") 


Metabo = read.csv("nadline_full_ms1.csv", sep = ",", h = T)


####Phylogeny ####################
#setwd("D:/Affaires Nadline/Master thesis/Statistiques/Metabolomic")
#setwd("D:/nadline_process/festuca phylogeny/festuca phylogeny/new")

tree <- read.tree("ITS_FestucaSwissIQTREE.contree")

plot(tree)
### climatic niche 

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Metabolomic")
#setwd("D:/nadline_process") 


Taxon_niche = read.csv("Taxon_Predictors_Ranges.csv", sep = ",", h = T)

#fix(Metabo)

#Supprimer les colonnes row.m.z et row.retention.time------------------------

Metabo2 <- Metabo[ , ! colnames(Metabo) %in% c("row.m.z", "row.retention.time")]
dim(Metabo2)
Metabo3 <- t(Metabo2)
class(Metabo3)
dim(Metabo3)

Metabo3 <- Metabo3[-1,] #Supprimer le nom de la ligne maintenant que la matrice a ?t? transpos?e (= chang?e de sens)


row.names(Metabo3) <- gsub("X201222_PFC_metabolomics_DDA_pos_", "",row.names(Metabo3))
row.names(Metabo3) <- gsub("X201221_PFC_metabolomics_DDA_pos_", "",row.names(Metabo3))
row.names(Metabo3) <- gsub(".mzXML.Peak.height", "",row.names(Metabo3)) ## area for nadline_gnps_quant / height full ms1


#Importer tableau avec les vrais ID ------


Rubra = read.csv("FestucaRubra_Metabo2.csv", sep = ";", h = T)

row.names(Rubra) <- Rubra$ID_Metabo

row.names(Rubra) <- gsub("201222_PFC_metabolomics_DDA_pos_", "",row.names(Rubra))
row.names(Rubra) <- gsub(".mzXML Peak height", "",row.names(Rubra))


#Merger les 2 matrices --------------------

Metabo4 = merge(Rubra,Metabo3, by = "row.names")
row.names(Metabo4) <- Metabo4$Row.names


####################################################################################
####################################################################################
####################################################################################
######################## 

Metabo_species <- Metabo4 # include rubra elevation gradient et celles du JB
fix(Metabo4)
#Que JB
Metabo_species <- Metabo4[is.na(Metabo4$Altitude),] # exclude rubra r?colt?es sur un gradient altitudinal
colnames(Metabo4)
Metabo_species2 <- Metabo_species[,29:ncol(Metabo_species)] #Nom des mol?cules commence ? partir de colonne 29 incluse jusqu'? la fin
Metabo_species2[Metabo_species2 > 0] <- 1
div_species <- rowSums(Metabo_species2)

boxplot(div_species~Metabo_species$Plant.Species,las = 2)

reg = lm(div_species~Metabo_species$Plant.Species)
anova(reg)

reg = lm(div_species~Metabo_species$Elevation+Metabo_species$Plant.Species)
anova(reg)

ggplot(data = Metabo_species, aes(x= Plant.Species, y = div_species))+
  geom_boxplot(aes(fill = Elevation))+
  theme(axis.text.x = element_text(size=11, angle=90))+
  labs(title="Phytochemistry richness according to plant species",
       x="Plant Species",
       y="Phytochemistry richness")

ggplot(data = Metabo_species, aes(x= Elevation, y = div_species))+
  geom_boxplot(aes(fill = Elevation))+
  theme(axis.text.x = element_text(size=11, angle=90))+
  labs(title="Phytochemistry richness according to plant species",
       x="Plant Species",
       y="Phytochemistry richness")

Metabo_species$Elevation<-as.factor(Metabo_species$Elevation)
Metabo_species$Elevation

library(plyr)
Metabo_species$Elevation<-revalue(Metabo_species$Elevation, c("High"="Sub-alpine", "Low"="Colline", "Mid"= "Mountain", "Very high" = "Alpine"))

Metabo_species$Elevation<-factor(Metabo_species$Elevation, levels = c("Colline","Mountain", "Sub-alpine",  "Alpine"))

divsp<-ggplot(data = Metabo_species, aes(x= Elevation, y = div_species))+
  geom_boxplot(aes(fill = Elevation))+
  theme(axis.text.x = element_text(size=11, angle=90))+
  labs(subtitle = "Festuca spp.",
       x="Elevation",
       y="Phytochemical richness")+
  scale_fill_brewer(palette="Blues")+
  theme(legend.position = "none")

divsp

######Boxplot pour Rubra###################################################

Metabo_rubra <- Metabo4[!is.na(Metabo4$Altitude),] #exclus les plantes non r?colt?es sur un gradient altitudinal
Metabo_rubra2 <- Metabo_rubra[,29:ncol(Metabo_species)] #Nom des mol?cules commence ? partir de colonne 29 incluse jusqu'? la fin
Metabo_rubra2[Metabo_rubra2 > 0] <- 1
div_rubra <- rowSums(Metabo_rubra2)

boxplot(div_rubra~Metabo_rubra$Altitude.2,las = 2)
cor.test(Metabo_rubra$Altitude.2,div_rubra)

fit = lm(div_rubra~Altitude.2, data=Metabo_rubra)
anova(fit)

  



Metabo_rubra$Elevation<-as.factor(Metabo_rubra$Elevation)
Metabo_rubra$Elevation<-revalue(Metabo_rubra$Elevation, c("High"="Sub-alpine", "Low"="Colline", "Mid"= "Mountain", "Very High" = "Alpine"))

Metabo_rubra$Elevation<-factor(Metabo_rubra$Elevation, levels = c("Colline","Mountain", "Sub-alpine",  "Alpine"))

divr<-ggplot(data = Metabo_rubra, aes(x= Elevation, y = div_rubra))+
  geom_boxplot(aes(fill = Elevation))+
  theme(axis.text.x = element_text(size=11, angle=90))+
  labs(subtitle = "Festuca rubra",
       x="Elevation",
       y="Phytochemical richness")+
  scale_fill_brewer(palette="Blues")+
  theme(legend.position = "none")

####
library(ggpubr)

ggarrange(divsp +theme(plot.subtitle=element_text(size=12,  face="italic", color="black"))+
            theme(axis.title.x=element_blank()),
          divr +theme(plot.subtitle=element_text(size=12,  face="italic", color="black")),
            #theme(axis.title.x=element_blank(),
            #      axis.text.x=element_blank()),
          nrow = 2,labels = c("A","B"))


###############JB BRAY Distance #############################################################


Metabo_species <- Metabo4[is.na(Metabo4$Altitude),] # exclude rubra gelevation radient
Metabo_species <- Metabo4[!is.na(Metabo4$Altitude),] #exclus les plantes non r?colt?es sur un gradient altitudinal
Metabo_species <- Metabo4 # include rubra gelevation radient
Metabo_species2 <- Metabo_species[,29:ncol(Metabo_species)]
Metabo_species2[Metabo_species2 > 0] <- 1
div_species <- rowSums(Metabo_species2)

boxplot(div_species~Metabo_species$Plant.Species,las = 2)
cor.test(Metabo_rubra$Altitude,div_rubra)

Metabo_species2_env <- Metabo_species[,3:28] #on garde les colonnes 3 ? 28 pour cr?er la matrice environnementale
Metabo_species2_env$div_sp <- div_species #Rajouter la colonne div_sp qui contient les ?l?ments de div_species (voir au dessus)
Metabo_species2_env$taxon <- gsub(" ","_",Metabo_species2_env$Plant.Species) #gsub(pattern = "o" , replacement = "O") donc ici
#pour remplacer l'espace entre les noms par un tiret
colnames(Metabo_species2_env)

Climat <- data.frame(Metabo_species2_env[,c(1, 6, 27)])

#### add pcao metric 

#calculer une matrice de distance en m?thode Bray curtis ----
#Bray curtis = structures chimiques proches, pas richesse !!
Metabo_species2 <- Metabo_species[,29:ncol(Metabo_species)]

spe.bray <- vegdist(Metabo_species2) #Calculer la distance de Bray entre les mol?cules et on verra si corr?l? aussi avec
#variables environnementales


#Multidimensional scaling (MDS) pour projeter la matrice de distance
#en plusieurs axes ------------------------------------------

spe.b.pcoa <- cmdscale(spe.bray, k=3) #k=3 car on veut une PCOA a 3 axes (dimensions)
spe.b.pcoa <- data.frame(spe.b.pcoa) #Mise en tableau des donn?es brutes (donc en 3 axes distincts qu'on peut voir dans
#"l'environnement" (=panneau ? droite de R))


Metabo_full_sp = merge(Metabo_species2_env,spe.b.pcoa, by = "row.names")
#Regarder s'il y a un effet de la matrice environnementale sur cette distance, on rajoute les coordonn?es (les 3 axes) aux
#variables environnementales. On ajoute des coordonn?es (de Bray donc de distance) ? mes f?tuques.
plot(Metabo_full_sp$X1,Metabo_full_sp$X2,col=as.factor(Metabo_full_sp$Plant.Species),pch=16,cex=2)

plot(Metabo_full_sp$X1,Metabo_full_sp$div_sp)
#PC1Rubra <- data.frame(Metabo_full_sp$X1, Metabo_full_sp$div_sp, Metabo_full_sp$PlantID)
#write.table(PC1Rubra, "PC1Rubra.csv", sep = ";")

##########Rubra BRAY Distance #############################################################################################

Metabo_rubra2_env <- Metabo_rubra[,3:28]
Metabo_rubra2_env$div_sp <- div_rubra


Metabo_rubra2 <- Metabo_rubra[,29:ncol(Metabo_rubra)]

spe.bray <- vegdist(Metabo_rubra2)
spe.b.pcoa <- cmdscale(spe.bray, k=3)
spe.b.pcoa <- data.frame(spe.b.pcoa)
Metabo_full_sp2 = merge(Metabo_rubra2_env,spe.b.pcoa, by = "row.names")
plot(Metabo_full_sp2$X1,Metabo_full_sp2$X3,col=as.factor(Metabo_full_sp$Altitude.2),pch=16,cex=2)
text(Metabo_full_sp2$X1,Metabo_full_sp2$X3, labels = Metabo_full_sp$Altitude.2,cex=0.9)

############################################################################33

###Phylogeny

Plant_Species <- unique(Metabo_full_sp$Plant.Species)
Plant_Species <- gsub(" ", "_",Plant_Species)
Plant_Species <- gsub("Festuca_valesciaca", "Festuca_valesiaca",Plant_Species)


##### species niche #S?lectionner les Festuca sp. dans le tableau complet de "Taxon_niche" et creation d'un data frame
matt_niche = NULL
for ( i in c (1:length(Plant_Species))) {
  
  matt_inter <- Taxon_niche[grep(Plant_Species[i],Taxon_niche$taxon),]
  matt_inter <-  matt_inter[1,]
  matt_inter$taxon <- Plant_Species[i]
  matt_niche <-  rbind(matt_niche,matt_inter)
}

matt_niche_clim_sp <- matt_niche[,c(6,14,18,21,22)]
#matt_niche_clim_sp <- matt_niche[,c(6, 9, 10, 11, 12, 14,18,21,22)] #contient le ndvi
#write.table(matt_niche_clim_sp, "ndvispecies.csv", sep = ";")
#write.table(matt_niche_clim_sp, "climspecies.csv", sep = ";")
row.names(matt_niche_clim_sp) <- matt_niche$taxon #  DEM_Mean precyy_Mean  sradyy_Mean taveyy_Mean   taveyy_Q05

########################################################################################
#PCA clim
library(ade4)


dudi1 <- dudi.pca(matt_niche_clim_sp, scale = TRUE, scan = FALSE, nf = 2)

#prendre la PC1 (axe x) pour les sites (1?re colonne de dudi$li donc axe 1 = PC1)
climat_sp <- dudi1$li[,1]*-1
taxon <- row.names(dudi1$li)
matt_clim <- data.frame(taxon, climat_sp)
matt_niche_clim_sp2 <- matt_niche[,c(3,6,14,18,21,22)] 
matt_niche = merge(matt_niche_clim_sp2,matt_clim, by = "taxon") 

Metabo_full_sp1 = merge(Metabo_full_sp,matt_niche, by = "taxon") 

#ClimatNouveau <- data.frame(Metabo_full_sp1[,c(1, 3, 29, 38)])

########################################################################################

Metabo_stat <- Metabo_full_sp1[!is.na(Metabo_full_sp1$Altitude),] ## delet ! to exclude rubra field data  # ! to keep rubra
#S?lectionner tous ceux qui n'ont pas de "NA" dans la colonne "Altitude" car ! signifie "diff?rent de"

plot(Metabo_stat$Altitude,Metabo_stat$div_sp)
cor.test(Metabo_stat$Altitude,Metabo_stat$div_sp)


sp_div_mean <- tapply(Metabo_stat$div_sp ,Metabo_stat$Plant.Species, mean) ## calculer la moyenne de phytodiversit? par esp?ces
names(sp_div_mean) <- gsub(" ", "_",names(sp_div_mean))
names(sp_div_mean) <- gsub("Festuca_valesciaca", "Festuca_valesiaca",names(sp_div_mean))


#Mettre de la couleur et du texte:

myColorRamp <- function(colors, values) { 
  v <- (values - min(values))/diff(range(values)) 
  x <- colorRamp(colors)(v) 
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255) 
} 


cols <- myColorRamp(c("darkred","orange","darkgreen"),Metabo_stat$Altitude) ### colour en fonction elevation

plot(Metabo_stat$Altitude,Metabo_stat$div_sp,col=cols,pch=16,cex=3)
text(Metabo_stat$Altitude,Metabo_stat$div_sp,labels = Metabo_stat$Site)

####################################################################
#####Div_sp en fonction des altitudes pour chaque transect##########

Metabo_rubra_env_plot <- Metabo_stat[grep("Zermatt", Metabo_stat$Transect),] #  "Chasseron" "Chasseral" "Zermatt"   "Salgesch"  "Lavey"

plot(Metabo_rubra_env_plot$Altitude,Metabo_rubra_env_plot$div_sp, pch=16,cex=3)
text(Metabo_rubra_env_plot$Altitude,Metabo_rubra_env_plot$div_sp,labels = Metabo_rubra_env_plot$Site) 


#JB--div_sp en fonction de Optimum de niche---------------------------------------------------------------

Metabo_stat <- Metabo_full_sp1[is.na(Metabo_full_sp1$Altitude),] ## delet ! to exclude rubra field data  # ! to keep rubra
#S?lectionner tous ceux qui ont des "NA" dans la colonne "Altitude"

plot(Metabo_stat$DEM_Mean,Metabo_stat$div_sp)
cor.test(Metabo_stat$DEM_Mean,Metabo_stat$div_sp)

##Matrice de distance en 3D###########################################

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#open3d(box=FALSE, axes=TRUE)
#plot3d(spe.b.pcoa[,1],spe.b.pcoa[,2], spe.b.pcoa[,3], type="p", xlab ="PC0A2",ylab ="PC0A1",zlab ="" ,box=FALSE, axes=TRUE,zlim=c(0,0),size=3)#,alpha=0.1) 
#
# a simple white background
bg3d("white")


sizer <- (range01(div_species)+5)

for (i in c(1:nrow(Metabo_stat))) {
  
  rgl.points(Metabo_stat[i,30:32],col="black",size=sizer[i],box=FALSE, axes=TRUE) #pris les colonnes X1, X2, X3 pour la PCA
  
}



cols <- myColorRamp(c("darkred","orange","darkgreen"),Metabo_stat$Altitude) ### colour

rgl.texts(Metabo_stat[,30:32],text=Metabo_stat$Plant.Species,col=cols,cex= 0.9)#,color=as.numeric(data_pool_plot$runer))


#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
### full species 

############## phylogeny ####################################################################
#Esp?ces dans l'ordre de leur optimum de niche avec la moyenne de leur diversit? phytochimique.

Metabo_stat <- Metabo_full_sp1[is.na(Metabo_full_sp1$Altitude),] ## delet ! to exclude rubra field data  # ! to keep rubra


ggplot(Metabo_stat, aes(x=reorder(DEM_Mean, DEM_Mean), y = div_sp )) + 
  geom_boxplot(aes(x=reorder(DEM_Mean, DEM_Mean)))+ scale_color_hue(l=40, c=35)+
  theme(axis.text.x = element_text(angle=70,hjust=1,size=10),panel.background = element_rect(fill = "white", colour = "grey50"),
        plot.margin = unit(c(0.5,3.2,0.2,0.2), "cm")) + 
  labs(y = "Chemical family consensus",x = "elevation")+
  stat_smooth(method="lm", formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", aes(group=1), size=1.5)
#stat_smooth(n= 40, se=TRUE, color="red", aes(group=1), size=1.5)


#############################################################################################
#Ploter les esp?ces en fonction de leur optimum de niche (points noirs) vs ploter les F.rubra du gradient (points rouges)
Metabo_stat_full_sp <- Metabo_full_sp1[is.na(Metabo_full_sp1$Altitude),]
Metabo_stat_rubra_sp <- Metabo_full_sp1[!is.na(Metabo_full_sp1$Altitude),]
#Metabo_stat_rubra_sp <- Metabo_stat_rubra_sp[grep("Zermatt", Metabo_stat_rubra_sp$Transect),] #  "Chasseron" "Chasseral" "Zermatt"   "Salgesch"  "Lavey"

ggplot(Metabo_stat, aes(x= DEM_Mean, y = div_sp))+
  geom_point()+
  stat_smooth( method="lm", formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", aes(group=1), size=1.5) 


plot(Metabo_stat$DEM_Mean,Metabo_stat$div_sp,pch=16,col=1)
Metabo_stat_rubra_sp <- Metabo_full_sp[!is.na(Metabo_full_sp$Altitude),]
points(Metabo_stat_rubra_sp$Altitude,Metabo_stat_rubra_sp$div_sp,pch=16,col="red")


#########################################################################################################
############## phylogeny ################################################################################

#creation arbre phylo
#Parmi l'arbre phylog?n?tique en entier, on s?lectionne seulement les noms des esp?ces qui apparaissent
#dans notre tableau sous "Plant_Species"
tip.labels=tree$tip.label
sub.tips<-setdiff(tip.labels, Plant_Species)
new_tree<-drop.tip(tree, sub.tips)
tr <- root(new_tree, 1)
plot(tr)


### phylo signal chimie richness (mis ensemble en fonction de la diversit? compos?s) // ARBRE PHYLOGENETIQUE
sp_div_mean <- tapply(Metabo_stat$div_sp ,Metabo_stat$Plant.Species, mean) ## calculer la moyenne de phytodiversit? par esp?ces
names(sp_div_mean) <- gsub(" ", "_",names(sp_div_mean))
names(sp_div_mean) <- gsub("Festuca_valesciaca", "Festuca_valesiaca",names(sp_div_mean))

obj<-contMap(tr,sp_div_mean,plot=FALSE)
obj<-setMap(obj,invert=TRUE)

plot(obj,type="fan",outline=FALSE)
phylosig(tr,sp_div_mean,method="lambda",test=TRUE) #YEY!
#+ c'est rouge, + il y a de compos?s.
#Ici signal phylog?n?tique = structur? par rapport ? phylog?nie.


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
#Signifie probablement : Les points proches sont proches phylog?n?tiquement. Ici on voit que les
#points proches n'ont pas forc?ment la m?me niche (car couleur correspond ? l'?l?vation de la niche)
#ni forc?ment la m?me richesse mol?culaire (car pas ? la m?me hauteur sur l'axe sp_div_mean).
#C'est donc pour voir si : Si structurellement proches, phylog?n?tiquement aussi ?
#Ici on voit que c'est surtout l'axe 2 (X2) qui montre un motif malgr? une partie bruit mais dans une
#partie il y a un signal donc c'est pas la phylog?nie qui structure le + (sinon il y aurait un motif
#aussi avec l'axe x = X1) mais cela explique beaucoup de sous-parties.


### phylo signal chimie structure (arbre en fonction de la structure des compos?s)-----------------------

sp_struc <- tapply(Metabo_stat_full_sp$X1, Metabo_stat_full_sp$Plant.Species, mean) ## calculer la moyenne de structure chimique par esp?ces
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

dist_chimie <- Metabo_species[,29:ncol(Metabo_species)]
dist_chimie <- aggregate(dist_chimie,by=list(Metabo_species$Plant.Species),FUN=mean)
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

####Heatmap###################################################################################

phylo_Matrix = as.matrix(phylo.dist) #Transformer une variable en matrice
chimie_Matrix = as.matrix(dist_chimie)

heatmap(phylo_Matrix, chimie_Matrix)#mantel.test pour faire des corr?lations sur des distances

#graphes-----------------------------------------------

tanglegram(dend12,common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
           margin_inner=9,
           lwd=1,remove_nodePar = TRUE,rank_branches = TRUE,type = "triangle")


# dend12 %>% untangle %>% tanglegram
dend12 %>%  untangle(method = "ladderize") %>% tanglegram(common_subtrees_color_lines = TRUE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
                                                          margin_inner=9,
                                                          lwd=0.1,remove_nodePar = TRUE,rank_branches = TRUE,axes = FALSE,lab.cex = 1)

#Test for coloration but doesn't work

dend12 %>%  untangle(method = "ladderize") %>% tanglegram(common_subtrees_color_lines = TRUE, color_lines = Species$DEM_Mean, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
                                                          margin_inner=9,
                                                          lwd=0.1,remove_nodePar = TRUE,rank_branches = TRUE,axes = FALSE,lab.cex = 1)

setwd("D:/Affaires Nadline/Master thesis/Statistiques/R?f?rences") 

Species = read.csv("Tableau Festuca Complet.csv", sep = ";", h = T)

Species$cat <-cut(Species$DEM_Mean, breaks=4, right=FALSE, labels=c(1:4))


### phylo signal optimum de niche------------------------------------------------------------------------

opti_niche <- matt_niche$DEM_Mean
names(opti_niche) <- matt_niche$taxon

obj<-contMap(tr,opti_niche,plot=FALSE)
obj<-setMap(obj,invert=TRUE)


plot(obj,type="fan",outline=FALSE)
phylosig(tr,opti_niche,method="lambda",test=TRUE)

#Un peu structur? par optiumum de niche (= climat o? on a + de chance de trouver l'esp?ce) mais signal
#moins fort.

##correlation diversit? et structure des compos?s chimiques par esp?ce et optimum de niche

opti_niche2 <- opti_niche[order(names(opti_niche))]
sp_div2 <- sp_div_mean[order(names(sp_div_mean))]
sp_struc2 <- sp_struc[order(names(sp_struc))]

cor.test(opti_niche2,sp_div2)
cor.test(opti_niche2,sp_struc2) #TOOOOUT petit peu
cor.test(sp_div2, sp_struc2)


matt_mean_sp <- data.frame(opti_niche2,sp_div2)
ggplot() +
  geom_point(data=matt_mean_sp,
             aes(x=opti_niche2, y=sp_div2)) +
  geom_smooth(data=matt_mean_sp, 
              aes(x=opti_niche2, y=sp_div2), method="lm", 
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5)


#(2)
##Statiques sur Altitude et phytochimie###################################################
###########################################################################################

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

setwd("D:/Affaires Nadline/Master thesis/Statistiques/R?f?rences") 

Rubra = read.csv("Tableau Rubra.csv", sep = ";", h = T) ### merge data clim

Species = read.csv("Tableau Festuca Complet.csv", sep = ";", h = T)### merge data clim

#Pour F. rubra##############
model_null <- lm(Rubra$div_sp ~ 1)

Rubra$cat2 <-cut(Rubra$Altitude , breaks=4, right=FALSE, labels=c(1:4))
#A runer avant la normalisation car pas de normalisation sur
#les cat?gories
#Donc on fait toutes les stats avec cat 2 car we normalize the 
#data to bring all the variables to the same range.
#Ici on a une variable cat?gorique (cat) alors l'autre (div_sp)
#On ne le normalise pas non plus

reg <- lmer(div_sp~as.numeric(cat2)+(1|Transect), data = Rubra)
summary(reg)
anova(reg) #categorical variable numeric

anova(lm(Rubra$div_sp~Rubra$cat2)) #categorical variable

reg <- lmer(div_sp~cat2+(1|Transect), data = Rubra)
summary(reg)
anova(reg) #categorical variable + random factor
#Pour l'AIC:
deg1 <- lmer(div_sp~as.numeric(cat2)+(1|Transect), data = Rubra)
deg4 <- lm(Rubra$div_sp~Rubra$cat2)
deg5 <- lmer(div_sp~cat2+(1|Transect), data = Rubra)
#Ensuite on fait sur la variable continue
hist(Rubra$div_sp) 
bn <- bestNormalize(Rubra$div_sp)
bn
plot(bn, leg_loc = "bottomright")
Rubra$div_sp <- predict(bn)

hist(Rubra$Altitude) 
bn <- bestNormalize(Rubra$Altitude)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Altitude <- predict(bn)

#response~explanatory
anova(lm(Rubra$div_sp~Rubra$Altitude)) #continuous variable

reg <- lmer(div_sp~Altitude+(1|Transect), data = Rubra)
summary(reg)
anova(reg) #continuous variable + random factor

deg2 <- lm(Rubra$div_sp~Rubra$Altitude) #Meilleur mod?le selon AIC plus bas
deg3 <- lmer(div_sp~Altitude+(1|Transect), data = Rubra)
#AIC

AIC(model_null)
AIC(deg1)
AIC(deg2) #Meilleur mod?le
AIC(deg3)
AIC(deg4)
AIC(deg5)
#Pour les species###############
model_null <- lm(Species$div_sp ~ 1)
deg2 <- lm(Species$div_sp ~ poly(as.numeric(Species$cat3), 2, raw=TRUE))
deg1 <- lm(Species$div_sp ~ poly(as.numeric(Species$cat3), 1, raw=TRUE))
deg3 <- lm(div_sp ~ poly(DEM_Mean, 2), data = Species)
deg4 <- lm(div_sp ~ poly(DEM_Mean, 1), data = Species)
AIC(model_null)
AIC(deg1)
AIC(deg2)
AIC(deg3)
AIC(deg4)
anova(deg2)
anova(deg1)
anova(deg3)
anova(deg4)


#(3)
############Cat?gories photochimie##################################################
######################################################################################

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tmap)
library(raster)
library(lme4)
library(lmerTest)
library(bestNormalize)

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Corr?lations Rubra")

Rubra = read.csv("FestucaRubra_SLA_WC_Chemicals.csv", sep = ";", h = T)


#Normalisation des donn?es#################################################
#(1) Climat

#Check for normality using (I) qqplots and (II) histogram
qqnorm(Rubra$Climat)
qqline(Rubra$Climat) #qqplots

hist(Rubra$Climat) #Histogram
#bn <- bestNormalize(Rubra$Climat, allow_orderNorm = FALSE, out_of_sample = FALSE)
#Perform the best normalize function and r tells you which transformation would be the best
bn <- bestNormalize(Rubra$Climat)
bn
#Plot the possible transformations
plot(bn, leg_loc = "bottomright")

#Proceed to the transformation
Rubra$Climat <- predict(bn)

#Check if it is better now:
qqnorm(Rubra$Climat)
qqline(Rubra$Climat) 
hist(Rubra$Climat)

#(2) Altitude

library(bestNormalize)

qqnorm(Rubra$Altitude)
qqline(Rubra$Altitude) #Consid?r? comme distribution normale car sym?trique 

hist(Rubra$Altitude)

bn <- bestNormalize(Rubra$Altitude)
bn
#Plot the possible transformations
plot(bn, leg_loc = "bottomright")

#Proceed to the transformation
Rubra$Altitude <- predict(bn)

#(3) SLA

qqnorm(Rubra$SLA)
qqline(Rubra$SLA)
hist(Rubra$SLA)

bn <- bestNormalize(Rubra$SLA)
bn
plot(bn, leg_loc = "bottomright")
Rubra$SLA <- predict(bn)

#(4) Water Content

qqnorm(Rubra$Water.Content)
qqline(Rubra$Water.Content) #Normale distribution

#(5) div_sp

qqnorm(Rubra$div_sp)
qqline(Rubra$div_sp)
hist(Rubra$div_sp)

bn <- bestNormalize(Rubra$div_sp)
bn
plot(bn, leg_loc = "bottomright")
Rubra$div_sp <- predict(bn)

ggplot(Rubra, aes(y=Climat, x = Elevation))+
  geom_boxplot(aes(fill=Elevation))+
  geom_point()+
  labs(title="Climatic conditions at collection site according to elevation", y = "Climate (PCA Axis 1)")+
  scale_x_discrete(name ="Altitudinal gradient", 
                   limits=c("Low","Mid","High", "Very High"))+
  theme(plot.title = element_text(hjust = 0.5))

#Plot climat avec les bonnes couleur
Rubra$cat3 <-cut(Rubra$Altitude , breaks=4, right=FALSE, labels=c(1:4))
ggplot(Rubra, aes(y=Climat, x = cat3))+
  geom_boxplot(aes(fill=cat3))+
  scale_x_discrete(breaks=c("1","2","3", "4"),
                   labels=c("Low", "Mid", "High", "Very high"))+
  labs(title="Impact of elevation on climate", y = "Climate (PCA Axis 1)", x = "Altitudinal gradient")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_discrete(name = "Elevation", labels = c("Low", "Mid", "High", "Very high"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))


reg<-lmer(Climat~Altitude+(1|Transect), data = Rubra)
summary(reg)
anova(reg)
#Corr?lations entre continuous variables###############################################
#Test for a simple linear regression:
ggplot(Rubra, aes(x=Climat, y = Altitude))+
  geom_point(aes(col = Site))+
  geom_smooth(method = "lm")+
  labs(title="Climate according to elevation")

cor.test(Rubra$Climat, Rubra$Altitude, method="pearson")

ggplot(Rubra, aes(x=SLA, y = Water.Content))+
  geom_point(aes(col = Site))+
  geom_smooth(method = "lm")+
  labs(title="SLA according to WC")

cor.test(Rubra$SLA, Rubra$Water.Content, method="pearson")

#Regression models avec Altitude############################################################
#(1) SLA
WC_Rubra <- ggplot(Rubra, aes(x=Altitude, y = SLA))+
  geom_point(aes(col=Transect))+
  geom_smooth(method = "lm", formula=y~x+I(x^2))+
  labs(title="Water Content according to transect and elevation")
WC_Rubra

ggplot(Rubra, aes(x=Altitude.2, y = SLA))+
  geom_point(aes(col = Site))+
  geom_smooth(method = "lm")+
  labs(title="SLA according to WC")

ggplot(Rubra, aes(x = Altitude.2, y = SLA))+
  labs(title = "Specific leaf area (SLA) according to sites and altitude",
       x="Altitude (m)",
       y="SLA (m^2/kg)")+
  geom_boxplot(aes(x = Altitude.2, y = SLA, fill = Site))+
  scale_fill_manual(breaks = c("Chasseral1", "Chasseral2", "Chasseral3",
                               "Chasseron1", "Chasseron2", "Chasseron3",
                               "Lavey1", "Lavey2", "Lavey3", "Lavey4",
                               "Salgesch1", "Salgesch2", "Salgesch3", "Salgesch4",
                               "Zermatt1", "Zermatt2", "Zermatt3", "Zermatt4"),
                    values = c("red", "red", "red",
                               "orange", "orange", "orange",
                               "Yellow", "Yellow", "Yellow", "Yellow",
                               "green", "green", "green", "green",
                               "blue", "blue", "blue", "blue"))+
  geom_point()+
  facet_grid(.~Transect)+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2))

reg <- lmer(SLA~Altitude+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

cor.test(Rubra$SLA, Rubra$Altitude, method = "pearson")

#(2) Water content

Rubra$cat3 <-cut(Rubra$Altitude , breaks=4, right=FALSE, labels=c(1:4))

ggplot(Rubra, aes(x=as.numeric(cat3), y = Water.Content))+
  geom_boxplot(aes(fill = cat3, width=0.5))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^1))+
  labs(title = "Impact of elevation on leaf water content",
       x="Altitude (m)",
       y="Leaf water content (%)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Elevation", 
                   limits=c("Low", "Mid", "High", "Very high"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

reg <- lmer(Water.Content~Altitude+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

#(3) Div_sp
Rubra$cat3 <-cut(Rubra$Altitude , breaks=4, right=FALSE, labels=c(1:4))
ggplot(Rubra, aes(x=as.numeric(cat3), y = div_sp))+
  geom_boxplot(aes(fill = cat3))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^2))+
  labs(title="Phytochemistry richness according to elevation",
       x="Elevation",
       y="Phytochemistry richness (nb of molecules)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Elevation", 
                   limits=c("Low", "Mid", "High", "Very high"))


ggplot(Rubra, aes(x=as.numeric(cat3), y = div_sp))+
  geom_boxplot(aes(col = cat3))+
  geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^2))+
  labs(title="SLA according to WC")

ggplot(Rubra, aes(x = Altitude.2, y = div_sp))+
  labs(title = "Phytochemistry richness according to altitude at each site",
       x="Altitude (m)",
       y="Phytochemistry richness (nb of molecules)")+
  geom_boxplot(aes(x = Altitude.2, y = div_sp, fill = Site))+
  scale_fill_manual(breaks = c("Chasseral1", "Chasseral2", "Chasseral3",
                               "Chasseron1", "Chasseron2", "Chasseron 3",
                               "Lavey1", "Lavey2", "Lavey3", "Lavey4",
                               "Salgesch1", "Salgesch2", "Salgesch3", "Salgesch4",
                               "Zermatt1", "Zermatt2", "Zermatt3", "Zermatt4"),
                    values = c("red", "red", "red",
                               "orange", "orange", "orange",
                               "Yellow", "Yellow", "Yellow", "Yellow",
                               "green", "green", "green", "green",
                               "blue", "blue", "blue", "blue"))+
  geom_point()+
  facet_grid(.~Transect)+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2))


reg=lm(div_sp~Altitude.2*Transect, data=Rubra)
summary(reg)
anova(reg)

cor.test(Rubra$div_sp, Rubra$Altitude)

#(4) Soils Normalization of data
#4.1 CEC

qqnorm(Rubra$CEC)
qqline(Rubra$CEC)
hist(Rubra$CEC)
bn <- bestNormalize(Rubra$CEC)
bn
plot(bn, leg_loc = "bottomright")
Rubra$CEC <- predict(bn)

reg <- lmer(CEC~Altitude+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <-lm(CEC~Altitude, data=Rubra)
summary(reg)
anova(reg)

reg <- lmer(div_sp~CEC+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

qqnorm(Rubra$PAF)
qqline(Rubra$PAF)
hist(Rubra$PAF)
bn <- bestNormalize(Rubra$PAF)
bn
plot(bn, leg_loc = "bottomright")
Rubra$PAF <- predict(bn)

reg <- lmer(PAF~Altitude+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lmer(div_sp~PAF+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lmer(PAF~div_sp+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

qqnorm(Rubra$Leaf_C.N)
qqline(Rubra$Leaf_C.N)
hist(Rubra$Leaf_C.N)
bn <- bestNormalize(Rubra$Leaf_C.N)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Leaf_C.N <- predict(bn)

reg <- lmer(Leaf_C.N~Altitude+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

qqnorm(Rubra$Soil_C.N)
qqline(Rubra$Soil_C.N)
hist(Rubra$Soil_C.N)
bn <- bestNormalize(Rubra$Soil_C.N)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Soil_C.N <- predict(bn)

reg <- lmer(Soil_C.N~Altitude+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

cor.test(Rubra$Soil_C.N, Rubra$Altitude)

reg <- lmer(div_sp~Soil_C.N+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

Rubra$cat3 <-cut(Rubra$Soil_C.N, breaks=9, right=FALSE, labels=c(1:9))
coefs <- coef(lm(div_sp ~ cat3, data = Rubra))
coefs
ggplot(Rubra, aes(x = cat3, y = div_sp)) +
  geom_boxplot(col = "red") +
  geom_smooth(data=Rubra,
              aes(x=cat3, y=div_sp), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="red", size=1.5)+
  labs(title="Impact of soil C/N ratio on phytochemical richness",
       x="Soil C/N ratio",
       y="Phytochemical richness (nb of molecules)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(breaks=c("1","2","3", "4", "9"),
                   labels=c("Low+", "Low", "Mid", "High", "High+"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_abline(intercept = 554.46154, slope = -16, color = "red")+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))
#geom_abline(intercept = coefs[1], slope = coefs[2])

ggplot(Rubra, aes(x=Soil_C.N, y = div_sp))+
  geom_point()


qqnorm(Rubra$Phosphorus)
qqline(Rubra$Phosphorus)
hist(Rubra$Phosphorus)
bn <- bestNormalize(Rubra$Phosphorus)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Phosphorus <- predict(bn)

reg <- lmer(Phosphorus~Altitude+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lmer(div_sp~Phosphorus+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

qqnorm(Rubra$HR)
qqline(Rubra$HR)
hist(Rubra$HR)
bn <- bestNormalize(Rubra$HR)
bn
plot(bn, leg_loc = "bottomright")
Rubra$HR <- predict(bn)

reg <- lmer(HR~Altitude+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lmer(div_sp~HR+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

qqnorm(Rubra$Soil_Corg)
qqline(Rubra$Soil_Corg)
hist(Rubra$Soil_Corg)
bn <- bestNormalize(Rubra$Soil_Corg)
bn
plot(bn, leg_loc = "bottomright")
Rubra$Soil_Corg <- predict(bn)

reg <- lmer(Soil_Corg~Altitude+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

reg <- lmer(div_sp~Soil_Corg+(1|Transect), data=Rubra)
summary(reg)
anova(reg)

Rubra$cat3 <-cut(Rubra$Altitude, breaks=4, right=FALSE, labels=c(1:4))

ggplot(Rubra, aes(x=as.numeric(cat3), y = CEC))+
  geom_boxplot(aes(fill = cat3))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^2))+
  labs(title="Cationic exchange capacity (CEC) according to elevation",
       x="Altitude (m)",
       y="CEC (cmolc.kg-1)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Elevation", 
                   limits=c("Low", "Mid", "High", "Very high"))


ggplot(Rubra, aes(x=as.numeric(cat3), y = Leaf_C.N))+
  geom_boxplot(aes(fill = cat3))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^2))+
  labs(title = "Impact of elevation on leaf C/N ratio",
       x="Elevation",
       y="Leaf C/N ratio")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Elevation", 
                   limits=c("Low", "Mid", "High", "Very high"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))


ggplot(Rubra, aes(x=as.numeric(cat3), y = Phosphorus))+
  geom_boxplot(aes(fill = cat3))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^2))+
  labs(title = "Impact of elevation on bioavailable phosphorus",
       x="Altitude (m)",
       y="Bioavailable phosphorus (mg/g)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Elevation", 
                   limits=c("Low", "Mid", "High", "Very high"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

ggplot(Rubra, aes(x=as.numeric(cat3), y = Phosphorus))+
  geom_boxplot(aes(fill = cat3))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^2))+
  labs(title = "Bioavailble phosphorus according to transect and altitude",
       x="Altitude (m)",
       y="Bioavailable phosphorus (mg/g)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Elevation", 
                   limits=c("Low", "Mid", "High", "Very high"))


#PCA#################################################################

#(3) PCA des traits fonctionnels (SLA et WC) et du climat en fonction de l'altitude et par cat?gories
#d'?l?vation #########################################################################################

rm(list=ls())
library(tmap)
library(raster)
library(lme4)
library(lmerTest)
library(bestNormalize)

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Corr?lations Rubra")

Rubra = read.csv("FestucaRubra_SLA_WC_Chemicals.csv", sep = ";", h = T)

library(ggplot2)
library(dplyr)
library(vegan)
library(remotes)
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(calibrate)
library(ggrepel)


#creation d'un data.frame
Rubra2 = as.data.frame(cbind(Climat = Rubra$Climat, WC = Rubra$Water.Content, SLA = Rubra$SLA, Elevation = Rubra$Elevation, Altitude = Rubra$Altitude))

#conversion en numeric variables
Rubra2$Climat = as.numeric(Rubra2$Climat)
Rubra2$WC = as.numeric(Rubra2$WC)
Rubra2$SLA = as.numeric(Rubra2$SLA)
Rubra2$Altitude = as.numeric(Rubra2$Altitude)

#normaliser les donnees
rescale = function(Rubra2) (Rubra2-min(Rubra2))/(max(Rubra2)-min(Rubra2)) * 100
Rubra2 = rescale(Rubra2)
###adonis
fit = adonis(Rubra2[,1:3]~Altitude, distance="camberra", permutations=999, data=Rubra2) #je sais pas ce que c'est...
fit
#p-value = 0.034

#PCA avec transect, site, SLA, WC, Climat
Rubra.pca = prcomp(Rubra[,c(8, 19, 22)], center = TRUE, scale. = TRUE)
ggbiplot(Rubra.pca, obs.scale = 1, var.scale = 1, groups = Rubra$Elevation, labels = Rubra$Transect, ellipse = TRUE, circle = FALSE)

#(4) PCA des traits fonctionnels (SLA et WC), climat et div_sp en fonction de l'altitude et par cat?gories
#d'?l?vation#######################################################################################

Rubra3 = as.data.frame(cbind(Phytochemistry = Rubra$div_sp, WC = Rubra$Water.Content, SLA = Rubra$SLA, Leaf_C.N = Rubra$Leaf_C.N, Elevation = Rubra$Elevation, Altitude = Rubra$Altitude))

#Rubra3$Climat = as.numeric(Rubra3$Climat)
Rubra3$WC = as.numeric(Rubra3$WC)
Rubra3$SLA = as.numeric(Rubra3$SLA)
Rubra3$Phytochemistry = as.numeric(Rubra3$Phytochemistry)
Rubra3$Altitude = as.numeric(Rubra3$Altitude)
Rubra3$Leaf_C.N = as.numeric(Rubra3$Leaf_C.N)

rescale = function(Rubra3) (Rubra3-min(Rubra3))/(max(Rubra3)-min(Rubra3)) * 100
Rubra2 = rescale(Rubra3)

fit = adonis(Rubra3[,1:4]~Elevation, distance="camberra", permutations=999, data=Rubra3) #je sais pas ce que c'est...
fit
Rubra.pca = prcomp(Rubra[,c(11, 19, 22, 26)], center = TRUE, scale. = TRUE)
ggbiplot(Rubra.pca, obs.scale = 1, var.scale = 1, groups = Rubra$Elevation, labels = Rubra$Transect, ellipse = TRUE, circle = FALSE)+
  ggtitle("Principal Component Analysis (PCA)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "Principal Component Analysis (PCA)", subtitle = "Variables analysis along an elevational gradient")+
  labs(colour = "Elevation")+
  theme_classic()
#On voit une corr?lation n?gative entre div_sp et climat!
#Voir : https://www.soft-concept.com/surveymag/comment-lire-une-acp.html

#PCA avec les bonnes couleurs
Rubra3 = as.data.frame(cbind(Phytochemistry = Rubra$div_sp, WC = Rubra$Water.Content, SLA = Rubra$SLA, Leaf_C.N = Rubra$Leaf_C.N, Elevation = Rubra$Elevation, Altitude = Rubra$Altitude))
Rubra$cat3 <-cut(Rubra$Altitude, breaks=4, right=FALSE, labels=c(1:4))

#Renommer les axes de PCA
names(Rubra)[11]<-"Phytochemical richness" 
names(Rubra)[19]<-"Specific leaf area"
names(Rubra)[22]<-"Leaf water content"
names(Rubra)[26]<-"Leaf C/N ratio"

###adonis
fit = adonis(Rubra3[,1:4]~Rubra3$cat3, distance="euclidean", permutations=999, data=Rubra3) #je sais pas ce que c'est...
fit

#PCA avec transect, site, SLA, WC, Climat
Rubra.pca = prcomp(Rubra[,c(11, 19, 22, 26)], center = TRUE, scale. = TRUE)
ggbiplot(Rubra.pca, obs.scale = 1, var.scale = 1, groups = Rubra$cat3, labels = Rubra$Transect, ellipse = TRUE, circle = FALSE)+
  ggtitle("Principal Component Analysis (PCA)", subtitle = "Variables analysis along species niche optimum")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.subtitle = element_text(hjust = 0.5))+
  #labs(title = "Principal Component Analysis (PCA)", subtitle = "Variables analysis along species niche optimum")+
  labs(colour = "Elevation")+
  theme_classic()

#Extract axis 1 PCA:
library(ade4)
library(factoextra)
traits<-cbind(Rubra[,c(11, 19, 22, 26)])
dudi1 <- dudi.pca(traits, scale = TRUE, scan = FALSE, nf = 2)
p <- fviz_pca_biplot(dudi1)
print(p)
dudi1$li[,1]

Model <- lm(dudi1$li[,1]~Rubra$cat3)
summary(Model)
anova(Model)

traits2<-cbind(dudi1$li[,1], traits)
#write.table(traits2, "TraitsPC1.csv", sep = ";")

#X1 = dudi1$li[,1]
#X2 = dudi1$li[,2]
#X3 = Rubra$cat3
#manova(X1, X2, X3)

#Rubra.pca$X1[,1]
#PC1Traits <- data.frame(Rubra.pca$PC1, Rubra.pca$Phytochemical.richness)
#write.table(PC1Traits, "PC1Rubra.csv", sep = ";")
#PC1Traits <- data.frame(dudi1$li[,1], Rubra$Phytochemical.richness)
#PC1Traits <- data.frame(dudi1$li[,1])
#write.table(PC1Traits, "PC1Traits.csv", sep = ";")
#PC1Rubratraits <- data.frame(Rubra)
#write.table(PC1Rubratraits, "PC1Rubratraits.csv", sep = ";")

#Avec PC1 des traits

PC1Traits = read.csv("TraitsPC1Rubra.csv", sep = ";", h = T)

PC1Traits$cat3 <-cut(PC1Traits$Altitude , breaks=4, right=FALSE, labels=c(1:4))
ggplot(PC1Traits, aes(x=as.numeric(cat3), y = dudi1))+
  geom_boxplot(aes(fill = cat3, width=0.5))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^1))+
  labs(title = "Impact of elevation on plant functional traits",
       x="Altitude (m)",
       y="PC1 of leaf functional traits (Axis 1)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Elevation", 
                   limits=c("Low", "Mid", "High", "Very high"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

reg <- lmer(dudi1~Altitude+(1|Transect), data = PC1Traits)
summary(reg)
anova(reg)

#V?rification des donn?es pour faire les tests stat

qqnorm(PC1Traits$Altitude)
qqline(PC1Traits$Altitude)
hist(PC1Traits$Altitude)
bn <- bestNormalize(PC1Traits$Altitude)
bn
plot(bn, leg_loc = "bottomright")
PC1Traits$Altitude <- predict(bn)

qqnorm(PC1Traits$dudi1)
qqline(PC1Traits$dudi1)
hist(PC1Traits$dudi1)
bn <- bestNormalize(PC1Traits$dudi1)
bn
plot(bn, leg_loc = "bottomright")
PC1Traits$dudi1 <- predict(bn)

reg <- lmer(Phytochemical.richness~Altitude+(1|Transect), data = PC1Traits)
summary(reg)
anova(reg)

#(4) ################################################################################
#Corr?lations JB
######################################################################################

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tmap)
library(raster)
library(lme4)
library(lmerTest)
library(bestNormalize)
library(ape)
library(plot.matrix)
library(ade4)
library(heatmaply)
library(reshape2)
library(vegan)
library(rgl)
library(phytools)

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Interspecific")

JB = read.csv("FestucaJB_SLA.csv", sep = ";", h = T)
JB_CHN = read.csv("FestucaJB_CHN.csv", sep = ";", h = T)

#JB_CHN <- subset(JB_CHN, div_sp != 867)
#JB_CHN <- subset(JB_CHN, div_sp != 814)

qqnorm(JB$SLA)
qqline(JB$SLA)
hist(JB$SLA)
bn <- bestNormalize(JB$SLA)
bn
plot(bn, leg_loc = "bottomright")
JB$SLA <- predict(bn)

qqnorm(JB$DEM_Mean)
qqline(JB$DEM_Mean)
hist(JB$DEM_Mean)
bn <- bestNormalize(JB$DEM_Mean)
bn
plot(bn, leg_loc = "bottomright")
JB$DEM_Mean <- predict(bn)

qqnorm(JB$Water.Content)
qqline(JB$Water.Content)
hist(JB$Water.Content)
bn <- bestNormalize(JB$Water.Content)
bn
plot(bn, leg_loc = "bottomright")
JB$Water.Content <- predict(bn)

qqnorm(JB$div_sp)
qqline(JB$div_sp)
hist(JB$div_sp)
bn <- bestNormalize(JB$div_sp)
bn
plot(bn, leg_loc = "bottomright")
JB$div_sp <- predict(bn)

qqnorm(JB_CHN$div_sp)
qqline(JB_CHN$div_sp)
hist(JB_CHN$div_sp)
bn <- bestNormalize(JB_CHN$div_sp)
bn
plot(bn, leg_loc = "bottomright")
JB_CHN$div_sp <- predict(bn)

qqnorm(JB_CHN$Leaf_C.N)
qqline(JB_CHN$Leaf_C.N)
hist(JB_CHN$Leaf_C.N)

qqnorm(JB_CHN$DEM_Mean)
qqline(JB_CHN$DEM_Mean)
hist(JB_CHN$DEM_Mean)
bn <- bestNormalize(JB_CHN$DEM_Mean)
bn
plot(bn, leg_loc = "bottomright")
JB_CHN$DEM_Mean <- predict(bn)

qqnorm(JB_CHN$Climat)
qqline(JB_CHN$Climat)
hist(JB_CHN$Climat)
bn <- bestNormalize(JB_CHN$Climat)
bn
plot(bn, leg_loc = "bottomright")
JB_CHN$Climat <- predict(bn)

reg<-lm(Climat~DEM_Mean, data =JB_CHN)
summary(reg)
anova(reg) #Climat est le m?me pour chaque plante.species donc pas besoin de le mettre en facteur random

JB_CHN$cat3 <-cut(JB_CHN$DEM_Mean , breaks=4, right=FALSE, labels=c(1:4))
ggplot(JB_CHN, aes(y=Climat, x = cat3))+
  geom_boxplot(aes(fill=cat3))+
  scale_x_discrete(breaks=c("1","2","3", "4"),
                   labels=c("Low", "Mid", "High", "Very high"))+
  labs(title="Impact of elevation on climate at niche optimum", y = "Climate (PCA Axis 1)", x = "Altitudinal gradient")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_discrete(name = "Elevation", labels = c("Low", "Mid", "High", "Very high"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

reg<-lmer(Climat~DEM_Mean+(1|Plant.Species), data =JB_CHN)
summary(reg)
anova(reg)

#(1) SLA en fonction esp?ce

JB$cat3 <-cut(JB$DEM_Mean , breaks=4, right=FALSE, labels=c(1:4))
ggplot(JB, aes(x=as.numeric(cat3), y = SLA))+
  geom_boxplot(aes(fill = cat3))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^2))+
  labs(title="Phytochemistry richness according to elevation",
       x="Elevation",
       y="Phytochemistry richness (nb of molecules)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Elevation", 
                   limits=c("Low", "Mid", "High", "Very high"))

reg <- lmer(SLA~DEM_Mean+(1|Plant.Species), data=JB)
summary(reg)
anova(reg)

#(2) Water Content en fonction esp?ce

ggplot(JB, aes(x=Water.Content, y = Plant.Species, col = CatElev22))+
  geom_point(size=8)+
  labs(title="Leaf water content (WC) according to species",
       subtitle = "P-value = 2.2e-16",
       x="Leaf water content (%)",
       y="Plant Species")

reg <- lm(Water.Content~DEM_Mean, data=JB)
summary(reg)
anova(reg) #ATTENTION On ne met pas plante.species en random factor car on a seulement une valeur
#de water content par esp?ce et pas par individu.

#Leaf_C.N En fonction esp?ce

reg <- lmer(Leaf_C.N~DEM_Mean+(1|Plant.Species), data=JB_CHN)
summary(reg)
anova(reg)

#(3) div_sp en fonction esp?ce

JB_CHN$cat3 <-cut(JB_CHN$DEM_Mean , breaks=4, right=FALSE, labels=c(1:4))

ggplot(JB_CHN, aes(x=as.numeric(cat3), y = div_sp))+
  geom_boxplot(aes(fill = cat3))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^2))+
  labs(title="Phytochemistry richness according to elevation",
       x="Elevation",
       y="Phytochemistry richness (nb of molecules)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Elevation", 
                   limits=c("Low", "Mid", "High", "Very high"))

reg <- lm(div_sp~Plant.Species, data=JB_CHN)
summary(reg)
anova(reg)

library(ggplot2)
library(dplyr)
library(vegan)
library(remotes)
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(calibrate)
library(ggrepel)

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Interspecific")
JB = read.csv("FestucaJB_PCA.csv", sep = ";", h = T)

###adonis
fit = adonis(JB[,c(2, 18, 21)]~Plant.Species, distance="euclidean", permutations=999, data=JB, na.rm=TRUE) #je sais pas ce que c'est...
fit

#PCA avec SLA, WC, div_sp
JB.pca = prcomp(JB[,c(2, 18, 21)], center = TRUE, scale. = TRUE)
ggbiplot(JB.pca, obs.scale = 1, var.scale = 1, groups = JB$CatElev22, ellipse = TRUE, circle = FALSE)



fit = adonis(JB[,c(12, 15)]~Plant.Species, distance="euclidean", permutations=999, data=JB) #je sais pas ce que c'est...
#PCA avec transect, site, SLA, WC, Climat
JB.pca = prcomp(JB[,c(2, 18, 21)], center = TRUE, scale. = TRUE)
ggbiplot(JB.pca, obs.scale = 1, var.scale = 1, groups = JB$Plant.Species, ellipse = FALSE, circle = FALSE)


#####################################

setwd("D:/Affaires Nadline/Master thesis/Statistiques/Interspecific")
JB = read.csv("FestucaJB_PCA.csv", sep = ";", h = T)
JB2 = merge(JB, JB_CHN, by="div_sp")
JB2$cat3 <-cut(JB2$DEM_Mean.x, breaks=4, right=FALSE, labels=c(1:4))

#Renommer les axes de PCA
names(JB2)[1]<-"Phytochemical richness" 
names(JB2)[18]<-"Specific leaf area"
names(JB2)[21]<-"Leaf water content"
names(JB2)[44]<-"Leaf C/N ratio"

###adonis
fit = adonis(JB2[,c(1, 18, 21, 44)]~JB2$cat3, distance="euclidean", permutations=999, data=JB, na.rm=TRUE) #je sais pas ce que c'est...
fit

#PCA avec transect, site, SLA, WC, Climat
JB.pca = prcomp(JB2[,c(1, 18, 21, 44)], center = TRUE, scale. = TRUE)
ggbiplot(JB.pca, obs.scale = 1, var.scale = 1, groups = JB2$cat3, labels = JB2$Plant.Species.x, ellipse = TRUE, circle = FALSE)+
  ggtitle("Principal Component Analysis (PCA)", subtitle = "Variables analysis along species niche optimum")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.subtitle = element_text(hjust = 0.5))+
  #labs(title = "Principal Component Analysis (PCA)", subtitle = "Variables analysis along species niche optimum")+
  labs(colour = "Elevation")+
  theme_classic()

#Extract axis 1 PCA:
library(ade4)
library(factoextra)
traits<-cbind(JB2[,c(1, 18, 21, 44)])
dudi1 <- dudi.pca(traits, scale = TRUE, scan = FALSE, nf = 2)
p <- fviz_pca_biplot(dudi1)
print(p)
dudi1$li[,1]

Model <- lm(dudi1$li[,1]~JB2$cat3)
summary(Model)
anova(Model)

traits2<-cbind(dudi1$li[,1], traits)
write.table(traits2, "TraitsPC1.csv", sep = ";")

#X1 = dudi1$li[,1]
#X2 = dudi1$li[,2]
#X3 = Rubra$cat3
#manova(X1, X2, X3)

#Avec PC1 des traits

PC1Traits = read.csv("TraitsPC1Species.csv", sep = ";", h = T)

PC1Traits$cat3 <-cut(PC1Traits$DEM_Mean , breaks=4, right=FALSE, labels=c(1:4))
ggplot(PC1Traits, aes(x=as.numeric(cat3), y = dudi1))+
  geom_boxplot(aes(fill = cat3, width=0.5))+
  #geom_point(aes(col = (cat3), size=7))+
  geom_smooth(method = "lm", formula=y~x+I(x^1))+
  labs(title = "Impact of elevation on plant functional traits",
       x="Altitude (m)",
       y="PC1 of leaf functional traits (Axis 1)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name ="Elevation", 
                   limits=c("Low", "Mid", "High", "Very high"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

#V?rification des donn?es pour faire les tests stat

qqnorm(PC1Traits$DEM_Mean)
qqline(PC1Traits$DEM_Mean)
hist(PC1Traits$DEM_Mean)
bn <- bestNormalize(PC1Traits$DEM_Mean)
bn
plot(bn, leg_loc = "bottomright")
PC1Traits$DEM_Mean <- predict(bn)

qqnorm(PC1Traits$dudi1)
qqline(PC1Traits$dudi1)
hist(PC1Traits$dudi1)
bn <- bestNormalize(PC1Traits$dudi1)
bn
plot(bn, leg_loc = "bottomright")
PC1Traits$dudi1 <- predict(bn)

reg <- lmer(dudi1~DEM_Mean+(1|Plant.Species), data = PC1Traits)
summary(reg)
anova(reg)

#(5)
#######Graphes par cat?gories##################################################
###############################################################################

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
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

citation("factoextra")

#range01(s)

setwd("D:/Affaires Nadline/Master thesis/Statistiques/R?f?rences") 

Rubra = read.csv("Tableau Rubra.csv", sep = ";", h = T) ### merge data clim

climrubra = read.csv("climrubra.csv", sep = ";", h = T)

Species = read.csv("Tableau Festuca Complet.csv", sep = ";", h = T)### merge data clim

climoutputspecies = read.csv("climoutputspecies.csv", sep = ";", h = T)

ClimatRubra = merge(Rubra, climrubra, by = "Site")

ClimatSpecies = merge(Species, climoutputspecies, by = "Plant.Species")

#Leaf C/N with div_sp With boxplots

ClimatRubra$cat2 <-cut(ClimatRubra$Leaf_C.N, breaks=4, right=FALSE, labels=c(1:4))
ClimatSpecies$cat3 <-cut(ClimatSpecies$Leaf_C.N, breaks=4, right=FALSE, labels=c(1:4))


ggplot() +
  geom_boxplot(data=ClimatSpecies,
               aes(x = cat3, y = div_sp, group=cat3, col="cat3"), fill = "dodgerblue",alpha=0.2,width=0.5) +
  geom_boxplot(data=ClimatRubra,
               aes(x=cat2, y = div_sp, group = cat2, col="cat2"), fill = "gold",alpha=0.2,width=0.5) +
  
  geom_smooth(data=ClimatSpecies,
              aes(x=as.numeric(cat3), y=div_sp), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="dodgerblue", size=1.5) +
  geom_smooth(data=ClimatRubra,
              aes(x=as.numeric(cat2), y=div_sp), method="lm",
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

reg <- lm(div_sp~as.numeric(cat3), data = ClimatSpecies)
summary(reg)
anova(reg) #Significatif

#Leaf C/N with altitude

ClimatRubra$cat2 <-cut(ClimatRubra$Altitude.x , breaks=4, right=FALSE, labels=c(1:4))
ClimatSpecies$cat3 <-cut(ClimatSpecies$DEM_Mean.x, breaks=4, right=FALSE, labels=c(1:4))


ggplot() +
  geom_boxplot(data=ClimatSpecies,
               aes(x = cat3, y = Leaf_C.N, group=cat3, col = "cat3"), fill = "dodgerblue",alpha=0.2,width=0.5) +
  geom_boxplot(data=ClimatRubra,
               aes(x=cat2, y = Leaf_C.N, group = cat2, col = "cat2"), fill = "gold",alpha=0.2,width=0.5) +
  
  geom_smooth(data=ClimatSpecies,
              aes(x=as.numeric(cat3), y=Leaf_C.N), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="dodgerblue", size=1.5) +
  geom_smooth(data=ClimatRubra,
              aes(x=as.numeric(cat2), y=Leaf_C.N), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="gold", size=1.5,alpha = 0.1)+
  theme_classic()+
  labs(title="Leaf C/N according to elevation",
       x="Elevation",
       y="Leaf C/N (ratio)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name="Elevation",
                   breaks=c("1", "2", "3", "4"),
                   labels=c("Low", "Mid", "High", "Very high"))+
  labs(colour = "Plant type")+
  scale_color_manual(labels = c("F. rubra", "Festuca sp."), values = c("gold", "dodgerblue"))

reg <- lmer(Leaf_C.N~as.numeric(cat2)+(1|Transect.x), data = ClimatRubra)
summary(reg)
anova(reg) #Significatif

reg <- lm(Leaf_C.N~as.numeric(cat3), data = ClimatSpecies)
summary(reg)
anova(reg) #Non significatif

anova(lm(ClimatSpecies$Leaf_C.N ~ poly(as.numeric(ClimatSpecies$cat3), 2, raw=TRUE)))#Non significatif

#anova(lm(ClimatSpecies$Leaf_C.N ~ poly((as.numeric(ClimatSpecies$cat3)), 2, raw=TRUE)))
#anova(lm(ClimatRubra$Leaf_C.N ~ poly((as.numeric(ClimatRubra$cat2)), 2, raw=TRUE)))

#Pour SLA
SpeciesSLAWC = read.csv("Tableau Festuca SLA and WC.csv", sep = ";", h = T)

ggplot(ClimatRubra, aes(x = SLA, y = div_sp)) +
  geom_point(col = "red") +
  geom_smooth(data=ClimatRubra,
              aes(x=SLA, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5) +
  geom_point(data=SpeciesSLAWC,
             aes(x=SLA, y=div_sp), color="turquoise", size=1.5) +
  geom_smooth(data=SpeciesSLAWC,
              aes(x=SLA, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5)

ClimatRubra$cat2 <-cut(ClimatRubra$SLA, breaks=4, right=FALSE, labels=c(1:4))
SpeciesSLAWC$cat3 <-cut(SpeciesSLAWC$SLA, breaks=4, right=FALSE, labels=c(1:4))

ggplot() +
  geom_boxplot(data=SpeciesSLAWC,
               aes(x = cat3, y = div_sp, group=cat3, col="cat3"), fill = "dodgerblue",alpha=0.2,width=0.5) +
  geom_boxplot(data=ClimatRubra,
               aes(x=cat2, y = div_sp, group = cat2, col="cat2"), fill = "gold",alpha=0.2,width=0.5) +
  
  geom_smooth(data=SpeciesSLAWC,
              aes(x=as.numeric(cat3), y=div_sp), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="dodgerblue", size=1.5) +
  geom_smooth(data=ClimatRubra,
              aes(x=as.numeric(cat2), y=div_sp), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="gold", size=1.5,alpha = 0.1)+
  coord_cartesian(ylim=c(200,850)) +
  theme_classic()+
  labs(title="Specific leaf area according to phytochemistry diversity",
       x="Specific leaf area (unit)",
       y="Phytochemistry diversity (nb of molecules)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name="SLA (m^2/kg)",
                   breaks=c("1", "2", "3", "4"),
                   labels=c("Small", "Mid", "Large", "Very large"))+
  labs(colour = "Plant type")+
  scale_color_manual(labels = c("F. rubra", "Festuca sp."), values = c("gold", "dodgerblue"))

reg <- lmer(div_sp~SLA+(1|Transect.x), data=ClimatRubra)
summary(reg)
anova(reg)#Pas significatif

reg <- lm(div_sp~SLA, data=SpeciesSLAWC)
summary(reg)
anova(reg)#Pas significatif

#SLA avec altitude (si on remplace SLA par Water.Content on a pour water content).

ClimatRubra$cat2 <-cut(ClimatRubra$Altitude.x , breaks=4, right=FALSE, labels=c(1:4))
SpeciesSLAWC$cat3 <-cut(SpeciesSLAWC$DEM_Mean.x, breaks=4, right=FALSE, labels=c(1:4))


ggplot() +
  geom_boxplot(data=ClimatSpecies,
               aes(x = cat3, y = SLA, group=cat3, col = "cat3"), fill = "dodgerblue",alpha=0.2,width=0.5) +
  geom_boxplot(data=ClimatRubra,
               aes(x=cat2, y = SLA, group = cat2, col = "cat2"), fill = "gold",alpha=0.2,width=0.5) +
  
  geom_smooth(data=ClimatSpecies,
              aes(x=as.numeric(cat3), y=Leaf_C.N), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="dodgerblue", size=1.5) +
  geom_smooth(data=ClimatRubra,
              aes(x=as.numeric(cat2), y=Leaf_C.N), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="gold", size=1.5,alpha = 0.1)+
  theme_classic()+
  labs(title="Leaf C/N according to elevation",
       x="Elevation",
       y="Leaf C/N (ratio)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name="Elevation",
                   breaks=c("1", "2", "3", "4"),
                   labels=c("Low", "Mid", "High", "Very high"))+
  labs(colour = "Plant type")+
  scale_color_manual(labels = c("F. rubra", "Festuca sp."), values = c("gold", "dodgerblue"))

reg <- lmer(Leaf_C.N~as.numeric(cat2)+(1|Transect.x), data = ClimatRubra)
summary(reg)
anova(reg) #Significatif

reg <- lm(Leaf_C.N~as.numeric(cat3), data = ClimatSpecies)
summary(reg)
anova(reg)

#Pour WC
SpeciesSLAWC = read.csv("Tableau Festuca SLA and WC.csv", sep = ";", h = T)

ggplot(ClimatRubra, aes(x = Water.Content, y = div_sp)) +
  geom_point(col = "red") +
  geom_smooth(data=ClimatRubra,
              aes(x=Water.Content, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5) +
  geom_point(data=SpeciesSLAWC,
             aes(x=Water.Content, y=div_sp), color="turquoise", size=1.5) +
  geom_smooth(data=SpeciesSLAWC,
              aes(x=Water.Content, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5)

ClimatRubra$cat2 <-cut(ClimatRubra$Water.Content, breaks=4, right=FALSE, labels=c(1:4))
SpeciesSLAWC$cat3 <-cut(SpeciesSLAWC$Water.Content, breaks=4, right=FALSE, labels=c(1:4))

ggplot() +
  geom_boxplot(data=SpeciesSLAWC,
               aes(x = cat3, y = div_sp, group=cat3, col = "cat3"), fill = "dodgerblue",alpha=0.2,width=0.5) +
  geom_boxplot(data=ClimatRubra,
               aes(x=cat2, y = div_sp, group = cat2, col = "cat2"), fill = "gold",alpha=0.2,width=0.5) +
  
  geom_smooth(data=SpeciesSLAWC,
              aes(x=as.numeric(cat3), y=div_sp), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="dodgerblue", size=1.5, alpha = 0.1) +
  geom_smooth(data=ClimatRubra,
              aes(x=as.numeric(cat2), y=div_sp), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="gold", size=1.5,alpha = 0.1)+
  coord_cartesian(ylim=c(200,850)) +
  theme_classic()+
  labs(title="Leaf water content according to phytochemistry diversity",
       x="Specific leaf area (unit)",
       y="Phytochemistry diversity (nb of molecules)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name="Leaf water content (%)",
                   breaks=c("1", "2", "3", "4"),
                   labels=c("Small", "Mid", "Large", "Very large"))+
  labs(colour = "Plant type")+
  scale_color_manual(labels = c("F. rubra", "Festuca sp."), values = c("gold", "dodgerblue"))

reg <- lmer(div_sp~Water.Content+(1|Transect.x), data=ClimatRubra)
summary(reg)
anova(reg)#Pas significatif

reg <- lm(div_sp~Water.Content, data=SpeciesSLAWC)
summary(reg)
anova(reg)#Pas significatif

#Pour Climat

ggplot(ClimatRubra, aes(x = Climat, y = div_sp)) +
  geom_point(col = "red") +
  geom_smooth(data=ClimatRubra,
              aes(x=Climat, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5) +
  geom_point(data=ClimatSpecies,
             aes(x=Climat, y=div_sp), color="turquoise", size=1.5) +
  geom_smooth(data=SpeciesSLAWC,
              aes(x=Climat, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5)

ClimatRubra$cat2 <-cut(ClimatRubra$Climat, breaks=4, right=FALSE, labels=c(1:4))
SpeciesSLAWC$cat3 <-cut(SpeciesSLAWC$Climat, breaks=4, right=FALSE, labels=c(1:4))

ggplot() +
  geom_boxplot(data=SpeciesSLAWC,
               aes(x = cat3, y = div_sp, group=cat3, col = "cat3"), fill = "dodgerblue",alpha=0.2,width=0.5) +
  geom_boxplot(data=ClimatRubra,
               aes(x=cat2, y = div_sp, group = cat2, col="cat2"), fill = "gold",alpha=0.2,width=0.5) +
  
  geom_smooth(data=SpeciesSLAWC,
              aes(x=as.numeric(cat3), y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="dodgerblue", size=1.5) +
  geom_smooth(data=ClimatRubra,
              aes(x=as.numeric(cat2), y=div_sp), method="lm",
              formula = y ~ poly(x,1), n= 40, se=TRUE, color="gold", size=1.5,alpha = 0.1)+
  coord_cartesian(ylim=c(200,850)) +
  theme_classic()+
  labs(title="Impact of climate on phytochemical richness",
       x="Climate (PC1)",
       y="Phytochemical richness (nb of molecules)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(name="Climate (PC1)",
                   breaks=c("1", "2", "3", "4"),
                   labels=c("Warm+", "Warm", "Cold", "Cold+"))+
  labs(colour = "Plant type")+
  scale_color_manual(labels = c("F. rubra", "Festuca sp."), values = c("gold", "dodgerblue"))+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

reg <- lmer(div_sp~Climat+(1|Transect.x), data=ClimatRubra)
summary(reg)
anova(reg)#Pas significatif

reg <- lm(div_sp~Climat, data=ClimatRubra)
summary(reg)
anova(reg)

reg <- lmer(div_sp~as.numeric(cat2)+(1|Transect.x), data=ClimatRubra)
summary(reg)
anova(reg)#Pas significatif

anova(lm(ClimatSpecies$div_sp ~ poly(as.numeric(ClimatSpecies$cat3), 2, raw=TRUE))) #Significatif

#Sans cat?gories (points)

ggplot() +
  geom_point(data=Species,
             aes(x = DEM_Mean, y = div_sp), col = "turquoise") +
  geom_point(data=Rubra,
             aes(x=Altitude, y = div_sp), col = "red")+
  geom_smooth(data=Species,
              aes(x=DEM_Mean, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5) +
  geom_smooth(data=Rubra,
              aes(x=Altitude, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5)+
  labs(title="Effect of elevation on phytochemistry richness",
       x="Elevation",
       y="Phytochemistry richness (nb of molecules)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

#Sans cat?gories mais moyenne des points (points)
Rubra$cat2 <-cut(Rubra$Altitude , breaks=4, right=FALSE, labels=c(1:4))
Species$cat3 <-cut(Species$DEM_Mean , breaks=4, right=FALSE, labels=c(1:4))

ggplot() +
  geom_point(data=Species,
             aes(x = DEM_Mean, y = Mean), col = "turquoise") +
  geom_point(data=Rubra,
             aes(x=Altitude, y = Mean), col = "red")+
  geom_smooth(data=Species,
              aes(x=DEM_Mean, y=Mean), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5) +
  geom_smooth(data=Rubra,
              aes(x=Altitude, y=Mean), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5)+
  labs(title="Effect of elevation on mean phytochemistry richness",
       x="Elevation",
       y="Mean phytochemistry richness (nb of molecules)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

#Pour avoir le nom des esp?ces
ggplot(Species, aes(x= DEM_Mean, y= div_sp, colour="green", label=Plant.Species)) +
  geom_point(size = 2,alpha = 0.6) +
  theme_bw()+
  geom_text(aes(label=Plant.Species),hjust=0, vjust=0)


#Graphe au propre

Rubra$cat2 <-cut(Rubra$Altitude , breaks=4, right=FALSE, labels=c(1:4))
Species$cat3 <-cut(Species$DEM_Mean , breaks=4, right=FALSE, labels=c(1:4))
new <- data.frame(Species[,c(3, 6, 23)])

ggplot() +
  geom_boxplot(data=Species,
               aes(x = cat3, y = div_sp, group=cat3, col = "cat3"), fill = "dodgerblue",alpha=0.2,width=0.5) +
  geom_boxplot(data=Rubra,
               aes(x=cat2, y = div_sp, group = cat2, col = "cat2"), fill = "gold",alpha=0.2,width=0.5) +
  
  geom_smooth(data=Species,
              aes(x=as.numeric(cat3), y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=FALSE, color="dodgerblue", size=1.5) +
  geom_smooth(data=Rubra,
              aes(x=as.numeric(cat2), y=div_sp), method="lm",
              formula = y ~ poly(x,1), n= 40, se=FALSE, color="gold", size=1.5,alpha = 0.1)+
  coord_cartesian(ylim=c(200,850)) +
  theme_classic()+
  labs(title="Impact of elevation on phytochemical richness",
       x="Elevation",
       y="Phytochemical richness (nb of molecules)")+
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

#Graphe boxplot sur la moyenne
Rubra$cat2 <-cut(Rubra$Altitude , breaks=4, right=FALSE, labels=c(1:4))
Species$cat3 <-cut(Species$DEM_Mean , breaks=4, right=FALSE, labels=c(1:4))


ggplot() +
  geom_boxplot(data=Species,
               aes(x = cat3, y = Mean, group=cat3, col = "cat3"), fill = "dodgerblue",alpha=0.2,width=0.5) +
  geom_boxplot(data=Rubra,
               aes(x=cat2, y = Mean, group = cat2, col = "cat2"), fill = "gold",alpha=0.2,width=0.5) +
  
  geom_smooth(data=Species,
              aes(x=as.numeric(cat3), y=Mean), method="lm",
              formula = y ~ poly(x,2), n= 40, se=FALSE, color="dodgerblue", size=1.5) +
  geom_smooth(data=Rubra,
              aes(x=as.numeric(cat2), y=Mean), method="lm",
              formula = y ~ poly(x,1), n= 40, se=FALSE, color="gold", size=1.5,alpha = 0.1)+
  coord_cartesian(ylim=c(200,850)) +
  theme_classic()+
  labs(title="Impact of elevation on mean phytochemical richness",
       x="Altitudinal gradient",
       y="Mean phytochemical richness (nb of molecules)")+
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
#scale_color_hue(labels = c("F. rubra", "Festuca sp."))

anova(lm(Rubra$Mean~Rubra$Altitude))

reg <- lmer(Mean~as.numeric(cat2)+(1|Transect), data = Rubra)
summary(reg)
anova(reg)

anova(lm(Species$Mean ~ poly(as.numeric(Species$cat3), 2, raw=TRUE)))


#-------------------------------------------------------------------------------------------------------

#Axe 1 PCA (PC1) structuration chimique

PC1Rubra = read.csv("PC1Rubra.csv", sep = ";", h = T) ### merge data clim

PC1Species = read.csv("PC1Species.csv", sep = ";", h = T)

PC1Species$X1 <- range01(PC1Species$X1)
PC1Rubra$X1 <- range01(PC1Rubra$X1)

ggplot(PC1Rubra, aes(x = X1, y = div_sp)) +
  geom_point(col = "red") +
  geom_smooth(data=PC1Rubra,
              aes(x=X1, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5) +
  geom_point(data=PC1Species,
             aes(x=X1, y=div_sp), color="turquoise", size=1.5) +
  geom_smooth(data=PC1Species,
              aes(x=X1, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5)

#Seulement pour les esp?ces

PC1Species = read.csv("PC1Species.csv", sep = ";", h = T)

PC1Species$cat3 <-cut(PC1Species$X1 , breaks=4, right=FALSE, labels=c(1:4))

ggplot() +
  geom_point(data=PC1Species,
             aes(x=X1, y=div_sp), color="turquoise", size=1.5) +
  geom_smooth(data=PC1Species,
              aes(x=X1, y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5)

#En boxplots

PC1Species$cat3 <-cut(PC1Species$X1 , breaks=4, right=FALSE, labels=c(1:4))

ggplot() +
  geom_boxplot(data=PC1Species,
               aes(x=cat3, y=div_sp, group=cat3), color="turquoise", size=1.5) +
  geom_smooth(data=PC1Species,
              aes(x=as.numeric(cat3), y=div_sp), method="lm",
              formula = y ~ poly(x,2), n= 40, se=TRUE, color="turquoise", size=1.5)+
  theme_classic()+
  labs(title="Impact of chemical structure on phytochemical richness",
       x="PC1 chemical structure",
       y="Phytochemical richness (nb of molecules)")+
  theme(plot.title = element_text(hjust = 0.5, size = 25))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))

reg <- lm(div_sp~X1, data=PC1Species)
summary(reg)
anova(reg)

qqnorm(PC1Species$X1)
qqline(PC1Species$X1)
hist(PC1Species$X1)
bn <- bestNormalize(PC1Species$X1)
bn
plot(bn, leg_loc = "bottomright")
PC1Species$X1 <- predict(bn)

qqnorm(PC1Species$div_sp)
qqline(PC1Species$div_sp)
hist(PC1Species$div_sp)
bn <- bestNormalize(PC1Species$div_sp)
bn
plot(bn, leg_loc = "bottomright")
PC1Species$div_sp <- predict(bn)

