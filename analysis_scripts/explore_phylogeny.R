################################################################################
#### Project: Festuca metabolites
#### Title:   Drake plan
#### Author:  Tom Walker (thomas.walker@unine.ch)
#### Date:    08 March 2022
#### ---------------------------------------------------------------------------


#### PROLOGUE ------------------------------------------------------------------

## Options ----
# remove objects from global environment
rm(list = ls())
# R session options (no factors, bias against scientific #s)
options(
  stringsAsFactors = F,
  scipen = 6
)

## Libraries ----
# standard library set
library(caper)
library(data.table)
library(tidyverse)

#### DATA ----------------------------------------------------------------------

drake::loadd()


#### CLIMATE PCA ---------------------------------------------------------------

## Do PCA ----
crossPCA <- prcomp(cross$clim, scale = T, center = T)
rubraPCA <- prcomp(rubra$clim, scale = T, center = T)

## Plot PCA ----
# variance explained
par(mfrow = c(1, 2))
plot(crossPCA)
plot(rubraPCA)
# biplots
par(mfrow = c(1, 2))
biplot(crossPCA)
biplot(rubraPCA)

## Extract scores ----
crossClimScores <- crossPCA$x %>%
  as.data.frame %>%
  select(PC1:PC4)
rubraClimScores <- rubraPCA$x %>%
  as.data.frame %>%
  select(PC1:PC4)

## Assemble scores ----
crossData <- bind_cols(
  cross$meta, 
  cross$clim,
  crossClimScores,
  cross$mvData
)
rubraData <- bind_cols(
  rubra$meta, 
  rubra$clim,
  rubraClimScores,
  rubra$mvData
)


#### MANTEL TEST ---------------------------------------------------------------

## Build datasets ----
# metabolite full matrix
mtbsMatrix <- mtbsCross$relAbun %>%
  mutate(sampleID = mtbsCross$sampleID) %>%
  left_join(cross$meta, ., by = c("mtbsSampleID" = "sampleID")) %>%
  select(-species, -plantID, -mtbsSampleID, -elevDisc)
# metabolite composition
mtbs <- mtbsCross$relAbun %>%
  mutate(sampleID = mtbsCross$sampleID) %>%
  left_join(cross$meta, ., by = c("mtbsSampleID" = "sampleID")) %>%
  select(-plantID, -mtbsSampleID, -elevDisc) %>%
  group_by(species) %>%
  summarise(across(everything(), mean, na.rm = T)) %>%
  ungroup %>%
  arrange(species) %>%
  select(-species) %>%
  as.matrix %>%
  parallelDist::parDist(mtbs, "bray")
# climate
clim <- cross$meta %>%
  select(species) %>%
  bind_cols(., cross$clim) %>%
  group_by(species) %>%
  summarise(across(everything(), mean)) %>%
  ungroup %>%
  arrange(species) %>%
  select(-species) %>%
  dist
# phylogeny
phyl <- phylogeny$phyloDist

## Build model ----
model <- phytools::multi.mantel(Y = mtbs, X = list(phyl, clim))
model

## PERMANOVA ----
# drop NAs
crossMtbsNA <- drop_na(mtbsMatrix)
crossDataNA <- crossData[-which(is.na(mtbsMatrix[, 1])), ]

pm <- adonis(
  crossMtbsNA ~ elevDisc + species,
  crossDataNA,
  na.rm = T
)
pm


nmds <- metaMDS(crossMtbsNA)
plot(nmds, display = "sites")
with(crossMtbsNA, ordiellipse(nmds, crossDataNA$species, label = T))
plot(nmds, display = "sites")
with(crossMtbsNA, ordiellipse(nmds, crossDataNA$elevDisc, label = T))

nmds2 <- metaMDS(trymax = 1000, mtbsRubra$relAbun)
plot(nmds, display = "sites")
with(crossMtbsNA, ordiellipse(nmds, crossDataNA$species, label = T))
plot(nmds, display = "sites")
with(mtbsRubra$relAbun, ordiellipse(nmds2, cut(rubraData$endoRich, 5), label = T))
orditorp(nmds2, display="sites", col="red", air=0.01)



plot(phylogeny$phyloTree)



m1 <- lme(
  mtbsPCoA1 ~ elevDisc,
  random = ~ 1 | species,
  crossData,
  method = "ML",
  na.action = "na.exclude"
)
drop1(m1, test = "Chisq")

m1 <- gls(mtbsPCoA1 ~ elevDisc + species, crossData, na.action = "na.exclude")
anova(m1)


#### MODEL ---------------------------------------------------------------------

## Aggregate data at the species level ----
crossSpecies <- crossData %>%
  select(-mtbsSampleID, -plantID) %>%
  select(-contains("phylo")) %>%
  group_by(species, elevDisc) %>%
  summarise(across(everything(), mean, na.rm = T)) %>%
  ungroup %>%
  as.data.frame
crossAgg <- comparative.data(
  phy = phylogeny$phyloTree,
  data = crossSpecies,
  names.col = species,
  vcv = T,
  vcv.dim = 3
)

## PGLS ----
# build base models
pgls1 <- pgls(
  mtbsRich ~ srad + tave + prec + soilF + soilR,
  crossAgg
)
lme1 <- lme(
  mtbsRich ~ srad + tave + prec + soilF + soilR, 
  random = ~ 1 | transect,
  rubraData,
  method = "ML"
)
# test significance
anova(pgls1, update(pgls1, ~.+ elevDisc))
anova(lme1, update(lme1, ~.+ elevCont))

anova(pgls1)
drop1(lme1, test = "Chisq")


par(mfrow = c(1, 2))
plot(mtbsRich ~ srad, crossData)
abline(lm(mtbsRich ~ srad, crossData))

phylosig(
  phylogeny$phyloTree,
  x = crossAgg$data$mtbsPCoA2,
  method = "K",
  test = T
)





