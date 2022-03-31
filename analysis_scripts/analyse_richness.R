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
library(vegan)
library(nlme)
library(data.table)
library(tidyverse)
# load mm2in converter
source("processing_scripts/mm2in.R")

#### DATA ----------------------------------------------------------------------

drake::loadd()


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
# build base models with climate plus structure
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
# add elevation to see if it additionally controls richness
pgls2 <- update(pgls1, ~.+ elevDisc)
lme2 <- update(lme1, ~.+ elevCont)
# test significance
anova(pgls1, pgls2)
anova(lme1, lme2)


#### PLOT ----------------------------------------------------------------------

## Build summary dataset for cross-species ----
crossPlotData <- crossAgg$data %>%
  mutate(elevLet = toupper(substr(elevDisc, 1, 1))) %>%
  mutate(elevLet = ifelse(elevLet == "V", "VH", elevLet)) %>%
  mutate(elevLet = factor(elevLet, levels = c("L", "M", "H", "VH"))) %>%
  group_by(elevLet) %>%
  summarise(
    mean = mean(mtbsRich, na.rm = T),
    se = sd(mtbsRich, na.rm = T)/sqrt(n())
  ) %>%
  ungroup

## Build plots ----
crossPlot <- ggplot(crossPlotData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(300, 700)) +
  aes(x = elevLet, y = mean, ymax = mean + se, ymin = mean - se) +
  geom_bar(stat = "identity", col = "black", fill = "#F67280") +
  geom_errorbar(width = 0) +
  labs(x = "", y = "Phytochemical richness")
rubraPlot <- ggplot(rubraData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = elevCont, y = mtbsRich) +
  geom_point(shape = 21, fill = "#F67280") +
  geom_smooth(method = "lm", se = F, col = "black", size = 0.5) +
  labs(x = "Elevation (m a.s.l)", y = "Phytochemical richness")

## Write to file ----
postscript(
  file = "./figure_builds/richness.eps",
  width = mm2in(120),
  height = mm2in(60)
)
cowplot::plot_grid(crossPlot, rubraPlot)
dev.off()
