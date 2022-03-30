








crossTest <- bind_cols(cross$metaData, cross$bigData) %>%
  select(species, mtbsPCoA1, mtbsPCoA2, mtbsRich) %>%
  group_by(species) %>%
  summarise(across(everything(), mean, na.rm = T))

mtbsRichVector <- crossTest$mtbsPCoA1
names(mtbsRichVector) <- crossTest$species


fanTree <- phylogeny$phyloTree %>%
  contMap(., mtbsRichVector, plot = F) %>%
  setMap(., invert = T)
plot(
  fanTree,
  outline = F
  )
phylosig(
  phylogeny$phyloTree,
  mtbsRichVector,
  method = "lambda",
  test = T
)



sp_div_mean <- tapply(Metabo_stat$div_sp ,Metabo_stat$Plant.Species, mean) ## calculer la moyenne de phytodiversit? par esp?ces
names(sp_div_mean) <- gsub(" ", "_",names(sp_div_mean))
names(sp_div_mean) <- gsub("Festuca_valesciaca", "Festuca_valesiaca",names(sp_div_mean))

obj<-contMap(tr,sp_div_mean,plot=FALSE)
obj<-setMap(obj,invert=TRUE)

plot(obj,type="fan",outline=FALSE)
phylosig(tr,sp_div_mean,method="lambda",test=TRUE) #YEY!


## To do ----
- multi mantel test between phylo, climate on mtbs
- PGLS of climate on mtbs richness (takes into account phylogeny)
- shannon entropy calculations
- can we extract festuca from community datasets (from PNAS) in co-existing festuca?
  
  
  




rubra <- readd(rubra)
cross <- readd(cross)

anova(lm(mtbsRich ~ phyloPCoA1 + phyloPCoA2 + phyloPCoA3, cross$bigData))
anova(lm(mtbsPCoA1 ~ phyloPCoA1 + phyloPCoA2 + phyloPCoA3, cross$bigData))
anova(lm(mtbsPCoA2 ~ phyloPCoA1 + phyloPCoA2 + phyloPCoA3, cross$bigData))

par(mfrow = c(3, 3))
plot(mtbsRich ~ phyloPCoA1, cross$bigData)
plot(mtbsRich ~ phyloPCoA2, cross$bigData)
plot(mtbsRich ~ phyloPCoA3, cross$bigData)
plot(mtbsPCoA1 ~ phyloPCoA1, cross$bigData)
plot(mtbsPCoA1 ~ phyloPCoA2, cross$bigData)
plot(mtbsPCoA1 ~ phyloPCoA3, cross$bigData)
abline(lm(mtbsPCoA1 ~ phyloPCoA3, cross$bigData))
plot(mtbsPCoA2 ~ phyloPCoA1, cross$bigData)
plot(mtbsPCoA2 ~ phyloPCoA2, cross$bigData)
abline(lm(mtbsPCoA2 ~ phyloPCoA2, cross$bigData))
plot(mtbsPCoA2 ~ phyloPCoA3, cross$bigData)

par(mfrow = c(2, 3))

anova(lm(mtbsRich ~ phyloPCoA2, cross$bigData))
anova(lm(mtbsRich ~ phyloPCoA3, cross$bigData))
plot(mtbsRich ~ phyloPCoA1, cross$bigData)
plot(mtbsRich ~ phyloPCoA2, cross$bigData)
plot(mtbsRich ~ phyloPCoA3, cross$bigData)

crossClim <- cross$climate %>%
  select(contains("Mean", ignore.case = F)) %>%
  select(-contains("Forest")) %>%
  select(-TPI_Mean, -TRI_Mean, -TWI_Mean)
crossClimPCA <- prcomp(crossClim, scale = T, center = T)



testData <- cbind(rubra$bigData, rubra$metaData)
m1 <- lm(mtbsRich ~ elevCont + transect, testData)
anova(m1)


climPCAcross <- prcomp(cross$climate, scale = T, center = T)
summary(climPCAcross)
climPCArubra <- prcomp(rubra$climate, scale = T, center = T)
summary(climPCArubra)
plot(climPCArubra)

rubraTest <- bind_cols(rubra$metaData, rubra$bigData) %>%
  mutate(climPC1 = climPCArubra$x[, 1]) %>%
  mutate(climPC2 = climPCArubra$x[, 2])
crossTest <- bind_cols(cross$metaData, cross$bigData) %>%
  mutate(climPC1 = climPCAcross$x[, 1]) %>%
  mutate(climPC2 = climPCAcross$x[, 2]) %>%
  mutate(DEM = cross$climate$DEM_Mean) %>%
  mutate(elevDisc = factor(elevDisc, levels = c("Low", "Mid", "High", "Very high")))


m1 <- lm(mtbsRich ~ climPC1 + climPC2 + species, crossTest)
anova(m1)
m1a <- lme(
  mtbsRich ~ climPC1 + climPC2,
  random = ~ 1 | species,
  method = "ML",
  data = crossTest,
  na.action = "na.exclude"
)
par(mfrow = c(1, 2))
plot(residuals(m1a, type = "normalized") ~ fitted(m1a))
hist(residuals(m1a, type = "normalized"))

anova(m1a, update(m1a, ~.- climPC2))



m2 <- lm(mtbsRich ~ elevCont + transect, rubraTest)
anova(m1)
anova(m2)
emmeans(m1, pairwise ~ elevDisc)

crossSumm <- crossTest %>%
  group_by(elevDisc) %>%
  summarise(
    mean = mean(mtbsRich, na.rm = T),
    se = sd(mtbsRich, na.rm = T)/sqrt(n())
    ) %>%
  ungroup


myPal <- wes_palette("Darjeeling1")
crossPlot <- ggplot(crossSumm) +
  aes(x = elevDisc, y = mean, ymax = mean + se, ymin = mean - se) +
  theme_bw() +
  coord_cartesian(ylim = c(300, 700)) +
  theme(panel.grid = element_blank(), legend.title = element_blank()) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "grey") +
  labs(x = "Elevation", y = "Metabolite richness")
rubraPlot <- ggplot(rubraTest) +
  aes(x = elevCont, y = mtbsRich, fill = transect) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  scale_fill_manual(values = myPal) +
  geom_point(shape = 21) +
  geom_abline(slope = -0.03868759, intercept = 607.45350148, col = myPal[1]) +
  geom_abline(slope = -0.03868759, intercept = 607.45350148 + 9.37143085, col = myPal[2]) +
  geom_abline(slope = -0.03868759, intercept = 607.45350148 - 65.7628491, col = myPal[3]) +
  geom_abline(slope = -0.03868759, intercept = 607.45350148 - 28.4616116, col = myPal[4]) +
  geom_abline(slope = -0.03868759, intercept = 607.45350148 - 54.9415939, col = myPal[5]) +
  labs(x = "Elevation (m ASL)", y  = "Metabolite richness")

cowplot::plot_grid(crossPlot, rubraPlot)

