################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Load phylogeny
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    22 February 2022
#### ---------------------------------------------------------------------------

build_phylogeny <- function(rawPlants){
  ## Load phylogeny ----
  # load
  rawPhylo <- read.tree("./data/ITS_FestucaSwissIQTREE.contree")
  # replace underscores in species names
  rawPhylo$tip.label <- gsub("_", " ", rawPhylo$tip.label)
  ## Organise species ----
  # extract species lists
  phyloSpecies <- rawPhylo$tip.label
  studySpecies <- unique(rawPlants$cross$species)
  # extract species not present in study
  notPresent <- setdiff(phyloSpecies, studySpecies)
  ## Format tree ----
  # drop species absent from study
  cutPhylo <- drop.tip(rawPhylo, notPresent)
  # force ultrametric
  ultraTree <- chronopl(
    cutPhylo, 
    lambda = 0.5, 
    age.max = 10,
    node = "root", 
    tol = 1e-8,
    CV = FALSE, 
    eval.max = 500, 
    iter.max = 500
  )
  ## Return ----
  return(ultraTree)
}
