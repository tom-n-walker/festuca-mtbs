################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Load phylogeny
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    22 February 2022
#### ---------------------------------------------------------------------------

build_phylogeny <- function(rawPlants){
  ## Build species list ----
  species <- rawPlants$cross %>%
    select(species) %>%
    mutate(genus = "Festuca") %>%
    mutate(family = "Poaceae") %>%
    distinct(species, .keep_all = T)
  ## Build tree ----
  tree <- phylo.maker(
    sp.list = species, 
    scenarios = "S3"
  )
  ## Format tree ----
  formatTree <- tree$scenario.3 %>%
    # solve dichotomies/trichotomies
    multi2di(., random = FALSE) %>%
    # make unique node labels
    makeNodeLabel
  # add 0.001 to all edge lengths to eliminate zeros
  formatTree$edge.length <- formatTree$edge.length + 0.001
  # final formatting of labels
  formatTree$tip.label <- gsub("_", " ", formatTree$tip.label)
  ## Return ----
  return(formatTree)
}
