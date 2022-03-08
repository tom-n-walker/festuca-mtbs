################################################################################
#### Project: Festuca metabolites
#### Title:   Drake plan
#### Author:  Tom Walker (thomas.walker@unine.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------


#### PROLOGUE ------------------------------------------------------------------

## Options ----
# remove objects from global environment
rm(list = ls())

# configure default R session options (no factors, bias against scientific #s)
options(
  stringsAsFactors = F,
  scipen = 6
)

## Libraries ----
source("packages.r")

## Code ----
sapply(list.files("./processing_scripts",full.names = T), source)


#### PLANS ---------------------------------------------------------------------

## Load data plan ----
loadData <- drake_plan(
  rawMtbs = load_mtbs(),
  rawPlants = load_plants(),
  rawASVs = load_endophytes(),
  rawPhylogeny = build_phylogeny(rawPlants = rawPlants),
  rawClimate = load_climate(rawPlants = rawPlants)
)

## Process data plan ----
processData <- drake_plan(
  # metabolite data
  mtbsCross = process_mtbs(
    rawMtbs = rawMtbs,
    rawPlants = rawPlants,
    study = "cross"
  ),
  mtbsRubra = process_mtbs(
    rawMtbs = rawMtbs,
    rawPlants = rawPlants,
    study = "rubra"
  ),
  # ASV data
  asvsCross = process_asvs(
    rawASVs = rawASVs,
    rawPlants = rawPlants,
    study = "cross"
  ),
  asvsRubra = process_asvs(
    rawASVs = rawASVs,
    rawPlants = rawPlants,
    study = "rubra"
  ),
  # phylogeny
  phylogeny = process_phylo(rawPhylogeny = rawPhylogeny)
)

## Collate data plan ----
collateData <- drake_plan(
  rubra = collate_rubra(
    rawPlants = rawPlants,
    asvsRubra = asvsRubra,
    mtbsRubra = mtbsRubra,
    rawClimate = rawClimate
  ),
  cross = collate_cross(
    rawPlants = rawPlants,
    asvsCross = asvsCross,
    mtbsCross = mtbsCross,
    rawClimate = rawClimate,
    phylogeny = phylogeny
  )
)


#### MAKE ----------------------------------------------------------------------

## Assemble plan ----
thePlan <- bind_rows(
  loadData,
  processData,
  collateData
)

## Make ----
vis_drake_graph(thePlan)
make(thePlan)
