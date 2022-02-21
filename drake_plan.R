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
loadDataPlan <- drake_plan(
  rawMtbs = load_mtbs(),
  rawTraits = load_plants()
)

#### MAKE ----------------------------------------------------------------------

## Assemble plan ----
thePlan <- bind_rows(
  loadDataPlan
)

## Make ----
make(thePlan)
