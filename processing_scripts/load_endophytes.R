################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Load endophyte data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

load_endophytes <- function(){
  ## Load raw data ----
  dataRaw <- fread(
    "./data/Festuca_endophytes_Illumina_ITS_1285332544.abundant_ASVs.tsv",
    data.table = F
  )
  ## Format data ----
  # generate clean dataset
  dataClean <- dataRaw %>%
    # select samples only
    select(contains("Festuca")) %>%
    # transpose and make data frame
    t %>%
    as.data.frame %>%
    # rename with ASV identifiers
    rename_with(function(x){gsub("V", "asv_", x)}) %>%
    # make plant ID column by splitting IDs
    rownames_to_column("plantID") %>%
    mutate(
      plantID = do.call(
        c, 
        lapply(
          str_split(plantID, "_"), 
          function(x) as.numeric(x[3])
        )
      )
    ) %>%
    # arrange by plantID
    arrange(plantID)
  ## Return ----
  return(dataClean)
}
