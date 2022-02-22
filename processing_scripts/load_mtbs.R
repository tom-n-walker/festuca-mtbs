################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Load metabolite data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

load_mtbs <- function(rawPlants){
  ## Load metabolite data ----
  mtbsRaw <- fread("./data/nadline_full_ms1.csv", data.table = F)
  ## Generate metabolite ID map ----
  mtbsID <- mtbsRaw %>%
    # isolate ID columns, rename ID
    select(contains("row")) %>%
    rename(mtbID = `row ID`, MZ = `row m/z`, RT = `row retention time`) %>%
    # make rowname compatible version
    mutate(mtbID = paste0("mtb_", mtbID))
  ## Clean-up peak area data ----
  mtbsClean <- mtbsRaw %>%
    # remove redundant rows
    select(-V141) %>%
    select(!contains("row")) %>%
    # transpose and make data frame
    t %>%
    as.data.frame %>%
    # create column names from mtb IDs
    rename_with(.fn = function(x){x <- mtbsID$mtbID}) %>%
    # reformat sample IDs
    rownames_to_column("sampleID") %>%
    mutate(sampleID = gsub("201222_PFC_metabolomics_DDA_pos_sample_", "", sampleID)) %>%
    mutate(sampleID = gsub("201221_PFC_metabolomics_DDA_pos_sample_", "", sampleID)) %>%
    mutate(sampleID = gsub(".mzXML.Peak.height", "", sampleID)) %>%
    # remove final sample IDs
    filter(!grepl("PFC", sampleID))
  ## Collate data and return ----
  allData <- list(
    mtb_table = mtbsID,
    feature_table = mtbsClean
  )
  return(allData)
}
