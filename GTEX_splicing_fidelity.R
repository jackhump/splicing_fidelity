library(stringr)
library(dplyr)
library(purrr)
library(data.table)
options(echo = TRUE)
source("extract_gtex_data.R")
source("splicing_fidelity.R")

hg19_introns <- "gencode_hg19_all_introns.bed.gz"
intron_db <- prepareIntrons(intronList=hg19_introns)


# for each tissue:
all_tissues <- unique(age_ranges$tissue)
for( tissue in all_tissues){
  print(tissue)
  # create metadata table for each available sample
  # with get_tissue_data
  metadata <- get_tissue_data(tissue)
  # annotate junctions
  annotated <- select(metadata, junc_file, accession) %>%
    rename( file = junc_file, file_id = accession) %>%
    purrr::pmap(annotateJunctions, intron_db = intron_db)
  # create proportions
  annotated <- purrr::map(annotated, createSimpleProportions)
  
  # get out proportions
  proportions <- purrr::map_df( annotated, "proportions")
  
  proportions <- metadata %>% left_join(proportions, by = "accession" )
  
  # write out proportion table
  tissueName <- str_to_lower(tissue) %>% 
    gsub("\\(|\\)|- ", "",.) %>%
    gsub(" ", "_", .) 
  outFile <- paste0("tissue_splicing_fidelity/", tissueName, "_splicing_fidelity.tab")
  write.table(proportions, "RBPs/all_RBPs_fidelity.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  save(annotated, file = gsub(outFile, ".tab", ".Rdata"))
}

