# ENCODE RBPs
library(data.table)
library(dplyr)
library(purrr)
setwd("/Users/Jack/SAN/HuRNASeq/GTEX")

source("splicing_fidelity.R")

hg38_introns <- "gencode_hg38_all_introns.bed.gz"

# test on hnRNPK first
files <- "RBPs/hnrnpk_files.txt"

files <- "RBPs/ALS_RBP_list.txt"
files <- read.table(files, header=TRUE, stringsAsFactors=FALSE)

#test <- "../ENCODE/HNRNPK/HepG2_ENCSR853ZJS/processed/control_1/control_1SJ.out.tab"

intron_db <- prepareIntrons(intronList=hg38_introns)

#anno <- annotateJunctions(files[1], intron_db, "control1")

#file_args <- list(file = files, file_id = files )

anno <- files[,1:2] %>%
  purrr::pmap(annotateJunctions, intron_db = intron_db)

anno <- purrr::map(anno, createSimpleProportions)

# get out proportions

proportions <- do.call(rbind, lapply( seq_along(anno), FUN = function(x){ anno[[x]]$proportions }) )

# merge together on file_id

names(proportions)[1] <- "file_id"
proportions <- proportions %>% left_join(files )

write.table(proportions, "RBPs/ALS_RBPs_fidelity.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# for transfer to Rstudio

proportions <- read.table("RBPs/ALS_RBPs_fidelity.tab",header=TRUE, stringsAsFactors=FALSE) 

plots <- proportions %>% 
  split(cell_type) %>%
  map( ~ggplot())
  
ggplot( proportions,
        aes( x = prop_unique, y = prop_sum, colour = condition )) +
  geom_text(aes(label = condition)) +
  facet_wrap(~cell_type)

# give the mean value for the two replicates
proportions_clean <- proportions %>%
  mutate(file_id = gsub("_[0-9]$", "", file_id) ) %>%
  group_by( file_id ) %>%
  dplyr::summarise( prop_unique = mean(prop_unique), 
             prop_sum = mean(prop_sum),
             sd_unique = sd(prop_unique),
             sd_sum = sd(prop_sum),
             condition = first(condition),
             cell_type = first(cell_type)) %>%
  mutate( dataset = gsub("_[^_]+$", "", file_id) )


library(plotly)
proportions_clean %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
          aes( x = prop_unique, y = prop_sum, colour = condition )) +
  #geom_text(aes(label = condition)) +
  geom_line( aes(group = dataset), colour = "black", linetype = 2) +
  geom_point( ) +
  geom_errorbar(aes(ymin = prop_sum - sd_sum,
                    ymax = prop_sum + sd_sum)) + #,
                    #xmin = prop_unique - sd_unique,
                    #xmax = prop_unique + sd_unique))
  facet_wrap(~cell_type)

ggplotly()

