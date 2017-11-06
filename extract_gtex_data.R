# working with GTEX
#setwd("/Users/Jack/SAN/HuRNASeq/GTEX/")
# I created an SQLite database consisting of three tables
	# patients - data on all 570 humans used by the GTEX project
  # samples - each of the 12,000 gtex samples extracted from the 570 patients
  # run_table - the SRA ids that link each GTEX sample to the RNA-seq junctions
library(RSQLite)
library(stringr)
library(dplyr)

junction_dir <- "/Users/Jack/SAN/HuRNASeq/GTEX/2016-11-11-fullrun1/"
junction_dir <- "/SAN/vyplab/HuRNASeq/GTEX/2016-11-11-fullrun1/"

junction_suffix <- ".leafcutter_op"

#support_out <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/GTEX_FC_AGE/leafcutter_support.tab"

# extract data from SQL database

con <- dbConnect(drv=RSQLite::SQLite(), dbname="GTEX.db")

tables <- dbListTables(con)

lDataFrames <- vector("list", length=length(tables))

for (i in seq(along=tables)) {
  lDataFrames[[i]] <- dbGetQuery(conn=con, statement=paste("SELECT * FROM '", tables[[i]], "'", sep=""))
}
# create objects named after the variable names
for( i in 1:length(tables) ){
  assign( tables[i], lDataFrames[[i]] )
}

samples$AGE <- patients$AGE[ match( 
  str_split_fixed(samples$SAMPID, "-",6 )[,2],
  str_split_fixed(patients$SUBJID, "-", 2)[,2] 
) ]

age_ranges <- dplyr::group_by( samples, SMTSD, AGE) %>%
  summarise( n = n() ) %>%
  rename(tissue = SMTSD, age = AGE) 

# ridge plot of proportions of ages in each tissue
# library(ggridges)
# ggplot(
#   age_ranges, 
#   aes( x = age, y = tissue,height = n) ) + 
#   geom_point( aes(size = n))


# all Frontal Cortex samples - what are the age distributions?
#fc <- samples[ samples$SMTSD == "Brain - Frontal Cortex (BA9)", ]

# samples$SAMPID consists of GTEX-patients$SUBJID-other-stuff
# ages <- patients$AGE[ match( 
#         str_split_fixed(fc$SAMPID, "-",6 )[,2],
#         str_split_fixed(patients$SUBJID, "-", 2)[,2] 
#         ) ]
# 
# gender <- patients$GENDER[ match( 
#         str_split_fixed(fc$SAMPID, "-",6 )[,2],
#         str_split_fixed(patients$SUBJID, "-", 2)[,2] 
#         ) ]

# Create metadata of patients for a particular tissue

get_tissue_data <- function( tissue_type ){
  # subset the sample table
  df <- samples[ samples$SMTSD == tissue_type, ]
  # match the df with the patients table
  tissue_data <- patients[ match( 
    str_split_fixed(df$SAMPID, "-",6 )[,2],
    str_split_fixed(patients$SUBJID, "-", 2)[,2] ), ]
  # add columns
  tissue_data$SAMPID <- df$SAMPID[ match(
    str_split_fixed(tissue_data$SUBJID, "-",6 )[,2],
    str_split_fixed(df$SAMPID, "-", 6)[,2] 
  )]
  tissue_data$accession <- run_table$run_id[ match( tissue_data$SAMPID, run_table$sample_id)]
  tissue_data$junc_file <- paste0( junction_dir, tissue_data$accession, junction_suffix)
  
  # remove_duplicate rows and files that don't exist
  tissue_data <- distinct(tissue_data) %>%
                 filter( file.exists( junc_file) ) %>%
                 mutate( tissue = tissue_type ) %>%
                 select( tissue, everything() )
  
  return(tissue_data)
}


# 
# 
# quit()
# 
# # for each tissue:
# all_tissues <- unique(age_ranges$tissue)
# for( tissue in all_tissues){
#   # create metadata table for each available sample
#   # with get_tissue_data
#   metadata <- get_tissue_data(tissue)
#   annotated <- select(metadata, junc_file, accession) %>%
#     purrr::pmap(annotateJunctions, intron_db = intron_db)
#   
# 
# # annotate junctions for that tissue
# 
# # create proportions
# 
# # save the complete object as an Rdata object
# }
# 
# heart <- filter(age_ranges, tissue == "Heart - Left Ventricle")
# sum(heart$n)
# 
# 
# # compare brain to heart
# 
# heart_df <- get_tissue_data("Heart - Left Ventricle")
# 
# # annotate heart junctions
# heart_annotated <- lapply( heart_df$junc_file[1:20], FUN = function(x){
#   annotateJunctions(file = x, intron_db)
# })
# # create simple proportions of novel to annotated - by number of junctions and by sums of junctions
# heart_prop <- createSimpleProportions(heart_annotated)
# 
# 
# tibial_df <- get_tissue_data("Nerve - Tibial")
# 
# tibial_annotated <- lapply( 1:length(tibial_df$junc_file), FUN = function(x){
#   annotateJunctions(file = tibial_df$junc_file[x], intron_db, file_id = tibial_df$accession[x])
# })
# 
# tibial_prop <- createSimpleProportions(tibial_annotated)
# 
# tibial_df <- left_join(tibial_df, tibial_prop, by = "accession") %>%
#              mutate( prop_unique = as.numeric(as.character(prop_unique)),
#                      prop_sum = as.numeric(as.character(prop_sum))
#                      )
# 
# ggplot( tibial_df, aes( x = AGE, y = as.numeric(prop_sum), group = AGE)) + geom_jitter() + geom_boxplot()
# ggplot( tibial_df, aes( x = AGE, y = as.numeric(prop_unique), group = AGE)) + geom_boxplot()
# 
# kruskal.test(x = tibial_df$prop_sum, g = as.factor(tibial_df$AGE)  )
# 
# testes_df <- get_tissue_data("Testis")
# testes_annotated <- lapply( testes_df$junc_file[1:20], FUN = function(x){
#   annotateJunctions( file = x, intron_db)
# })
# testes_prop <- createSimpleProportions(testes_annotated)
# 
# vagina_df <- get_tissue_data("Vagina")
# vagina_annotated <- lapply( vagina_df$junc_file[1:20], FUN = function(x){
#   annotateJunctions( file = x, intron_db)
# })
# vagina_prop <- createSimpleProportions(vagina_annotated)
# 
# 
# compare <- rbind(heart_prop, tibial_prop) %>%
#            rbind( testes_prop) %>%
#            rbind( vagina_prop) %>%
#            mutate( 
#              tissue = c( 
#              rep( "heart", nrow(heart_prop) ),
#              rep("tibial nerve", nrow(tibial_prop) ),
#              rep("testes", nrow(testes_prop)),
#              rep("vagina", nrow(vagina_prop))
#              ) )
# 
# ggplot( compare[ compare$tissue != "testes",], aes( x = prop_unique, y = prop_sum, colour = tissue )) + geom_point()
# 
# 
# fc_patients <- patients[ match( 
#         str_split_fixed(fc$SAMPID, "-",6 )[,2],
#         str_split_fixed(patients$SUBJID, "-", 2)[,2] 
#         ), ]
# 
# 
# 
# 
# 
# # the DTHHRDY is a integer scale encoding of the death circumstance
#     # "Death classification based on the 4-point Hardy Scale:
#     # 1) Violent and fast death Deaths due to accident, blunt force trauma or suicide, terminal phase estimated at < 10 min. 
#     # 2) Fast death of natural causes Sudden unexpected deaths of people who had been reasonably healthy, after a terminal phase estimated at < 1 hr (with sudden death from a myocardial infarction as a model cause of death for this category) 
#     # 3) Intermediate death Death after a terminal phase of 1 to 24 hrs (not classifiable as 2 or 4); patients who were ill but death was unexpected 
#     # 4) Slow death Death after a long illness, with a terminal phase longer than 1 day (commonly cancer or chronic pulmonary disease); deaths that are not unexpected 
#     # 0) Ventilator Case All cases on a ventilator immediately before death."
# 
# 
# # dig out the neuronal proportions estimation from CIBERSORT - do any of these factors explain the variance in neuronal proportion?
# # hypothesis - slow deaths and older age might predict lower neuronal proportions in the cortex
# 
# neuronal <- read.csv("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/Results/CIBERSORT.GTEX_frontal_cortex_Darmanis_single_cell_signature.csv",header=TRUE)
# 
# fc_patients$neuronal_prop <- neuronal$fpkm_neurons[ match(
#         str_split_fixed(fc_patients$SUBJID, "-",6 )[,2],
#         str_split_fixed(neuronal$Input.Sample, "-", 6)[,2] 
#         ) ]
# 
# fc_patients$microglial_prop <- neuronal$fpkm_microglia[ match(
#         str_split_fixed(fc_patients$SUBJID, "-",6 )[,2],
#         str_split_fixed(neuronal$Input.Sample, "-", 6)[,2] 
#         ) ]
# 
# # group samples by young (20 - 59) and old (60-79) - fairly arbitrary but it gives similar sample numbers on each side.
# 
# ages <- c("20-29", "30-39","40-49", "50-59", "60-69", "70-79")
# 
# fc_patients$age_bin <- ifelse( fc_patients$AGE %in% ages[1:3],
#       "YOUNG", "OLD")
# 
# 
# # fish out the full sample ID from fc and add an age bin
# fc_patients$SAMPID <- fc$SAMPID[ match(
#   str_split_fixed(fc_patients$SUBJID, "-",6 )[,2],
#   str_split_fixed(fc$SAMPID, "-", 6)[,2] 
#   )]
# 
# fc_patients$accession <- run_table$run_id[ match( fc_patients$SAMPID, run_table$sample_id)]
# 
# fc_patients$junc_file <- paste0( junction_dir, fc_patients$accession, junction_suffix)
# 
# table(file.exists( fc_patients$junc_file ))
# 
# # create support file for leafcutter
# 
# out_table <- select( fc_patients, junc_file, age_bin, GENDER, DTHHRDY) %>% arrange( desc(age_bin) )
# 
# write.table( out_table, support_out, col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE )
# 
