library(dplyr)
library(stringr)
library(data.table)
# for a given junction file:
#setwd("/Users/Jack/SAN")
#setwd("/SAN/vyplab/HuRNASeq/GTEX/")
#file <- "HuRNASeq/GTEX/2016-11-11-fullrun1/SRR1413562.leafcutter_op"
#introns <- "gencode_hg19_all_introns.bed.gz"
#introns <- "/Users/Jack/google_drive/Work/PhD_Year_3/leafcutter/leafviz/annotation_codes/gencode_hg19/gencode_hg19_all_introns.bed.gz"

# read in intron data
prepareIntrons <- function(intronList){
  introns <- fread(paste( "zcat < ", intronList) , data.table=FALSE)
  intron_db <- introns %>%
    rename(chr = V1, 
         start = V2, 
         end = V3, 
         gene = V4, 
         EnsemblID = V5, 
         strand = V6, 
         transcript = V7, 
         dot = V8, 
         transcript_type = V9,
         transcript_annotation = V10) %>%
    select( -dot ) %>%
    mutate( end = end - 1 ) %>% # due to 0-base/1-base shenanigans
    filter( !duplicated( paste(chr, start, end) ) ) # remove duplicate entries
  return(intron_db)
}


annotateJunctions <- function(file, intron_db, file_id){
  # function to annotate a list of junctions according to a set of introns
  # presumably from GENCODE
  
  if( !file.exists(file) ){
    message("file doesnt exist")
    return(NULL)
  }
  
  junctions <- fread(file, data.table=FALSE)
  
  # work out whether a leafcutter junction file or a STAR junction file
  fileType <- "unknown"
  if( 
    ncol(junctions) == 9 &
    all( sapply(junctions[,2:ncol(junctions)], is.integer) ) 
    ){
    fileType <- "STAR"
    
    sorted <- junctions %>% 
      rename( chr = V1, 
              start = V2, 
              end = V3, 
              strand = V4, 
              intron_motif = V5, 
              annotated = V6,
              count = V7,
              multiCount = V8,
              maxOverhang = V9) %>%
      dplyr::arrange(chr, start) %>%
      select( chr, start, end, count, strand ) %>%
      filter( strand != 0) %>%
      mutate( start = start - 1, 
              strand = ifelse( strand == 1, "+", "-" ))
  }
  if(ncol(junctions) == 6){
  fileType <- "leafcutter"
  sorted <- junctions %>% 
    dplyr::arrange(V1, V2) %>%
    mutate( V1 = paste0("chr", V1)) %>%
    rename( chr = V1, start = V2, end = V3, dot = V4, count = V5, strand = V6) %>%
    select( chr, start, end, count, strand )
  }
  
  if( fileType == "unknown"){
    message("file not recognised")
    return(NULL)
  }
  
  # find exact matches of both splice sites
  intersect_both <- sorted %>%
    left_join(intron_db, by = c("chr","start", "end", "strand")) %>%
    select( chr, start,end, strand, count, transcript )
  
  # get just annotated junctions out
  annotated <- intersect_both %>%
    filter( !is.na(transcript) ) %>%
    mutate( type = "annotated") %>%
    select( chr, start, end, strand, count, type)
  
  # by definition any junction that can't be found in the GENCODE intron table is novel
  novel_junctions <- filter(intersect_both, is.na(transcript)) %>%
    select( chr, start, end, strand, count)
  
  
  
  # get out different types of novel junction
  # skiptics - both ends are annotated but separately
  # anchored cryptics - only one end is annotated
  # cryptic_unanchored cryptics - neither end are annotated 
  
  # semi join only keeps rows in X that match in Y
  # so only keep novel junctions where start and end match separately
  
  skiptic <- novel_junctions %>%
    semi_join( intron_db, by = c("chr", "start", "strand") ) %>%
    semi_join( intron_db, by = c( "chr", "end", "strand" ) ) %>%
    arrange( chr,start )  %>%
    mutate( type = "skiptic")
  
  anchored_start <- novel_junctions %>%
    semi_join( intron_db, by = c("chr", "start", "strand") ) 
  
  anchored_end <- novel_junctions %>%
    semi_join( intron_db, by = c("chr", "end", "strand") ) 
  
  cryptic_anchored <- rbind( anchored_start, anchored_end) %>%
    arrange( chr,start ) %>%
    anti_join( skiptic, by = c("chr", "start", "end", "count", "strand") ) %>%
    filter( !duplicated( paste(chr, start, end) ) )  %>%
    mutate( type = "cryptic_anchored")
  
  cryptic_unanchored <- novel_junctions %>%
    anti_join(skiptic, by = c("chr", "start", "end", "count", "strand")) %>%
    anti_join(cryptic_anchored, by = c("chr", "start", "end", "count", "strand")) %>%
    mutate( type = "cryptic_unanchored")
  
  # bind all together
  all_junctions <- rbind( annotated, skiptic, cryptic_unanchored, cryptic_anchored ) %>%
    arrange(chr,start) 
  
  summary <- group_by(all_junctions, type) %>%
    summarise( n_unique = n(), sum_counts = sum(count)) %>%
    mutate( prop_unique = n_unique / sum(n_unique),
            prop_counts = sum_counts / sum(sum_counts))
  
  #SRR_code <- gsub(".leafcutter_op", "", basename(file) )
  
  # test that everything worked
  if( 
    ( nrow(cryptic_unanchored) + nrow(cryptic_anchored) + nrow(skiptic) == nrow(novel_junctions) ) &
    nrow(all_junctions) == nrow(intersect_both) ){
    return( 
      list(ID = file_id,
           counts = summary,
           all = as.tbl(all_junctions)
           )
      )
  }else{
    message("Error! The sums don't add up")
    return("ERROR")
  }

}


# system.time({
#   annotateJunctions( file, introns )
# })
# 
# # test - all ~140 frontal cortex samples
# 
# system.time({
#   frontal_cortex_annotated <- lapply( fc_patients$junc_file, FUN = function(x){
#     annotateJunctions(file = x, intron_db)
#   })
# })

# function to create simple proportion of annotated vs non annotated
createSimpleProportions <- function( annotated_df ){
    sample <- annotated_df$counts
    ID <- annotated_df$ID
    prop_unique_anchored <- sample$n_unique[2] / sum(sample$n_unique)
    prop_unique_unanchored <- sample$n_unique[3] / sum(sample$n_unique)
    prop_unique_skiptic <- sample$n_unique[4] / sum(sample$n_unique)
    prop_unique <- 1 - (sample$n_unique[1] / sum(sample$n_unique) )
    prop_sum <- 1 - (sample$sum_counts[1] / sum(sample$sum_counts) )
    prop_sum_anchored <- sample$sum_counts[2] / sum(sample$sum_counts)
    prop_sum_unanchored <- sample$sum_counts[3] / sum(sample$sum_counts)
    prop_sum_skiptic <- sample$sum_counts[4] / sum(sample$sum_counts)
    annotated_df$proportions <- data.frame(accession = ID,
                                           prop_unique = prop_unique,
                                           prop_sum = prop_sum,
                                           prop_unique_anchored = prop_unique_anchored,
                                           prop_unique_unanchored = prop_unique_unanchored,
                                           prop_unique_skiptic = prop_unique_skiptic,
                                           prop_sum_anchored = prop_sum_anchored,
                                           prop_sum_unanchored = prop_sum_unanchored,
                                           prop_sum_skiptic = prop_sum_skiptic)
  return(annotated_df)
}

#with(all_props, plot( prop_unique, prop_sum) )
# 
# chr_sizes <- intron_db %>%
#         group_by(chr) %>%
#         summarise( end = max(end) ) 
# 
# 
# # what is the distribution of novel junctions per chromosome?
# chr_sizes <- intron_db %>%
#   group_by(chr) %>%
#   summarise( end = max(end) ) 
# 
# novel_junctions_per_chr <- novel_junctions %>% 
#         group_by( chr ) %>%
#         summarise( count = sum(count) ) %>%
#         left_join( chr_sizes, by = "chr" ) %>%
#         filter( !grepl( "GL|MT", chr) ) %>%
#         mutate(prop = (count /  end) * 1E6,
#                chr = factor(chr, levels = paste0("chr", c(1:22,"X","Y") ) ) ) %>%
#         arrange( chr )



