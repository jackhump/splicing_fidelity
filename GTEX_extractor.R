library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
options(echo=T)
iFolder <- "/SAN/vyplab/HuRNASeq/opthalmology_work"
plotFolder <- "/SAN/vyplab/HuRNASeq/opthalmology_work/plots/"
if(!dir.exists(plotFolder)){dir.create(plotFolder) }

outFolder <- "/SAN/vyplab/HuRNASeq/GTEX/"
if(!dir.exists(outFolder)){dir.create(outFolder) }
organ <- "Brain"
organFolder <- paste0(outFolder,organ)
if(!dir.exists(organFolder)){dir.create(organFolder) }



annotation_file="/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab"
gtex_samples=paste(iFolder, "GTEx_Data_V6_Annotations_SampleAttributesDS.txt",sep="/")
gtex_subjects=paste(iFolder, "GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", sep = "/")
# GTEx data 
gsample <- fread(gtex_samples, header=T )
gtex <- fread(gtex_subjects, header=T )

# get the sample IDs
gtex_brain_samples <- subset(gsample, gsample$SMTS == "Brain" )$SAMPID

# subset out just those subjects from the main gtex database
gtex_brain <- gtex[, which(names(gtex) %in% c("Name", "Description",gtex_brain_samples) ) ]
# using data.table syntax
gtex_brain <- gtex[, gtex_brain_samples %in% names(gtex), with = FALSE ]

# get the brain region for each sample
gtex_brain_region <- subset(gsample, gsample$SMTS == "Brain" )$SMTSD
# remove "Brain - " from each entry
gtex_brain_region <- str_sub( gtex_brain_region, start = 9)
# remove bracketed part of region
gtex_brain_region <- gsub(" \\(.*\\)", "", gtex_brain_region)

# separate by region
for( region in unique(gtex_brain_region) ){
  region_samples <- gtex_brain_samples[ gtex_brain_region == region ]
  # extract the region specific samples - switch to data.frame 
  data <- as.data.frame( gtex_brain[, region_samples %in% names(gtex), with = FALSE ] )
  # set row names and remove the other columns
  row.names(data) <- data$Name
  data[,1:2] <- NULL
  # create summary statistics
  data.mean <- apply(data, MAR = 1, FUN = mean )
  data.sd <- apply(data, MAR = 1, FUN = sd)
  data.sem <- apply(data, MAR = 1, FUN = function(x) sd(x) / sqrt( length(x) ))
  # create new frame with the row names of the old frame with the summary statistics
  data <- data.frame(EnsemblID = row.names(data), mean = data.mean, sd = data.sd, sem = data.sem, row.names = row.names(data))
  # write out this table 
  region.out <- str_to_lower(region)
  region.out <- gsub(" ", "_", region.out)

  table.out <- paste0( organFolder, "/", region.out, "_expression.tab" )

  write.table( data, table.out, sep = "\t", row.names = FALSE, col.names = T, quote = FALSE )
  # append to list of dataframes for each region

# 

}

gtex_brain <- rbind(gtex_brain_region, gtex_brain)

gtex_brain$geneID <- str_split_fixed(gtex_brain$Name, "\\.", 2)[,1]


# compute average FPKM per gene


