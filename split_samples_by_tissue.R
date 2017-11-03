#library(data.table)
library(stringr)
options(echo=TRUE)
# read in SRA run IDs to sample IDs
sra <- read.csv("/SAN/vyplab/HuRNASeq/GTEX_public/gtex_condensed_run_table.csv", header =TRUE, stringsAsFactors = FALSE)

# read in GTEX sample table
samples <- read.table("/SAN/vyplab/HuRNASeq/GTEX/GTEx_Data_V6_Annotations_SampleAttributesDS.txt", header =TRUE, sep = "\t", stringsAsFactors =FALSE) 


samples$sra_run <- sra$run_id[ match( samples$SAMPID, sra$sample_id) ]

# remove samples with no matching SRA run 

samples <- samples[ !is.na(samples$sra) & samples$SMTS != "" ,]

tissues <- unique(samples$SMTSD)
tissues <- tissues[ tissues != "" ]

tissue_format <- function(x){
	x <- str_to_lower(gsub( "\\(.*\\)", "", gsub(" ","_",x)))
	x <- gsub("_-_", "_", x)
	x <- gsub("_$","", x)
	return(x)
}
tissues_formatted <- tissue_format(tissues)
print(tissues_formatted)

tissue_list <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/GTEX/tissue_list.txt"

writeLines(tissues_formatted, tissue_list)

#quit()

# for each tissue, create folder and write list of junctions

outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/GTEX"

filename <- "/SAN/vyplab/HuRNASeq/GTEX/2016-11-11-fullrun1/SRA_ID.leafcutter_op"

for(tis in tissues){
	print(tis)
	tis_format <- tissue_format(tis)
	tisFolder <- paste0(outFolder,"/",tis_format)
	if( ! dir.exists(tisFolder) ){ dir.create(tisFolder) }

	junction_list <- paste0(tisFolder,"/junction_list.txt")

	accessions <- samples[ samples$SMTSD == tis, ]$sra_run

	accessions <- sapply(accessions, FUN = function(x){
				gsub( "SRA_ID", x, filename )} )

	accessions <- accessions[ file.exists(accessions) ]

	print( paste0(length(accessions), " samples found") )
	writeLines( accessions, junction_list )
}	

quit()

