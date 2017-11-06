# ENCODE RBPs
library(data.table)
library(dplyr)
library(purrr)
library(ggrepel)
setwd("/Users/Jack/SAN/HuRNASeq/GTEX")

source("splicing_fidelity.R")

hg38_introns <- "gencode_hg38_all_introns.bed.gz"

# test on hnRNPK first
files <- "RBPs/hnrnpk_files.txt"

files <- "RBPs/ALS_RBP_list.txt"
files <- "RBPs/all_ENCODE_junctions.txt"

files <- read.table(files, header=TRUE, stringsAsFactors=FALSE)

#test <- "../ENCODE/HNRNPK/HepG2_ENCSR853ZJS/processed/control_1/control_1SJ.out.tab"

intron_db <- prepareIntrons(intronList=hg38_introns)

#anno <- annotateJunctions(files[1], intron_db, "control1")

#file_args <- list(file = files, file_id = files )

anno <- files[,1:2] %>%
  purrr::pmap(annotateJunctions, intron_db = intron_db)

anno <- purrr::map(anno, createSimpleProportions)

# get out counts
counts <- purrr::map_df(anno, "counts")

# get out proportions
proportions <- purrr::map_df( anno, "proportions")

# merge together on file_id

names(proportions)[1] <- "file_id"
proportions <- proportions %>% left_join(files )

write.table(proportions, "RBPs/all_RBPs_fidelity.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# for transfer to Rstudio

#proportions2 <- read.table("RBPs/ALS_RBPs_fidelity.tab",header=TRUE, stringsAsFactors=FALSE) 

proportions <- read.table("RBPs/all_RBPs_fidelity.tab",header=TRUE, stringsAsFactors=FALSE) 

# first test plot
# ggplot( proportions,
#         aes( x = prop_unique, y = prop_sum, colour = condition )) +
#   geom_text(aes(label = condition)) +
#   facet_wrap(~cell_type)

# give the mean value for the two replicates
proportions_clean <- proportions %>%
  mutate(file_id = gsub("_[0-9]$", "", file_id) ) %>%
  group_by( file_id ) %>%
  dplyr::summarise(
             prop_unique = mean(prop_unique), 
             prop_sum = mean(prop_sum),
             prop_unique_anchored = mean(prop_unique_anchored),
             prop_unique_unanchored = mean(prop_unique_unanchored),
             prop_unique_skiptic = mean(prop_unique_skiptic),
             prop_sum_anchored = mean(prop_sum_anchored),
             prop_sum_unanchored = mean(prop_sum_unanchored),
             prop_sum_skiptic = mean(prop_sum_skiptic),
             condition = first(condition),
             cell_type = first(cell_type)) %>%
  mutate( dataset = gsub("_[^_]+$", "", file_id) ) %>%
  arrange(desc(prop_sum), desc(prop_unique))



proportions_clean %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
          aes( x = prop_unique, y = prop_sum, colour = condition )) +
  #geom_text(aes(label = condition)) +
  geom_line( aes(group = dataset), colour = "black", linetype = 2) +
  #geom_point( ) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) +
  geom_text( aes(label = condition), size = 2) +
  theme_bw() +
  ylab("Proportion of total novel junctions") +
  xlab("Proportion of unique novel junctions") +
  labs(title = "Each knockdown linked to its control")
  # geom_text_repel( data = head(proportions_clean, 40),
  #                  aes( x = prop_unique, y = prop_sum,
  #                       label = condition),
  #                  colour = "black", size = 2)

ggsave("RBPs/all_RBPs_scatter_fidelity.pdf")

# explore the plot with plotly!
library(plotly)
ggplotly()


# can I link each knockdown to its control and plot the relative change instead?

relative <- proportions_clean %>%
                 mutate(experiment = str_split_fixed(file_id, "_", 4)[,3]) %>%
                 split( .$condition == "control" )    

relative <- left_join(relative[[1]],
                      relative[[2]], 
                      by = c("experiment","cell_type", "dataset"),
                      suffix = c("_kd", "_ctl")) %>%
            mutate( 
              rel_prop_unique = prop_unique_kd - prop_unique_ctl,
              rel_prop_sum = prop_sum_kd - prop_sum_ctl,
              rel_sum_anchored = prop_sum_anchored_kd - prop_sum_anchored_ctl,
              rel_sum_unanchored = prop_sum_unanchored_kd - prop_sum_unanchored_ctl,
              rel_sum_skiptic = prop_sum_skiptic_kd - prop_sum_skiptic_ctl,
              rel_unique_anchored = prop_unique_anchored_kd - prop_unique_anchored_ctl,
              rel_unique_unanchored = prop_unique_unanchored_kd - prop_unique_unanchored_ctl,
              rel_unique_skiptic = prop_unique_skiptic_kd - prop_unique_skiptic_ctl,
              condition = condition_kd)

relative %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_prop_unique, y = rel_prop_sum, colour = condition )) +
  geom_text( aes(label = condition_kd), size = 2) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) + 
  theme_bw() +
  ylab("Proportion of total novel junctions\nrelative to control") +
  xlab("Proportion of unique novel junctions\nrelative to control") +
  labs(title = "All ENCODE RBPs")

ggsave("RBPs/all_RBPs_scatter_fidelity_relative.pdf")


# match in the specific control - is there a batch effect in the relative directions?
meta <- read.table("RBPs/all_ENCODE_metadata.txt",header=TRUE, sep = "\t", stringsAsFactors = FALSE)

relative_controls <- meta %>%
          select( Accession, Controls) %>%
          left_join(x = relative, y = ., by = c( "experiment" = "Accession"))
# colour by control
# relative_controls %>%
#   #filter( condition != "HNRNPC" ) %>%
#   ggplot(
#     aes( x = rel_prop_unique, y = rel_prop_sum, colour = Controls )) +
#   #geom_text(aes(label = condition)) +
#   #geom_line( aes(group = dataset), colour = "black", linetype = 2) +
#   geom_point( aes(alpha = 0.2) ) +
#   guides(colour = FALSE) +
#   facet_wrap(~cell_type) + 
#   theme_bw()

top_controls <- group_by(relative_controls, Controls) %>%
                summarise( n = n() ) %>%
                arrange(desc(n)) %>%
                head(20)

relative_controls %>%
  filter( Controls %in% top_controls$Controls) %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_prop_unique, y = rel_prop_sum, colour = Controls )) +
  #geom_text(aes(label = condition)) +
  #geom_line( aes(group = dataset), colour = "black", linetype = 2) +
  #geom_point( aes(alpha = 0.2) ) +
  geom_text( aes( label = condition_kd), size = 2) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) + 
  theme_bw() +
  ylab("Proportion of total novel junctions\nrelative to control") +
  xlab("Proportion of unique novel junctions\nrelative to control") +
  labs(title = "Coloured by batch of control samples")

ggsave("RBPs/all_RBPs_scatter_fidelity_relative_batch.pdf")

ggplotly()

# separate by different classes of novel junction

proportions_clean %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = prop_unique_anchored, y = prop_sum_anchored, colour = condition )) +
  #geom_text(aes(label = condition)) +
  geom_line( aes(group = dataset), colour = "gray", linetype = 2) +
  #geom_point( ) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) +
  geom_text( aes(label = condition), size = 2) +
  theme_bw() +
  ylab("Proportion of total anchored novel junctions") +
  xlab("Proportion of anchored novel junctions") +
  labs(title = "Anchored novel")

proportions_clean %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = prop_unique_unanchored, y = prop_sum_unanchored, colour = condition )) +
  #geom_text(aes(label = condition)) +
  geom_line( aes(group = dataset), colour = "gray", linetype = 2) +
  #geom_point( ) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) +
  geom_text( aes(label = condition), size = 2) +
  theme_bw() +
  ylab("Proportion of total novel junctions") +
  xlab("Proportion of unanchored novel junctions") +
  labs(title = "Unanchored novel")

proportions_clean %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = prop_unique_skiptic, y = prop_unique_sum, colour = condition )) +
  #geom_text(aes(label = condition)) +
  geom_line( aes(group = dataset), colour = "gray", linetype = 2) +
  #geom_point( ) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) +
  geom_text( aes(label = condition), size = 2) +
  theme_bw() +
  ylab("Proportion of total skiptic junctions") +
  xlab("Proportion of skiptic junctions") +
  labs(title = "Skiptic")

 # and the relative plots:

rel_anchored <-relative %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_unique_anchored, y = rel_sum_anchored, colour = condition )) +
  geom_text( aes(label = condition_kd), size = 2) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) + 
  theme_bw() +
  ylab("Proportion of total anchored novel junctions\nrelative to control") +
  xlab("Proportion of unique anchored novel junctions\nrelative to control") +
  labs(title = "Relative anchored novel")

rel_unanchored <- relative %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_unique_unanchored, y = rel_sum_unanchored, colour = condition )) +
  geom_text( aes(label = condition_kd), size = 2) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) + 
  theme_bw() +
  ylab("Proportion of total unanchored novel junctions\nrelative to control") +
  xlab("Proportion of unique unanchored novel junctions\nrelative to control") +
  labs(title = "Relative unanchored ")

rel_skiptic <- relative %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_unique_skiptic, y = rel_sum_skiptic, colour = condition )) +
  geom_text( aes(label = condition_kd), size = 2) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) + 
  theme_bw() +
  ylab("Proportion of total skiptic junctions\nrelative to control") +
  xlab("Proportion of unique skiptic junctions\nrelative to control") +
  labs(title = "Relative Skiptic")

ggplotly()

pdf("RBPs/all_RBPs_scatter_fidelity_relative_by_type.pdf")
print(rel_anchored)
print(rel_unanchored)
print(rel_skiptic)
dev.off()
  