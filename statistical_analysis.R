library(tidyverse)
library(ape) # tree based operations
library(survival) # For coin package
library(coin) # For permuation test
library(effsize) # for cohen's d statistic
library(vcfR) # For parwise SNP 
library(cluster)
library(igraph)
library(visNetwork)

setwd("/home/dfgmrtc/NC_ANALYSIS/North_Central_analysis_2023")

# Load Metadata

analysis_meta_info <- read_tsv(file = "./Metadata/metadata_2.tsv")



########################## 1. Distribution of terminal branch lengths ########################################

## Load tree file

north_central <- read.tree("./Tree_data/nc_tree.nhx")

south_west <- read.tree("./Tree_data/sw2_tree.nhx")

north_we <- read.tree("./Tree_data/nwe2_tree.nhx")

## Retrieve TBL data using cophenetic.phylo function

# NC
nc_tbl_pairwise_dist <- cophenetic.phylo(north_central)
nc_tbl_pairwise_dist <- nc_tbl_pairwise_dist[upper.tri(nc_tbl_pairwise_dist)]
nc_tbl_pairwise_dist <- nc_tbl_pairwise_dist[nc_tbl_pairwise_dist > 0]

# SW

sw_tbl_pairwise_dist <- cophenetic.phylo(south_west)
sw_tbl_pairwise_dist <- sw_tbl_pairwise_dist[upper.tri(sw_tbl_pairwise_dist)]
sw_tbl_pairwise_dist <- sw_tbl_pairwise_dist[sw_tbl_pairwise_dist > 0]

# NWE

nwe_tbl_pairwise_dist <- cophenetic.phylo(north_we)
nwe_tbl_pairwise_dist  <- nwe_tbl_pairwise_dist [upper.tri(nwe_tbl_pairwise_dist )]
nwe_tbl_pairwise_dist  <- nwe_tbl_pairwise_dist [nwe_tbl_pairwise_dist > 0]


# Convert tbl data to data_frame
max_length <- max(length(nc_tbl_pairwise_dist), length(sw_tbl_pairwise_dist), length(nwe_tbl_pairwise_dist))

north_central_v <- c(nc_tbl_pairwise_dist, rep(NA, max_length - length(nc_tbl_pairwise_dist)))
south_west_v <- c(sw_tbl_pairwise_dist, rep(NA, max_length - length(sw_tbl_pairwise_dist)))
north_we_v <- c(nwe_tbl_pairwise_dist, rep(NA, max_length - length(nwe_tbl_pairwise_dist)))

tbl_dist_df <- data.frame("north_central" = north_central_v, "south_west" = south_west_v, "north_we" = north_we_v)

# Re-construct data_frame for box plots

boxplot_tbl_dist_df <- tbl_dist_df %>%
  gather("north_central", "south_west", "north_we", key = "region", value = "terminal_branch_length")

# Plot rough box plots to visualize distribution of tbl

ggplot(boxplot_tbl_dist_df, aes(x = region, y = terminal_branch_length )) +
  geom_boxplot()


## Testing for normality
# Visual inspection
par(mfrow = c(1,3))
hist(nc_tbl_pairwise_dist, main = "North_central")
hist(sw_tbl_pairwise_dist, main = "south_west")
hist(nwe_tbl_pairwise_dist, main = "North_east_west_2")

# q-q plot
qqnorm(nc_tbl_pairwise_dist, main = "North_central")
qqline(nc_tbl_pairwise_dist)

qqnorm(sw_tbl_pairwise_dist, main = "south_west")
qqline(sw_tbl_pairwise_dist)

qqnorm(nwe_tbl_pairwise_dist, main = "North_east_west")
qqline(nwe_tbl_pairwise_dist)

# Shapiro test
shapiro.test(nc_tbl_pairwise_dist)
shapiro.test(sw_tbl_pairwise_dist)
shapiro.test(nwe_tbl_pairwise_dist)

# All the tests for normality indicates that the distributions are not normally distributed

# SO check if the distribution of terminal branch lengths are different (Wilcox text sensitive to different observation length)
nc_sw_2 <- wilcox.test(nc_tbl_pairwise_dist, sw_tbl_pairwise_dist)
nc_nwe_2 <- wilcox.test(nc_tbl_pairwise_dist, nwe_tbl_pairwise_dist)
sw_2_nwe_2 <- wilcox.test(sw_tbl_pairwise_dist, nwe_tbl_pairwise_dist)


# Normalize terminal branch lengths to adjust for sample size difference and see if there is still going to be a statistically significant difference
normalize_branch_lengths <- function(terminal_lengths) {
  normalized_lengths <- terminal_lengths / sum(terminal_lengths)  # Normalize by total sum
  return(normalized_lengths)
}

normalized_nc_tbl_pairwise_dist <- normalize_branch_lengths(nc_tbl_pairwise_dist)

normalized_sw_tbl_pairwise_dist <- normalize_branch_lengths(sw_tbl_pairwise_dist)

normalized_nwe_tbl_pairwise_dist <- normalize_branch_lengths(nwe_tbl_pairwise_dist)

nc_sw_n <- wilcox.test(normalized_nc_tbl_pairwise_dist, normalized_sw_tbl_pairwise_dist)
nc_nwe_n <- wilcox.test(nc_tbl_pairwise_dist, nwe_tbl_pairwise_dist)
sw_n_nwe_n <- wilcox.test(sw_tbl_pairwise_dist, nwe_tbl_pairwise_dist)


## Boostrap_resampling: Rather than simply normalizing the  original values, perform a boostrap resampling with 10,000 replicates and subsample to adjust for sample size
bootstrap_resampling <- function(terminal_lengths, size) {
  replicate_distributions <- sample(terminal_lengths, size = size, replace = TRUE)
  return(replicate_distributions)
}


boostrap_nc_tbl_pairwise_dist <- bootstrap_resampling(nc_tbl_pairwise_dist, 10000)
boostrap_sw_tbl_pairwise_dist <- bootstrap_resampling(sw_tbl_pairwise_dist, 10000)
boostrap_nwe_tbl_pairwise_dist <- bootstrap_resampling(nwe_tbl_pairwise_dist, 10000)


# Adjust the sample size of largest terminal branch lengths distribution to match for the smallest
samallest_size <- sum(lengths(boostrap_nwe_tbl_pairwise_dist))
unlisted_nc_bs <- unlist(boostrap_nc_tbl_pairwise_dist)
adjusted_nc <- list(sample(unlisted_nc_bs, samallest_size, replace = TRUE))
adjusted_nc <- unlist(adjusted_nc)

adjusted_sw <- sample(boostrap_sw_tbl_pairwise_dist, samallest_size, replace = TRUE)
adjusted_nwe <- sample(boostrap_nwe_tbl_pairwise_dist, samallest_size, replace = TRUE)

nc_sw_b <- wilcox.test(adjusted_nc, adjusted_sw)
nc_nwe_b <- wilcox.test(adjusted_nc, adjusted_nwe)
sw_nwe_b <- wilcox.test(adjusted_sw, adjusted_nwe)


## Use permutation test that do not rely on the assumption of equal sample sizes. 

## NC vs SW
# Combine the data into a single vector
combined_nc_sw <- c(nc_tbl_pairwise_dist, sw_tbl_pairwise_dist)

# Create group labels indicating the origin of each observation
nc_sw_group_labels <- as.factor(rep(c("NC", "SW"), c(length(nc_tbl_pairwise_dist), length(sw_tbl_pairwise_dist))))

# Perform permutation test for one-way ANOVA
nc_sw_perm_test_result <- oneway_test(combined_nc_sw ~ nc_sw_group_labels, distribution = approximate(nresample = 10000))

## Calculate effect size for the permutation test to quantify the magnitude of differences between the terminal branch length distributions of the trees

# Calculate the effect size (Cohen's d)
nc_sw_effect_size <- cohen.d(combined_nc_sw ~ nc_sw_group_labels)

## NC vs NWE
combined_nc_nwe <- c(nc_tbl_pairwise_dist, nwe_tbl_pairwise_dist)
nc_nwe_group_labels <- as.factor(rep(c("NC", "NWE"), c(length(nc_tbl_pairwise_dist), length(nwe_tbl_pairwise_dist))))
nc_nwe_perm_test_result <- oneway_test(combined_nc_nwe ~ nc_nwe_group_labels, distribution = approximate(nresample = 10000))
nc_nwe_effect_size <- cohen.d(combined_nc_nwe ~ nc_nwe_group_labels)

## SW vs NWE
combined_sw_nwe <- c(sw_tbl_pairwise_dist, nwe_tbl_pairwise_dist)
sw_nwe_group_labels <- as.factor(rep(c("SW", "NWE"), c(length(sw_tbl_pairwise_dist), length(nwe_tbl_pairwise_dist))))
sw_nwe_perm_test_result <- oneway_test(combined_sw_nwe ~ sw_nwe_group_labels, distribution = approximate(nresample = 10000))
nc_nwe_effect_size <- cohen.d(combined_sw_nwe ~ sw_nwe_group_labels)



# 2. PAIRWISE SNP DISTANCE
# load SNPs
#library(poppr)
nc_merged_vcfs <- read.vcfR("./VCFs/nc_merged_filtered_vcfs_output.vcf.gz")
sw_merged_vcfs <- read.vcfR("./VCFs/sw_merged_filtered_vcfs_output.vcf")
nwe_merged_vcfs <- read.vcfR("./VCFs/nwe_merged_filtered_vcfs_output.vcf")
nc_epi_linked_vcfs <- read.vcfR("./VCFs/nc_epi_linked.vcf")




# Obtain genelight object
nc_genlight <- vcfR2genlight(nc_merged_vcfs)  
sw_genlight <- vcfR2genlight(sw_merged_vcfs)
nwc_genlight <- vcfR2genlight(nwe_merged_vcfs)
nc_epi_linked_genlight <- vcfR2genlight(nc_epi_linked_vcfs)



# Using genelight object, get SNP pairwise distance matrix
# NC
nc_matrix <- as.matrix(nc_genlight)
nc_pairwise <- dist.gene(nc_matrix, method = "pairwise")
nc_pairwise_dist_matrix <- as.matrix(nc_pairwise)

#NC Epilinked
nc_epilinked_matrix <- as.matrix(nc_epi_linked_genlight) 
nc_epilinked_pairwise <- dist.gene(nc_epilinked_matrix, method = "pairwise")
nc_epilinked_dist_matrix <- as.matrix(nc_epilinked_pairwise)

# SW
sw_matrix <- as.matrix(sw_genlight)
sw_pairwise <- dist.gene(sw_matrix, method = "pairwise")
sw_pairwise_dist_matrix <- as.matrix(sw_pairwise)

# NWE
nwe_matrix <- as.matrix(nwc_genlight)
nwe_pairwise <- dist.gene(nwe_matrix, method = "pairwise")
nwe_pairwise_dist_matrix <- as.matrix(nwe_pairwise)

## Descriptive statistics on Pairwise SNP matrix
nc_pairwise_dist_distribution <- nc_pairwise_dist_matrix[upper.tri(nc_pairwise_dist_matrix)]
nc_pairwise_dist_distribution <- nc_pairwise_dist_distribution[nc_pairwise_dist_distribution > 0]

sw_pairwise_dist_distribution <- sw_pairwise_dist_matrix[upper.tri(sw_pairwise_dist_matrix)]
sw_pairwise_dist_distribution <- sw_pairwise_dist_distribution[sw_pairwise_dist_distribution >0]

nwe_pairwise_dist_distribution <- nwe_pairwise_dist_matrix[upper.tri(nwe_pairwise_dist_matrix)]
nwe_pairwise_dist_distribution <- nwe_pairwise_dist_distribution[nwe_pairwise_dist_distribution > 0]

# mean and median
mean(nc_pairwise_dist_distribution)
median(nc_pairwise_dist_distribution)

mean(sw_pairwise_dist_distribution)
median(sw_pairwise_dist_distribution)

mean(nwe_pairwise_dist_distribution)

median(nwe_pairwise_dist_distribution)


###  Permutation test on SNP pairwise distribution
# NC vs SW
combined_nc_sw_snps <- c(nc_pairwise_dist_distribution, sw_pairwise_dist_distribution)
nc_sw_snp_group_labels <- as.factor(rep(c("NC", "SW"), c(length(nc_pairwise_dist_distribution), length(sw_pairwise_dist_distribution))))
nc_sw_snp_perm_test_result <- oneway_test(combined_nc_sw_snps ~ nc_sw_snp_group_labels, distribution = approximate(nresample = 10000))
nc_sw_snp_effect_size <- cohen.d(combined_nc_sw_snps ~ nc_sw_snp_group_labels)

# NC vs NWE
combined_nc_nwe_snps <- c(nc_pairwise_dist_distribution, nwe_pairwise_dist_distribution)
nc_nwe_snp_group_labels <- as.factor(rep(c("NC", "NWE"), c(length(nc_pairwise_dist_distribution), length(nwe_pairwise_dist_distribution))))
nc_nwe_snp_perm_test_result <- oneway_test(combined_nc_nwe_snps  ~ nc_nwe_snp_group_labels, distribution = approximate(nresample = 10000))
nc_nwe_snp_effect_size <- cohen.d(combined_nc_nwe_snps  ~ nc_nwe_snp_group_labels)

# SW vs NWE
combined_sw_nwe_snps <- c(sw_pairwise_dist_distribution, nwe_pairwise_dist_distribution)
sw_nwe_snp_group_labels <- as.factor(rep(c("SW", "NWE"), c(length(sw_pairwise_dist_distribution), length(nwe_pairwise_dist_distribution))))
sw_nwe_snp_perm_test_result <- oneway_test(combined_sw_nwe_snps  ~ sw_nwe_snp_group_labels, distribution = approximate(nresample = 10000))
sw_nwe_snp_effect_size <- cohen.d(combined_sw_nwe_snps  ~ sw_nwe_snp_group_labels)


### Clustering ANalysis

# Calculate SNP Clustering Rate using a defined threshold
threshold_distance <- 15
pairs_below_threshold <- sum(nc_pairwise_dist_distribution < threshold_distance)
total_pairs <- length(nc_pairwise_dist_distribution)
clustering_rate <- pairs_below_threshold / total_pairs


# Perform clustering based on SNP distances and a SNP threshold of 10 (h=10)
clusters <- agnes(nc_pairwise_dist_matrix, diss = TRUE, method = "average")
clusters <- as.data.frame(cutree(as.hclust(clusters), h = 30))
colnames(clusters) <- "cluster_id"

# Get cluster_ids that are "duplicated" (meaning they have more than one patient)
clusterIDs <- unique(subset(clusters, duplicated(clusters$cluster_id))$cluster_id)
# Get samples within these "duplicated" clusters
clusters <- subset(clusters, cluster_id %in% clusterIDs)
# Add proper Sample name column
clusters <- cbind(Sample = rownames(clusters), clusters)
print(clusters)


#### Distribution of Nucleotide diversity
nc_nuc_div_dist <- read_tsv("./Results/nc_nuc_div_pi_only.txt")
nc_nuc_div_dist <- nc_nuc_div_dist$PI
sw_nuc_div_dist <- read_tsv("./Results/sw_nuc_div_pi_only.txt")
sw_nuc_div_dist <- sw_nuc_div_dist$PI
nwe_nuc_div_dist <- read_tsv("./Results/nwe_nuc_div_pi_only.txt")
nwe_nuc_div_dist <- nwe_nuc_div_dist$PI


## Permutation test to Compare distributions
# NC vs SW
combined_nc_sw_nuc_div <- c(nc_nuc_div_dist, sw_nuc_div_dist)
nc_sw_nuc_group_labels <- as.factor(rep(c("NC", "SW"), c(length(nc_nuc_div_dist), length(sw_nuc_div_dist))))
nc_sw_nuc_perm_test_result <- oneway_test(combined_nc_sw_nuc_div ~ nc_sw_nuc_group_labels, distribution = approximate(nresample = 10000)) # Sig, p-value < 1e-04
nc_sw_nuc_effect_size <- cohen.d(combined_nc_sw_nuc_div ~ nc_sw_nuc_group_labels) # d estimate: -0.320886 (small)

# NC vs NWE
combined_nc_nwe_nuc_div <- c(nc_nuc_div_dist, nwe_nuc_div_dist)
nc_nwe_nuc_group_labels <- as.factor(rep(c("NC", "NWE"), c(length(nc_nuc_div_dist), length(nwe_nuc_div_dist))))
nc_nwe_nuc_perm_test_result <- oneway_test(combined_nc_nwe_nuc_div ~ nc_nwe_nuc_group_labels, distribution = approximate(nresample = 10000)) # Sig, p-value < 1e-04
nc_nwe_nuc_effect_size <- cohen.d(combined_nc_nwe_nuc_div ~ nc_nwe_nuc_group_labels) # d estimate: -0.867706 (large)

# SW vs NWE
combined_sw_nwe_nuc_div <- c(sw_nuc_div_dist, nwe_nuc_div_dist)
sw_nwe_nuc_group_labels <- as.factor(rep(c("SW", "NWE"), c(length(sw_nuc_div_dist), length(nwe_nuc_div_dist))))
sw_nwe_nuc_perm_test_result <- oneway_test(combined_sw_nwe_nuc_div ~ sw_nwe_nuc_group_labels, distribution = approximate(nresample = 10000)) # Sig, p-value < 1e-04
sw_nwe_nuc_effect_size <- cohen.d(combined_sw_nwe_nuc_div ~ sw_nwe_nuc_group_labels) # d estimate: 0.5326853 (medium)




## Visualize cluster (NC EPI-Linked) in a minimum spanning tree

# Minimum spanning tree
nc_epilinked_graph_object <- graph_from_adjacency_matrix(as.matrix(nc_epilinked_dist_matrix), mode = "undirected", weighted = TRUE)
epi_linked_mst_tree <- mst(nc_epilinked_graph_object)
epi_edges_df <- as.data.frame(as_edgelist(epi_linked_mst_tree))
colnames(epi_edges_df) <- c("from", "to")

## Plot the Minimum spanning tree
# EPI DRAW

# raw draw
visNetwork::visNetwork(nodes = data.frame(id = V(epi_linked_mst_tree)$name),
                       edges = epi_edges_df) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)

# Annotated draw

epi_linked_mst_shot <- visNetwork::visNetwork(nodes = data.frame(id = V(epi_linked_mst_tree)$name, label = V(epi_linked_mst_tree)$name),  # Add node labels
                                              edges = epi_edges_df,  # Use the edges data frame
                                              #main = "Minimum Spanning Tree",  # Title
                                              width = "100%", height = "800px") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visNodes(color = list(border = "red", background = "#1e3664"),  # Node color
           shadow = list(enabled = TRUE),
           font = list(size = 20, face = "Arial", color = "black", bold = "true")) %>%
  visEdges(arrows = "to", shadow = list(enabled = TRUE)) 




### Other analysis
summary(analysis_meta_info)

ass_table <- table(analysis_meta_info$treatment_status, analysis_meta_info$resistance_status)
table(analysis_meta_info$hiv_status, analysis_meta_info$resistance_status)
gender_cont_table <- table(analysis_meta_info$sex, analysis_meta_info$gRIF)
analysis_meta_info$hiv_status
hiv_cont_table <- chisq.test(ass_table)
chisq.test(gender_cont_table)
fisher.test(gender_cont_table)


log_reg_data <- analysis_meta_info %>%
  dplyr::mutate(treatment_status = factor(treatment_status, levels = c("NEW", "F-UP")), 
                hiv_status = factor(hiv_status, levels = c("P", "N")),
                gender = factor(sex, levels = c("M", "F")),
                resistance_status = factor(resistance_status, levels = c("R", "S")))


## Logistic regression
log_reg_model <- glm(resistance_status ~ treatment_status + gender + hiv_status, data =log_reg_data, family = "binomial")

summary(log_reg_model)
