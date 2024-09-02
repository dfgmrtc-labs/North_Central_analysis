library(tidyverse)
library(readxl)
library(plotly)
library(vioplot)


source("./statistical_analysis.R")

###################  ANALYIS REPOS ########################
setwd("/home/dfgmrtc/NC_ANALYSIS/North_Central_analysis_2023")
### Load data
## Detailed resistance info
detailed_resistance_mutation_info <- read_tsv("./Results/detailed_drug_resistance_profile.tsv")
detailed_resistance_mutation_info <- detailed_resistance_mutation_info %>%
  dplyr::filter(conf_grading != "combo")

## extract drug_resistance variants only
detailed_resistance_only <- detailed_resistance_mutation_info %>%
  dplyr::filter(conf_grading == "Assoc_w_R" | conf_grading == "Assoc_w_R_Interim")


### load drug resistance variant frequency

dr_variant_freq <- read_xlsx("./Results/Resistant_Mutation_and_Genes.xlsx")

# Load Short drug_resistance info
short_resistance_var <- read_tsv("./Results/short_drug_resistance_profile.tsv")

# Load curated meta data
curated_meta_data <- read_tsv("./Metadata/metadata_2.tsv")

# Preliminary statistics
table(curated_meta_data$lineage, curated_meta_data$resistance_status)
chisq.test(curated_meta_data$lineage, curated_meta_data$resistance_status)
fisher.test(curated_meta_data$lineage, curated_meta_data$resistance_status)
no_NA_hiv <- curated_meta_data %>%
  filter(hiv_status != "NA")
fisher.test(no_NA_hiv$hiv_status, no_NA_hiv$gRIF)
length(as.vector(no_NA_hiv$hiv_status))
table(curated_meta_data$STATE, curated_meta_data$resistance_status)
state_rs <- fisher.test(curated_meta_data$STATE, curated_meta_data$resistance_status)

fisher.test(curated_meta_data$sex, curated_meta_data$resistance_status)
fisher.test(curated_meta_data$age, curated_meta_data$resistance_status)

state_rs$alternative


######## PLOTS ############
# 1. DOT PLOT with Plotly
# Reform short dr info for a dot plot

df_4_dot_plot <- short_resistance_var %>%
  filter(who_drug_res_profile != "Drug-Susceptible") %>%
  separate(ass_res_drugs, into = c("ETH", "INH", "PZA", "RIF", "STM"), sep = ",")

# df_4_dot_plot was saved to memory and manually curated to df_4_dot_plot_2 which is then loaded below

df_4_dot_plot_2 <- read_tsv("./Results/df_4_dot_plot_2.tsv")

df_4_dot_plot_2 <- df_4_dot_plot_2 %>%
  gather("EMB", "ETH", "INH", "PZA", "RIF", "STM", key = "drug", value = "resistance_status")


# Plot the first plot
custom_color <- c("#990000", "#fff3e7")
plot_ly(data = df_4_dot_plot_2, x = ~sample_name, y = ~drug, color = ~resistance_status, 
        type = "scatter", mode = "markers", marker = list(size = 20), colors = custom_color) %>%
  layout(
    xaxis = list(title = "Sample Name", tickangle = 45, tickfont = list(family = "Arial", weight = "bold")),
    yaxis = list(title = "Drug", tickfont = list(weight = "bold")))


# Plot the second plot
custom_color_2 <- c("#f47920","#ff0096", "#76a5af", "#990000", "#faab36")
plot_ly(data = df_4_dot_plot_2, x = ~sample_name, y = ~res_profile, color = ~res_profile, 
                       type = "scatter", mode = "markers", marker = list(size = 20), colors = custom_color_2) %>%
  layout(
    xaxis = list(title = "Sample Name", tickangle = 45, tickfont = list(family = "Arial", weight = "bold")),
    yaxis = list(title = "Resistance Profile", tickfont = list(weight = "bold")))


#### 2. Bar PLOTS
# ANti-tuberculosis drugs frequencies
# detailed_resistance_mutation_info was loaded to memory and manually curated for consistency sake to generate detailed_resistance_mutation_info_edited
write_tsv(detailed_resistance_mutation_info, file = "./Results/detailed_resistance_mutation_info.tsv")

detailed_resistance_mutation_info_edited <- read_tsv("./Results/detailed_resistance_mutation_info_edited.tsv")

custom_color_3 <- c("#f47920", "#ff0096", "#76a5af", "#e0b0ff", "#990000", "#ffd770")
res_legend_labels <- c("EMB: Ethambutol", "ETH: Ethionamide", "INH: Isoniazid", "PZA: Pyrazinamide", "RIF: Rifampicin", "SM: Streptomycin")

detailed_resistance_mutation_info_edited %>%
  filter(conf_grading == "Assoc_w_R" | conf_grading == "Assoc_w_R_Interim") %>%
  filter(variant != "c-15t" & variant != "L203L") %>%
  ggplot(aes(x = factor(associated_drug), fill = associated_drug)) + 
  geom_bar() +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.5, size = 7) +
  scale_fill_manual(values = custom_color_3,
                    labels = res_legend_labels) +
  labs(x = "Anti-tuberculosis drugs", y = "Total number of samples resistant to anti-tuberculosis drugs") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(face = "bold", size = 10),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(family = "Arial", size = 12, face = "bold"), 
        axis.title.y = element_text(family = "Arial", size = 15, face = "bold"),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        legend.text = element_text(family = "Arial", size = 10), legend.title = element_blank(), 
        plot.background = element_rect(fill = "white", colour = "white"))



## WHO resistance categories
# nc_shortresistance_var written to memory and manually curated for consistency sake
write_tsv(x = short_resistance_var, file = "./Results/nc_shortresistance_var.tsv")

short_resistance_var_edited <- read_tsv(file = "./Results/nc_shortresistance_var_edited.tsv")

custom_color_4 <- c("#f47920", "#ff0096", "#76a5af", "#e0b0ff", "#990000", "#ffd770")
who_legend_labels <- c("DS: Susceptible", "INH: Isoniazid", "MDR: Multidrug", "PZN: Pyrazinamide", "RR: Rifampicin", "SM: Streptomycin")

ggplot(short_resistance_var_edited, aes(x = factor(WHO_DR_STATUS), fill = WHO_DR_STATUS)) +
  geom_bar() +
  scale_fill_manual(values = custom_color_4,
                    labels = who_legend_labels) +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.5, size = 6) +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(face = "bold", size = 10),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(family = "Arial", size = 12, face = "bold"), 
        axis.title.y = element_text(family = "Arial", size = 15, face = "bold"),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        legend.text = element_text(family = "Arial", size = 10), legend.title = element_blank(), 
        plot.background = element_rect(fill = "white", colour = "white")) +
  xlab( "WHO RESISTANCE STATUS") +
  ylab("Frequency")



### Visualize TBL distribution using vioplot

summary(nc_tbl_pairwise_dist)
summary(sw_tbl_pairwise_dist)
summary(nwe_tbl_pairwise_dist)


boostrap_nc_tbl_pairwise_dist <- bootstrap_resampling(nc_tbl_pairwise_dist, 1000)
boostrap_sw_tbl_pairwise_dist <- bootstrap_resampling(sw_tbl_pairwise_dist, 1000)
boostrap_nwe_tbl_pairwise_dist <- bootstrap_resampling(nwe_tbl_pairwise_dist, 1000)

north_central_tblc_v <- as.vector(boostrap_nc_tbl_pairwise_dist)
south_west_tbl_v <- as.vector(boostrap_sw_tbl_pairwise_dist)
north_we_tbl_v <- as.vector(boostrap_nwe_tbl_pairwise_dist)

summary(north_central_tblc_v)
summary(south_west_tbl_v)
summary(north_we_tbl_v)


par(mfrow = c(1, 1))
medians <- reorder(boxplot_tbl_dist_df$region, boxplot_tbl_dist_df$terminal_branch_length, median)
vioplot(nc_tbl_pairwise_dist, sw_tbl_pairwise_dist, nwe_tbl_pairwise_dist, names = c("North_Central", "South_West", "North_West_East"), col = custom_color_5)


#### Distribution of Nucleotide diversity
summary(nc_nuc_div_dist)
summary(sw_nuc_div_dist)
summary(nwe_nuc_div_dist)

boostrap_nc_nuc_div_dist <- bootstrap_resampling(nc_nuc_div_dist, 1000)
boostrap_sw_nuc_div_dist <- bootstrap_resampling(sw_nuc_div_dist, 1000)
boostrap_nwe_nuc_div_dist <- bootstrap_resampling(nwe_nuc_div_dist, 1000)

length(boostrap_nc_nuc_div_dist) == length(boostrap_sw_nuc_div_dist)

north_central_nuc_v <- as.vector(boostrap_nc_nuc_div_dist)
south_west_nuc_v <- as.vector(boostrap_sw_nuc_div_dist)
north_we_nuc_v <- as.vector(boostrap_nwe_nuc_div_dist)

length(north_central_nuc_v) == length(south_west_nuc_v)
summary(north_central_nuc_v)
summary(south_west_nuc_v)
summary(north_we_nuc_v)
nuc_div_dist_df <- data.frame("North_Central" = north_central_nuc_v, "South_West" = south_west_nuc_v, "North_East_West" = north_we_nuc_v)



boxplot_nuc_dist_df <- nuc_div_dist_df %>%
  gather("North_Central", "South_West", "North_East_West", key = "region", value = "nucleotide_diversity")


# Calculate mean and median values of Nucleotide diversity
pop_means <- aggregate(nucleotide_diversity ~ region, boxplot_nuc_dist_df, mean)
pop_medians <- aggregate(nucleotide_diversity ~ region, boxplot_nuc_dist_df, median)

# Create the plot
custom_color_8 <- c("#f1bb93", "#8b0000", "#f77a20")
ggplot() +
  geom_vline(data = pop_means, aes(xintercept = region, color = region), linewidth = 1.5) +
  geom_vline(data = pop_medians, aes(xintercept = region, color = region), linewidth = 1.5) +
  geom_point(data = pop_means, aes(x = region, y = nucleotide_diversity, color = region), shape = 16, size = 6) +
  geom_point(data = pop_medians, aes(x = region, y = nucleotide_diversity, color = region), shape = 18, size = 7) +
  labs(x = "Region", y = "Nucleotide_diversity", size = 10) +
  scale_color_manual(values =  custom_color_8) +theme_bw() +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(face = "bold", size = 10),
        panel.grid.major.x = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.x = element_line(color = "gray", linewidth = 0.5),
        panel.grid.major.y = element_line(color = "gray", linewidth = 0.5),
        panel.grid.minor.y = element_line(color = "gray", linewidth = 0.5),
        #panel.grid.major.y = element_line(color = "gray", size = 0.5),
        #panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(family = "Arial", size = 12, face = "bold"), 
        axis.title.y = element_text(family = "Arial", size = 15, face = "bold"),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        legend.text = element_text(family = "Arial", size = 10, face = "bold"), legend.title = element_blank(), 
        plot.background = element_rect(fill = "white", colour = "white"))





