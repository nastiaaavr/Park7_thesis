
# Astrocytes 
# Load necessary libraries
library(ggplot2)
library(forcats) # For factor manipulation
library(readr) # For reading TSV files


setwd('Desktop/th')
# Reading the Astro datasets
data1_F_Astro <- read_tsv("2024-03-14_unskrunken_siDJ1_F_Astro_vs_siNEG_F_Astro.tsv")
data1_M_Astro <- read_tsv("2024-03-14_unskrunken_siDJ1_M_Astro_vs_siNEG_M_Astro.tsv")
data2_F_Astro <- read_tsv("2024-03-14_unskrunken_siNRF2_F_Astro_vs_siNEG_F_Astro.tsv")
data2_M_Astro <- read_tsv("2024-03-14_unskrunken_siNRF2_M_Astro_vs_siNEG_M_Astro.tsv")


# Identifying differentially expressed genes with the same criteria
diff_genes_data1_F_Astro <- data1_F_Astro[data1_F_Astro$padj < 0.05 & !is.na(data1_F_Astro$padj), ]
diff_genes_data1_M_Astro <- data1_M_Astro[data1_M_Astro$padj < 0.05 & !is.na(data1_M_Astro$padj), ]
diff_genes_data2_F_Astro <- data2_F_Astro[data2_F_Astro$padj < 0.05 & !is.na(data2_F_Astro$padj), ]
diff_genes_data2_M_Astro <- data2_M_Astro[data2_M_Astro$padj < 0.05 & !is.na(data2_M_Astro$padj), ]

# Creating a data frame for plotting
degs_counts_Astro <- data.frame(
  Dataset = c("data1_F", "data1_M", "data2_F", "data2_M"),
  DEGs = c(nrow(diff_genes_data1_F_Astro), nrow(diff_genes_data1_M_Astro), nrow(diff_genes_data2_F_Astro), nrow(diff_genes_data2_M_Astro)),
  CustomLabel = rep(c("siPark7", "siNfe2l2"), each = 2),
  Gender = rep(c("Female", "Male"), times = 2),
  Position = c(1, 1.5, 2.5, 3)
)

# Plotting for "Astro"
ggplot(degs_counts_Astro, aes(x = Position, y = DEGs, fill = Gender)) +
  geom_bar(stat = "identity", color = "black", width = 0.4) +
  geom_text(aes(label = DEGs), vjust = -0.3, size = 3.5, color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_blank()
  ) +
  labs(
    title = "DEGs in Astrocytes",
    y = "Number of DEGs"
  ) +
  scale_fill_manual(values = c("Male" = "gray", "Female" = "pink")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_continuous(breaks = c(1.25, 2.75), labels = c("siPark7", "siNfe2l2"))




#Oligodendrocytes 
# Load necessary libraries
library(ggplot2)
library(forcats) # For factor manipulation
library(readr) # For read_tsv

# Reading the "Oligo" datasets
data1_F_Oligo <- read_tsv("2024-03-14_unskrunken_siDJ1_F_Oligo_vs_siNEG_F_Oligo.tsv")
data1_M_Oligo <- read_tsv("2024-03-14_unskrunken_siDJ1_M_Oligo_vs_siNEG_M_Oligo.tsv")
data2_F_Oligo <- read_tsv("2024-03-14_unskrunken_siNRF2_F_Oligo_vs_siNEG_F_Oligo.tsv")
data2_M_Oligo <- read_tsv("2024-03-14_unskrunken_siNRF2_M_Oligo_vs_siNEG_M_Oligo.tsv")

# Identifying differentially expressed genes with the same criteria
diff_genes_data1_F_Oligo <- data1_F_Oligo[data1_F_Oligo$padj < 0.05 & !is.na(data1_F_Oligo$padj), ]
diff_genes_data1_M_Oligo <- data1_M_Oligo[data1_M_Oligo$padj < 0.05 & !is.na(data1_M_Oligo$padj), ]
diff_genes_data2_F_Oligo <- data2_F_Oligo[data2_F_Oligo$padj < 0.05 & !is.na(data2_F_Oligo$padj), ]
diff_genes_data2_M_Oligo <- data2_M_Oligo[data2_M_Oligo$padj < 0.05 & !is.na(data2_M_Oligo$padj), ]

# Creating a data frame for plotting
degs_counts_Oligo <- data.frame(
  Dataset = c("data1_F", "data1_M", "data2_F", "data2_M"),
  DEGs = c(nrow(diff_genes_data1_F_Oligo), nrow(diff_genes_data1_M_Oligo), nrow(diff_genes_data2_F_Oligo), nrow(diff_genes_data2_M_Oligo)),
  CustomLabel = rep(c("siPark7", "siNfe2l2"), each = 2),
  Gender = rep(c("Female", "Male"), times = 2),
  Position = c(1, 1.5, 2.5, 3)
)

# Plotting for "Oligo"
ggplot(degs_counts_Oligo, aes(x = Position, y = DEGs, fill = Gender)) +
  geom_bar(stat = "identity", color = "black", width = 0.4) +
  geom_text(aes(label = DEGs), vjust = -0.3, size = 3.5, color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_blank()
  ) +
  labs(
    title = "DEGs in Oligodendrocytes",
    y = "Number of DEGs"
  ) +
  scale_fill_manual(values = c("Male" = "gray", "Female" = "pink")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_continuous(breaks = c(1.25, 2.75), labels = c("siPark7", "siNfe2l2"))



#Microglia 
# Load necessary libraries
library(ggplot2)
library(forcats) # For factor manipulation
library(readr) # For read_tsv

# Reading the "Microglia" datasets
data1_F_Microglia <- read_tsv("2024-03-14_unskrunken_siDJ1_F_Microglia_vs_siNEG_F_Microglia.tsv")
data1_M_Microglia <- read_tsv("2024-03-14_unskrunken_siDJ1_M_Microglia_vs_siNEG_M_Microglia.tsv")
data2_F_Microglia <- read_tsv("2024-03-14_unskrunken_siNRF2_F_Microglia_vs_siNEG_F_Microglia.tsv")
data2_M_Microglia <- read_tsv("2024-03-14_unskrunken_siNRF2_M_Microglia_vs_siNEG_M_Microglia.tsv")

# Identifying differentially expressed genes with the same criteria
diff_genes_data1_F_Microglia <- data1_F_Microglia[data1_F_Microglia$padj < 0.05 & !is.na(data1_F_Microglia$padj), ]
diff_genes_data1_M_Microglia <- data1_M_Microglia[data1_M_Microglia$padj < 0.05 & !is.na(data1_M_Microglia$padj), ]
diff_genes_data2_F_Microglia <- data2_F_Microglia[data2_F_Microglia$padj < 0.05 & !is.na(data2_F_Microglia$padj), ]
diff_genes_data2_M_Microglia <- data2_M_Microglia[data2_M_Microglia$padj < 0.05 & !is.na(data2_M_Microglia$padj), ]

# Creating a data frame for plotting
degs_counts_Microglia <- data.frame(
  Dataset = c("data1_F", "data1_M", "data2_F", "data2_M"),
  DEGs = c(nrow(diff_genes_data1_F_Microglia), nrow(diff_genes_data1_M_Microglia), nrow(diff_genes_data2_F_Microglia), nrow(diff_genes_data2_M_Microglia)),
  CustomLabel = rep(c("siPark7", "siNfe2l2"), each = 2),
  Gender = rep(c("Female", "Male"), times = 2),
  Position = c(1, 1.5, 2.5, 3)
)

# Plotting for "Microglia"
ggplot(degs_counts_Microglia, aes(x = Position, y = DEGs, fill = Gender)) +
  geom_bar(stat = "identity", color = "black", width = 0.4) +
  geom_text(aes(label = DEGs), vjust = -0.3, size = 3.5, color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_blank()
  ) +
  labs(
    title = "DEGs in Microglia",
    y = "Number of DEGs"
  ) +
  scale_fill_manual(values = c("Male" = "gray", "Female" = "pink")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_continuous(breaks = c(1.25, 2.75), labels = c("siPark7", "siNfe2l2"))


#Volcno plots for Microglia 
library(ggplot2)
library(readr)
library(dplyr)

# Load the datasets
data1_F_Microglia <- read_csv("diff_genes_data1_F_Microglia.csv")
data1_M_Microglia <- read_csv("diff_genes_data1_M_Microglia.csv")


create_volcano_plot <- function(data, title) {
  # Filter to include only significant upregulated and downregulated genes for labeling
  labeled_genes <- data %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    select(gene_name, log2FoldChange, padj) # Ensure these are your correct column names
  
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = factor(ifelse(padj < 0.05 & abs(log2FoldChange) > 1, 
                                         ifelse(log2FoldChange > 1, "Upregulated", "Downregulated"), 
                                         "Not significant"))), size = 1.5) +
    geom_text(data = labeled_genes, aes(label = gene_name), vjust = 1.5, hjust = 0.5, check_overlap = TRUE, size = 2.5) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey"),
                       name = "", # You can give this a name or leave it blank for no title
                       labels = c("Not significant", "Downregulated", "Upregulated")) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10(Adjusted P-Value)", color = "Gene Expression\nCategory") +
    theme_minimal() +
    theme(legend.position = "right", 
          legend.title.align = 0.5,
          legend.background = element_blank(), 
          legend.box.background = element_rect(color="black", size=0.5),
          legend.text = element_text(size=10))
}


volcano_plot_F_Microglia <- create_volcano_plot(data1_F_Microglia, "Volcano Plot - Female Microglia")
volcano_plot_M_Microglia <- create_volcano_plot(data1_M_Microglia, "Volcano Plot - Male Microglia")

# Display the plots F for females and M for males 
print(volcano_plot_F_Microglia)
print(volcano_plot_M_Microglia)

#Directinality of qualitative chnage for Nfe2l2 knockdown 

library(ggplot2)
library(tidyverse)
library(ggrepel)

full_join(
  read_tsv("2024-03-14_unskrunken_siNRF2_F_Microglia_vs_siNEG_F_Microglia.tsv",
           show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  read_tsv("2024-03-14_unskrunken_siNRF2_M_Microglia_vs_siNEG_M_Microglia.tsv",
           show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  by = "gene_id", suffix = c("_F", "_M")
) |>  
  filter(padj_F < 0.05 | padj_M < 0.05 | gene_name_F == "Park7") |>
  ggplot(aes(x = log2FoldChange_F, y = log2FoldChange_M, label = gene_name_F)) +
  geom_point(data = \(x) filter(x, abs(log2FoldChange_F - log2FoldChange_M) > 0.4)) +
  geom_point(data = \(x) filter(x, !abs(log2FoldChange_F - log2FoldChange_M) > 0.4),
             alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_text_repel(data = \(x) filter(x, abs(log2FoldChange_F - log2FoldChange_M) > 0.4),
                  size = 3, max.overlaps = 40) +
  geom_point(data = \(x) filter(x, gene_name_F == "Park7"), color = "steelblue", size = 2) +
  geom_text_repel(data = \(x) filter(x, gene_name_F == "Nfe2l2"), color = "steelblue", size = 4) +
  labs(x = "Log2 fold change (female)", y = "Log2 fold change (male)",
       color = "Gene", title = "Significant genes M/F siNfe2l2 Microglia",
       caption = "Gene names for absolute diff of logFC > 0.4") +
  theme_minimal() 


filtered_genes_male <- full_join(
  read_tsv("2024-03-14_unskrunken_siNRF2_F_Microglia_vs_siNEG_F_Microglia.tsv", show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  read_tsv("2024-03-14_unskrunken_siNRF2_M_Microglia_vs_siNEG_M_Microglia.tsv", show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  by = "gene_id", suffix = c("_F", "_M")
) |>  
  filter((padj_F < 0.05 | padj_M < 0.05) & abs(log2FoldChange_F - log2FoldChange_M) > 0.4)

filtered_genes_male$gene_name_F


library(ggplot2)
library(tidyverse)
library(ggrepel)

full_join(
  read_tsv("2024-03-14_unskrunken_siNRF2_F_Microglia_vs_siNEG_F_Microglia.tsv",
           show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  read_tsv("2024-03-14_unskrunken_siNRF2_M_Microglia_vs_siNEG_M_Microglia.tsv",
           show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  by = "gene_id", suffix = c("_F", "_M")
) |>  
  filter(padj_F < 0.05 | padj_M < 0.05 | gene_name_F == "Park7") |>
  ggplot(aes(x = log2FoldChange_F, y = log2FoldChange_M, label = gene_name_F)) +
  geom_point(data = \(x) filter(x, abs(log2FoldChange_F - log2FoldChange_M) > 0.4)) +
  geom_point(data = \(x) filter(x, !abs(log2FoldChange_F - log2FoldChange_M) > 0.4),
             alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_text_repel(data = \(x) filter(x, abs(log2FoldChange_F - log2FoldChange_M) > 0.4),
                  size = 3, max.overlaps = 40) +
  geom_point(data = \(x) filter(x, gene_name_F == "Park7"), color = "steelblue", size = 2) +
  geom_text_repel(data = \(x) filter(x, gene_name_F == "Nfe2l2"), color = "steelblue", size = 4) +
  labs(x = "Log2 fold change (female)", y = "Log2 fold change (male)",
       color = "Gene", title = "Significant genes M/F siNfe2l2 Microglia",
       caption = "Gene names for absolute diff of logFC > 0.4") +
  theme_minimal() 


filtered_genes_male <- full_join(
  read_tsv("2024-03-14_unskrunken_siNRF2_F_Microglia_vs_siNEG_F_Microglia.tsv", show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  read_tsv("2024-03-14_unskrunken_siNRF2_M_Microglia_vs_siNEG_M_Microglia.tsv", show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  by = "gene_id", suffix = c("_F", "_M")
) |>  
  filter((padj_F < 0.05 | padj_M < 0.05) & abs(log2FoldChange_F - log2FoldChange_M) > 0.4)

filtered_genes_male$gene_name_F


#Directinality of qualitative chnage for Park7 knockdown 


full_join(
  read_tsv("2024-03-14_unskrunken_siDJ1_F_Microglia_vs_siNEG_F_Microglia.tsv",
           show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  read_tsv("2024-03-14_unskrunken_siDJ1_M_Microglia_vs_siNEG_M_Microglia.tsv",
           show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  by = "gene_id", suffix = c("_F", "_M")
) |>  
  filter(padj_F < 0.05 | padj_M < 0.05 | gene_name_F == "Park7") |>
  ggplot(aes(x = log2FoldChange_F, y = log2FoldChange_M, label = gene_name_F)) +
  geom_point(data = \(x) filter(x, abs(log2FoldChange_F - log2FoldChange_M) > 0.4)) +
  geom_point(data = \(x) filter(x, !abs(log2FoldChange_F - log2FoldChange_M) > 0.4),
             alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_text_repel(data = \(x) filter(x, abs(log2FoldChange_F - log2FoldChange_M) > 0.4),
                  size = 3, max.overlaps = 40) +
  geom_point(data = \(x) filter(x, gene_name_F == "Park7"), color = "steelblue", size = 2) +
  geom_text_repel(data = \(x) filter(x, gene_name_F == "Nfe2l2"), color = "steelblue", size = 4) +
  labs(x = "Log2 fold change (female)", y = "Log2 fold change (male)",
       color = "Gene", title = "Significant genes M/F siNfe2l2 Microglia",
       caption = "Gene names for absolute diff of logFC > 0.4") +
  theme_minimal() 


filtered_genes_male <- full_join(
  read_tsv("2024-03-14_unskrunken_siNRF2_F_Microglia_vs_siNEG_F_Microglia.tsv", show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  read_tsv("2024-03-14_unskrunken_siNRF2_M_Microglia_vs_siNEG_M_Microglia.tsv", show_col_types = FALSE) |> 
    select(log2FoldChange, padj, starts_with("gene")),
  by = "gene_id", suffix = c("_F", "_M")
) |>  
  filter((padj_F < 0.05 | padj_M < 0.05) & abs(log2FoldChange_F - log2FoldChange_M) > 0.4)

filtered_genes_male$gene_name_F


# analysis of counts 

library(dplyr)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)



#| eval: false

if (!is.data.frame(data)) {
  if (is.matrix(data)) {
    gene_ids <- rownames(data)
    data <- as.data.frame(data)
    data$gene_id <- gene_ids
    # Move gene_id column to the first position
    data <- data %>% select(gene_id, everything())
  } else {
    stop("Data is not in an expected format (data frame or matrix).")
  }
}


colnames(data) <- gsub("SM", "AM", colnames(data))
colnames(data) <- gsub("siNeg", "siNEG", colnames(data), ignore.case = TRUE)
colnames(data) <- gsub("siDj1", "siDJ1", colnames(data), ignore.case = TRUE)
colnames(data) <- gsub("siNrf2", "siNRF2", colnames(data), ignore.case = TRUE)


if ("gene_id" %in% colnames(data)) {
  # Revert any unintended replacement in gene_id column name itself
  colnames(data)[colnames(data) == "geneAM_id"] <- "gene_id"
}


saveRDS(data, "12024-02-26_counts_modified.rds")

data <- readRDS("12024-02-26_counts_modified.rds")

#Enter ENSEMBL_id of gene of interest 

filtered_data <- data %>%
  filter(gene_id %in% c("ENSMUSG00000015839", ""))
filtered_data


filtered_data %>%
  pivot_longer(cols = -gene_id, 
               names_to = "sample_name", 
               values_to = "count") |>
  separate(col = sample_name, into = c("condition", "cell_type_sex", "replicate"), sep = "_") %>%
  separate(col = cell_type_sex, into = c("cell_type", "sex"), sep = "(?<=^[AMO])(?=[MF])", extra = "merge") %>%
  mutate(replicate = as.numeric(replicate)) |>
  group_by(gene_id, condition, cell_type, sex) %>%
  summarise(mean_expression = mean(count), 
            sd_expression = sd(count), 
            n = n(),
            .groups = 'drop') -> summary_data
summary_data

#comparison of gene counts siNeg vs siPark7 for selected gene 
summary_data %>%
  left_join(gene_correspondences) %>%
  filter(condition != "siNRF2") %>% 
  mutate(condition = fct_relevel(condition, "siNEG")) %>% 
  ggplot(aes(x = interaction(sex, cell_type), y = mean_expression, fill = condition)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax =  mean_expression + sd_expression),
                position = position_dodge(width = 0.95),
                width = 0.2) +
  facet_wrap(vars(gene_name), scales = "free_y") +
  labs(title = "Gene Expression for siPark7 Condition",
       x = "Cell Type and Condition", y = "Mean Expression",
       caption = "Error bars = standard deviation to the mean") +
  theme_minimal()  +
  theme(strip.text = element_text(face = "bold", size = 13)) +
  scale_fill_manual(values = c("siNEG" = "gray", "siDJ1" = "steelblue", "siNRF2" = "pink"))


#Plotting gene expression,normalized to siNeg with baseline 1

library(dplyr)
library(ggplot2)
library(forcats)

# Join summary_data with gene_correspondences, normalize to siNEG, and plot only siDJ1/siNEG
summary_data %>%
  left_join(gene_correspondences, by = "gene_id") %>%
  filter(condition != "siNRF2") %>%
  # Normalize mean_expression to siNEG
  group_by(gene_id, sex, cell_type) %>%
  mutate(mean_expression_siNEG = mean_expression[condition == "siNEG"]) %>%
  ungroup() %>%
  mutate(
    normalized_expression = mean_expression / mean_expression_siNEG,
    sd_normalized = sd_expression / mean_expression_siNEG
  ) %>%
  filter(condition == "siDJ1") %>%  # Keep only the siDJ1 condition
  ggplot(aes(x = interaction(sex, cell_type), y = normalized_expression, fill = condition)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = normalized_expression - sd_normalized, 
                    ymax = normalized_expression + sd_normalized),
                position = position_dodge(width = 0.95),
                width = 0.2) +
  facet_wrap(vars(gene_name), scales = "free_y") +
  labs(title = "Normalized siDJ1 Expression Relative to siNEG",
       x = "Cell Type and Condition", 
       y = "Normalized Expression (siDJ1 / siNEG)",
       caption = "Error bars = standard deviation to the mean") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 13)) +
  scale_fill_manual(values = c("siDJ1" = "steelblue")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")  # Add baseline at 1 for siNEG



#New volcano for the thesis

library(ggplot2)
library(readr)
library(dplyr)

data1_F_Microglia <- read_csv("diff_genes_data1_F_Microglia.csv")
data1_M_Microglia <- read_csv("diff_genes_data1_M_Microglia.csv")


create_volcano_plot <- function(data, title) {
  labeled_genes <- data %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    select(gene_name, log2FoldChange, padj) # Ensure these are your correct column names
  
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = factor(ifelse(padj < 0.05 & abs(log2FoldChange) > 1, 
                                         ifelse(log2FoldChange > 1, "Upregulated", "Downregulated"), 
                                         "Not significant"))), size = 1.5) +
    geom_text(data = labeled_genes, aes(label = gene_name), vjust = 1.5, hjust = 0.5, check_overlap = TRUE, size = 3) +
    geom_hline(yintercept = -log10(0.05), color = "pink", linetype = "dashed") +
    geom_vline(xintercept = 1, color = "pink", linetype = "dashed") +
    geom_vline(xintercept = -1, color = "pink", linetype = "dashed") +
    scale_color_manual(values = c("Upregulated" = "pink", "Downregulated" = "steelblue", "Not significant" = "grey"),
                       name = "Gene expression", 
                       labels = c("Downregulated","Not significant",  "Upregulated")) +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(padj)", color = "Gene Expression\nCategory", caption = "Infertior to 5% threshold") +
    theme_minimal() +
    theme(legend.position = "right", 
          legend.title.align = 0.5,
          legend.background = element_blank(), 
          legend.box.background = element_rect(color="black", size=0.5),
          legend.text = element_text(size=10),
          plot.caption = element_text(hjust = 1), # Right-align the caption
          plot.caption.position = "plot") + # Position the caption at the bottom of the plot
    ylim(c(0, NA))
}

volcano_plot_F_Microglia <- create_volcano_plot(data1_F_Microglia, "Volcano Plot - Female siPark7 Microglia")
volcano_plot_M_Microglia <- create_volcano_plot(data1_M_Microglia, "Volcano Plot - Male siPark7 Microglia")


print(volcano_plot_F_Microglia)
print(volcano_plot_M_Microglia)



#Also, violin plots that show sample dispersion

dds_relabel <- read_rds("2024-02-26_dds_label_swapped.rds")
colData(dds_relabel) |>
  as_tibble() |>
  mutate(keep =  sex == "M" & cell_type == "Astro" & condition %in% c("siNEG", "siDJ1")) |> 
  pull(keep) -> keep

# subset_dds function from https://github.com/sysbiolux/Helgueta_et_al_2023/blob/b24ad6d05108643b270327a00bf464852692c2dd/RNA-seq_Sex-differences_8months_mice/Brain_sex_dispersion.qmd#L134

relabel_F_Microglia_DJ1 <- subset_dds(dds_relabel, sex == "F" & cell_type == "Microglia" & condition %in% c("siNEG", "siDJ1"))
relabel_M_Microglia_DJ1 <- subset_dds(dds_relabel, sex == "M" & cell_type == "Microglia"  & condition %in% c("siNEG", "siDJ1"))
relabel_F_Microglia_NRF2 <- subset_dds(dds_relabel, sex == "F" & cell_type == "Microglia" & condition %in% c("siNEG", "siNRF2"))
relabel_M_Microglia_NRF2 <- subset_dds(dds_relabel, sex == "M" & cell_type == "Microglia"  & condition %in% c("siNEG", "siNRF2"))
relabel_F_Microglia_NEG <- subset_dds(dds_relabel, sex == "F" & cell_type == "Microglia" & condition == "siNEG")
relabel_M_Microglia_NEG <- subset_dds(dds_relabel, sex == "M" & cell_type == "Microglia"  & condition == "siNEG")

bind_rows(
  M_Microglia_DJ1 = relabel_M_Microglia_DJ1,
  F_Microglia_DJ1 = relabel_F_Microglia_DJ1,
  M_Microglia_NRF2 = relabel_M_Microglia_NRF2,
  F_Microglia_NRF2 = relabel_F_Microglia_NRF2,
  M_Microglia_NEG = relabel_M_Microglia_NEG,
  F_Microglia_NEG = relabel_F_Microglia_NEG,
  .id = "id"
) |> 
  mutate(sex = if_else(str_detect(id, "^F"), "F", "M")) |>
  ggplot(aes(x = id, y = dispersion, colour = sex)) +
  ggbeeswarm::geom_quasirandom(color = "grey90", alpha = 0.1) +
  geom_boxplot(alpha = 0.6) +
  scale_y_log10() +
  theme_bw(base_size = 10) + 
  labs(x = NULL)
