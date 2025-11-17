#############################
# Script to create Figure 1 #
#############################

# Author: Dave Edison Rojas Calderon

# Loading required packages
library(tidyverse)
library(ggtree)
library(tidytree)
library(ggtreeExtra)
library(ape)
library(ggplot2)
library(aplot)
library(ggstar)
library(Polychrome)

# Load input data
# ------------------------------------------------------------------------------
# Define path of input files
tree_path <- file.path("../proc_data/phylo_tree/acinetobacter_tree.treefile")
names_path <- file.path("../raw_data/gtdb-search_filtered_for phylotree.xlsx")
checkm2_path <- file.path("../proc_data/checkM2/quality_report.tsv")

# Load Input files
unrooted_tree <- read.tree(tree_path) # Unrooted tree
# Names
names <- readxl::read_excel(names_path,
  sheet = "Sheet1"
) %>%
  # change accession of isolate for "isolate"
  mutate(accession = ifelse(is.na(accession), "isolate", accession)) %>%
  # Create New names column
  mutate(
    New_Names = ifelse(accession != "isolate", paste(ncbi_organism_name, accession), ncbi_organism_name),
    New_Names = str_replace(New_Names, "GCF\\_", "(GCF_"),
    New_Names = str_replace(New_Names, "\\.1", ".1)"),
    New_Names = str_replace(New_Names, "\\.2", ".2)")
  ) %>%
  select(New_Names, accession)
# CheckM2
checkm2 <- read_delim(checkm2_path,
  delim = "\t"
) %>%
  # Extract accession
  mutate(accession = ifelse(str_detect(Name, "(GC[A|F]\\_\\d+\\.\\d{1})"), str_extract(Name, "(GC[A|F]\\_\\d+\\.\\d{1})"), "isolate")) %>%
  # Select Important columns
  select(accession, Name, Completeness_Specific, Contamination, Genome_Size, Total_Contigs)

# Define metadata
metadata <- checkm2 %>%
  left_join(names,
    by = "accession"
  ) %>%
  # Renome New_Names
  rename(
    GenomeID = New_Names,
    Completeness = Completeness_Specific,
    Accession = accession,
    Filename = Name
  ) %>%
  # Reorder
  select(GenomeID, everything())

# Root tree
rooted_tree <- root(unrooted_tree,
  outgroup = "GCF_002080125.1_ASM208012v1_genomic.fna",
  resolve.root = TRUE
)

# Cleaning Tree Names
# ------------------------------------------------------------------------------
# Create a named vector mapping Name and New_Names
name_map <- metadata %>%
  select(Filename, GenomeID) %>%
  deframe() # Creates a names vector, where names = Name, values = New_Names
# Replace the tip labels in your tree
rooted_tree$tip.label <- name_map[rooted_tree$tip.label]

# Prune tree
# ------------------------------------------------------------------------------
# Get all descendant tips of each node
desc_132 <- phytools::getDescendants(rooted_tree, node = 132)
desc_102 <- phytools::getDescendants(rooted_tree, node = 102)

# Get their tip labels
tips_to_drop <- rooted_tree$tip.label[c(desc_102, desc_132)]

# Prune the tree
pruned_tree <- drop.tip(rooted_tree, tips_to_drop)

# Make Phylogenetic tree
# ------------------------------------------------------------------------------
# Make simple phylogenetic tree
pt <- ggtree(
  pruned_tree,
  # branch.length="none",
  # layout = "fan"
) +
  coord_cartesian(clip = "off") +
  geom_tiplab(aes(colour = (label == "Acinetobacter towneri RMS-02")),
    align = TRUE,
    size = 3,
    offset = 0.01,
    hjust = 0,
    width = 0.6
  ) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey50")) +
  theme_tree2() +
  xlim(c(0, 2.6)) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 1)
  )

# Get the tip labels in the plotted order
# ------------------------------------------------------------------------------
# Get the tree data from the plot
tree_data <- pt@data
# Get tip labels ordered
tip_labels_ordered <- tree_data[tree_data$isTip, ] |>
  dplyr::arrange(y) |>
  dplyr::pull(label)

# Create manually each independent plot and then arrange them together
# ------------------------------------------------------------------------------
# Remove pruned tips
metadata <- metadata %>%
  filter(!GenomeID %in% tips_to_drop)

# Assign levels to GenomeID
metadata$GenomeID <- factor(metadata$GenomeID, levels = tip_labels_ordered)

# Make second plot: Genome Quality
quality_p <- metadata %>%
  # Pivot longer Completeness and Contoamination
  pivot_longer(
    cols = c(Completeness, Contamination),
    names_to = "Genome_Quality",
    values_to = "Percentage"
  ) %>%
  # Mutate Genome_Quality values
  mutate(Genome_Quality = ifelse(Genome_Quality == "Completeness", "Completeness (%)", "Contamination (%)")) %>%
  ggplot(aes(
    x = Percentage,
    y = GenomeID,
    fill = Genome_Quality
  )) +
  geom_bar(position = "identity", stat = "identity", width = 0.6) +
  theme_bw() +
  # Manually assign colors
  scale_fill_manual(
    values = c(
      "Completeness (%)" = "#5d7a60",
      "Contamination (%)" = "#f4c465"
    ),
    guide = guide_legend(position = "bottom"),
  ) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  # Set legen title
  labs(fill = "Genome Quality") +
  # Personalize
  # Set theme
  theme(
    axis.title.y = element_blank(), # Remove title axis y
    axis.text.y = element_blank(), # Remove text axis y
    axis.ticks.y = element_blank(), # Remove text axis y
    axis.line.y = element_blank(), # Remove text axis y
    axis.title.x = element_blank(), # Remove title axis x
    axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 1)
  ) # Rotate x axis text

# Get the legend
legend_quality <- cowplot::get_legend(quality_p)

# Remove legend from plot
quality_p <- quality_p +
  theme(legend.position = "none")

# Make fourth and last plot: Genome Size
size_p <- metadata %>%
  ggplot(aes(
    x = Genome_Size / 1e6,
    y = GenomeID,
    fill = "Genome Size (Mb)"
  )) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_bw() +
  # Personalize color
  scale_fill_manual(
    name = "Genome Size",
    values = "#bdbdbd",
    guide = guide_legend(position = "bottom")
  ) +
  # Set theme
  theme(
    axis.title.y = element_blank(), # Remove title axis y
    axis.text.y = element_blank(), # Remove text axis y
    axis.ticks.y = element_blank(), # Remove text axis y
    axis.line.y = element_blank(), # Remove text axis y
    axis.title.x = element_blank(), # Remove title axis x
    axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 1)
  ) # Rotate x axis text

# Get legend
legend_size <- cowplot::get_legend(size_p)

# Remove legend from plot
size_p <- size_p +
  theme(legend.position = "none")

# Arrange Legends
# ------------------------------------------------------------------------------
combined_legends <- cowplot::plot_grid(legend_quality, legend_size,
  align = "v",
  rel_widths = c(2, 1)
)

# Arrange Plots
# ------------------------------------------------------------------------------
joined_plots <- cowplot::plot_grid(pt, quality_p, size_p,
  nrow = 1,
  ncol = 3,
  align = "h",
  rel_widths = c(10, 1, 1)
)

# Arrange Plots and Legends
# ------------------------------------------------------------------------------
final_plot <- cowplot::plot_grid(joined_plots, combined_legends,
  nrow = 2,
  align = "h",
  rel_heights = c(5, 0.5)
)
# Save Phylogenetic tree as pdf and png
ggsave(
  plot = final_plot,
  filename = "phylogenetic_tree.pdf",
  path = "plots/"
)
ggsave(
  plot = final_plot,
  filename = "phylogenetic_tree.png",
  path = "plots/"
)
