#### Script to make heatmap Figure 2

# Author: Ronja Marlonsdotter Sandholm

library(tidyverse)
library(readxl)
library(ggfortify)
library(ggstar)
library(rstatix)
library(ggrepel)
library(DT)

theme <- theme_classic() +
    theme(
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.72),
        axis.text = element_text(color = "black", size = 14),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(margin = margin(5, 0, 0, 0)),
        axis.text.y = element_text(margin = margin(0, 5, 0, 0)),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        # strip.text = element_blank(),
        # legend.position = "bottom",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )

# Growth curve
## Panel A

gc <- read_xlsx("data/AtRMS-02_gc.xlsx", sheet = 1) %>%
    rowwise() %>%
    mutate(
        mean = mean(c_across(contains("Rep")), na.rm = TRUE),
        sd = sd(c_across(contains("Rep")), na.rm = TRUE)
    ) %>%
    ungroup()

panel_a <- ggplot(gc, aes(x = Hours, y = mean)) +
    geom_line(linewidth = 1, aes(color = Substrate)) +
    geom_point(color = "black", fill = "black", size = 2.5, aes(shape = Substrate)) +
    scale_shape_manual(values = c(21, 23)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5) +
    scale_color_manual(values = c("#5d7a60", "#f5c464")) +
    # scale_x_continuous(breaks = seq(0, 250, 50), limits = c(0, 250)) +
    # scale_y_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3)) +
    xlab("Time (h)") +
    ylab("OD600") +
    theme +
    theme(legend.position = "inside")

# panel_a

# Proteomics
prot <- read_tsv("data/AtRMS-02_combined_protein.tsv") %>%
    select(Protein, contains("Intensity"), -contains("LFQ")) %>%
    pivot_longer(-Protein, names_to = "Sample", values_to = "Intensity") %>%
    mutate(
        Sample = str_remove(Sample, " Intensity"),
        Sample = str_replace(Sample, "NaS", "Succinate")
    )

prot_short <- prot %>%
    pivot_wider(names_from = Sample, values_from = Intensity)
datatable(prot_short)
source("impute_normal.R")

imputed <- prot %>%
    mutate(Intensity = log2(Intensity)) %>%
    mutate(Intensity = if_else(is.infinite(Intensity), NA_real_, Intensity)) %>%
    pivot_wider(names_from = Sample, values_from = Intensity) %>%
    filter(
        rowSums(is.na(select(., matches("^Succinate_[1-3]$")))) <= 1 &
            rowSums(is.na(select(., matches("^LMWPE_[1-3]$")))) <= 1 &
            rowSums(is.na(select(., matches("^LMWPE_biofilm_[1-3]$")))) <= 1
    ) %>%
    mutate(across(where(is.numeric), ~ impute_normal(.))) %>%
    mutate(across(where(is.numeric), ~ as.numeric(.)))

dram <- read_tsv("AtRMS-02_annotations.tsv") %>%
    rename(Protein = `...1`) %>%
    mutate(dram_annotation = case_when(rank == "E" ~ "unknown", rank == "D" ~ pfam_hits, rank == "C" ~ kegg_hit)) %>%
    mutate(dram_annotation = str_replace(dram_annotation, "(\\]; ).*", "]")) %>%
    select(Protein, dram_annotation, KO = ko_id)
library(DT)
datatable(dram)

ko <- dram %>%
    select(Protein, KO) %>%
    filter(!is.na(KO))

# write_tsv(ko, "/glittertind/home/ronjasan/scratch/atowneri/bioinformatics/ko_list.tsv")

kegg <- read_tsv("data/KEGG_pathways.tsv")


blast <- read_tsv("data/plastic_blast_annotations.tsv") %>%
    select(Protein = `Gene ID`, gene_blast = Gene)

annotations <- dram %>%
    left_join(kegg, by = "KO") %>%
    left_join(blast, by = "Protein") %>%
    mutate(
        gene = case_when(
            !is.na(gene_blast) ~ gene_blast,
            is.na(gene) & is.na(KO) & grepl("Cytochrome P450", dram_annotation) ~ "cyp",
            is.na(gene) & is.na(KO) & grepl("Rubredoxin", dram_annotation) ~ "Rub",
            is.na(gene) & is.na(KO) & grepl("Flavin-binding monooxygenase-like", dram_annotation) ~ "BVMO",
            is.na(gene) & is.na(KO) & grepl("Alcohol dehydrogenase ", dram_annotation) ~ "ADH",
            is.na(gene) & is.na(KO) & grepl("Aldehyde dehydrogenase family", dram_annotation) ~ "ALDH",
            is.na(gene) & is.na(KO) & grepl("Carboxylesterase family", dram_annotation) ~ "CES",
            is.na(gene) & is.na(KO) & grepl("AtRMS-02_contig_1_1259", Protein) ~ "CES",
            is.na(gene) & is.na(KO) & grepl("AtRMS-02_contig_1_1008", Protein) ~ "CES",
            TRUE ~ gene
        ),
        pathway = case_when(
            is.na(pathway) & grepl("CYP", gene) ~ "alkane",
            is.na(pathway) & grepl("BVMO", gene) ~ "ketone",
            is.na(pathway) & grepl("CES", gene) ~ "ketone",
            is.na(pathway) & grepl("ALDH", gene) ~ "alkane",
            is.na(pathway) & grepl("ADH", gene) ~ "alkane",
            is.na(pathway) & grepl("rub", gene) ~ "alkane",
            is.na(pathway) & grepl("alkMa", gene) ~ "alkane",
            is.na(pathway) & grepl("alkB", gene) ~ "alkane",
            is.na(pathway) & grepl("alkN", gene) ~ "transport",
            TRUE ~ pathway
        ),
        gene = str_remove(gene, "_mcp$")
    ) %>%
    mutate(
        protein = case_when(
            str_detect(gene, "[A-Z0-9]") ~ paste0(toupper(substr(gene, 1, 1)), substr(gene, 2, nchar(gene))),
            TRUE ~ toupper(gene)
        )
    ) %>%
    drop_na(pathway) %>%
    select(-gene_blast)

## Panel B

pca <- imputed %>%
    column_to_rownames("Protein") %>%
    t()

pca_meta <- pca %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    mutate(
        Condition = case_when(
            str_detect(Sample, "Succinate") ~ "Succinate",
            str_detect(Sample, "LMWPE_biofilm") ~ "LMWPE_biofilm",
            TRUE ~ "LMWPE"
        )
    ) %>%
    relocate(Condition, .after = Sample) %>%
    column_to_rownames("Sample")

pca_res <- prcomp(pca, scale. = TRUE)

pcaData <- as.data.frame(pca_res$x[, 1:2]) %>%
    cbind(pca_meta$Condition) %>%
    rename(Condition = `pca_meta$Condition`)

panel_b <- ggplot(pcaData, aes(x = PC1, y = PC2, fill = Condition)) +
    geom_point(size = 3, aes(shape = Condition)) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_manual(values = c("#5b795f", "#3a4d3c", "#f5c464")) +
    theme +
    theme(legend.position = "inside") +
    xlab(paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 2), "%)")) +
    ylab(paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 2), "%)"))

# Panel C + D
imputed_planktonic <- prot %>%
    filter(!str_detect(Sample, "biofilm")) %>%
    mutate(Intensity = log2(Intensity)) %>%
    mutate(Intensity = if_else(is.infinite(Intensity), NA_real_, Intensity)) %>%
    pivot_wider(names_from = Sample, values_from = Intensity) %>%
    filter(
        rowSums(is.na(select(., matches("^Succinate_[1-3]$")))) <= 1 &
            rowSums(is.na(select(., matches("^LMWPE_[1-3]$")))) <= 1
    ) %>%
    mutate(across(where(is.numeric), ~ impute_normal(.))) %>%
    mutate(across(where(is.numeric), ~ as.numeric(.)))

diff_planktonic <- imputed_planktonic %>%
    mutate(mean_LMWPE = rowMeans(select(., contains("LMWPE")), na.rm = TRUE)) %>%
    mutate(mean_Succinate = rowMeans(select(., contains("Succinate")), na.rm = TRUE)) %>%
    mutate(log2FC = mean_LMWPE - mean_Succinate)


impute_test <- prot %>%
    filter(!str_detect(Sample, "LMWPE_[1-3]")) %>%
    mutate(Intensity = log2(Intensity)) %>%
    mutate(Intensity = if_else(is.infinite(Intensity), NA_real_, Intensity)) %>%
    pivot_wider(names_from = Sample, values_from = Intensity) %>%
    filter(Protein == "AtRMS-02_contig_1_1006" | Protein == "AtRMS-02_contig_1_1015" | Protein == "AtRMS-02_contig_1_1019" | Protein == "AtRMS-02_contig_1_10" | Protein == "AtRMS-02_contig_1_751") %>%
    filter(
        rowSums(is.na(select(., matches("^LMWPE_biofilm_[1-3]$")))) <= 1 &
            rowSums(is.na(select(., matches("^Succinate_[1-3]$")))) <= 1
    )

imputed_nas_biofilm <- prot %>%
    filter(!str_detect(Sample, "LMWPE_[1-3]")) %>%
    mutate(Intensity = log2(Intensity)) %>%
    mutate(Intensity = if_else(is.infinite(Intensity), NA_real_, Intensity)) %>%
    pivot_wider(names_from = Sample, values_from = Intensity) %>%
    filter(
        rowSums(is.na(select(., matches("^LMWPE_biofilm_[1-3]$")))) <= 1 &
            rowSums(is.na(select(., matches("^Succinate_[1-3]$")))) <= 1
    ) %>%
    mutate(across(where(is.numeric), ~ impute_normal(.))) %>%
    mutate(across(where(is.numeric), ~ as.numeric(.)))
diff_nas_biofilm <- imputed_nas_biofilm %>%
    mutate(mean_Succinate = rowMeans(select(., contains("Succinate")), na.rm = TRUE)) %>%
    mutate(mean_biofilm = rowMeans(select(., contains("biofilm")), na.rm = TRUE)) %>%
    mutate(log2FC = mean_biofilm - mean_Succinate)

ttest_planktonic <- imputed_planktonic %>%
    pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Intensity") %>%
    filter(!str_detect(Sample, "biofilm")) %>%
    mutate(
        condition = case_when(
            str_detect(Sample, "Succinate") ~ "Succinate",
            str_detect(Sample, "LMWPE") ~ "LMWPE",
            TRUE ~ NA_character_
        ),
        replicate = str_extract(Sample, "[1-3]$")
    ) %>%
    group_by(Protein) %>%
    t_test(
        Intensity ~ condition,
        p.adjust.method = "fdr",
        paired = TRUE
    )

ttest_nas_biofilm <- imputed_nas_biofilm %>%
    pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Intensity") %>%
    mutate(
        condition = case_when(
            str_detect(Sample, "LMWPE_biofilm") ~ "LMWPE_biofilm",
            TRUE ~ "Succinate"
        ),
        replicate = str_extract(Sample, "[1-3]$")
    ) %>%
    group_by(Protein) %>%
    t_test(
        Intensity ~ condition,
        p.adjust.method = "fdr",
        paired = TRUE
    )


volcano_planktonic <- diff_planktonic %>%
    left_join(ttest_planktonic, by = "Protein") %>%
    left_join(annotations, by = "Protein") %>%
    mutate(
        significant = case_when(
            p <= 0.05 & log2FC >= 1 ~ "Upregulated",
            p <= 0.05 & log2FC <= -1 ~ "Downregulated",
            TRUE ~ "Not significant"
        ),
        label = case_when(
            significant != "Not significant" & !is.na(protein) ~ protein,
            TRUE ~ NA_character_
        ),
        log10p = -log10(p)
    )

panel_c <- volcano_planktonic %>%
    ggplot(aes(x = log2FC, y = log10p)) +
    geom_point(size = 1, aes(color = significant)) +
    geom_vline(xintercept = c(-1, 1), col = "darkgray", linetype = "dashed") +
    geom_hline(yintercept = 1.3, col = "darkgray", linetype = "dashed") +
    scale_color_manual(values = c("Upregulated" = "#5d7a60", "Downregulated" = "#f4c465", "Not significant" = "#bdbdbd")) +
    geom_label_repel(aes(label = label, fontface = "bold", segment.color = "black"), box.padding = 1, max.overlaps = Inf, min.segment.length = 0.1, size = 3.5, na.rm = TRUE) +
    labs(x = "log2(LMWPE/Succinate)", y = "-log10(p-value)") +
    theme +
    theme(legend.position = "none")

volcano_nas_biofilm <- diff_nas_biofilm %>%
    left_join(ttest_nas_biofilm, by = "Protein") %>%
    left_join(annotations, by = "Protein") %>%
    mutate(
        significant = case_when(
            p <= 0.05 & log2FC >= 1 ~ "Upregulated",
            p <= 0.05 & log2FC <= -1 ~ "Downregulated",
            TRUE ~ "Not significant"
        ),
        label = case_when(
            significant != "Not significant" & !is.na(protein) ~ protein,
            TRUE ~ NA_character_
        ),
        log10p = -log10(p)
    )


panel_d <- volcano_nas_biofilm %>%
    ggplot(aes(x = log2FC, y = log10p)) +
    geom_point(size = 1, aes(color = significant)) +
    geom_vline(xintercept = c(-1, 1), col = "darkgray", linetype = "dashed") +
    geom_hline(yintercept = 1.3, col = "darkgray", linetype = "dashed") +
    scale_color_manual(values = c("Upregulated" = "#5d7a60", "Downregulated" = "#f4c465", "Not significant" = "#bdbdbd")) +
    geom_label_repel(aes(label = label, fontface = "bold", segment.color = "black"), box.padding = 1, max.overlaps = Inf, min.segment.length = 0.1, size = 3.5, na.rm = TRUE) +
    labs(x = "log2(LMWPE_biofilm/Succinate)", y = "-log10(p-value)") +
    theme +
    theme(legend.position = "none")

# Panel E
## Heatmap of selected proteins

heat <- prot %>%
    mutate(Intensity = log2(Intensity)) %>%
    mutate(Intensity = if_else(is.infinite(Intensity), NA_real_, Intensity)) %>%
    pivot_wider(names_from = Sample, values_from = Intensity) %>%
    right_join(annotations, by = "Protein") %>%
    drop_na(gene) %>%
    mutate(
        mean_Succinate = rowMeans(select(., matches("^Succinate_[1-3]$")), na.rm = TRUE),
        mean_LMWPE = rowMeans(select(., matches("^LMWPE_[1-3]$")), na.rm = TRUE),
        mean_biofilm = rowMeans(select(., matches("^LMWPE_biofilm_[1-3]$")), na.rm = TRUE),
        protein = paste0(protein, "_", str_extract(Protein, "(?<=_)[^_]+$"))
    ) %>%
    filter(!(mean_Succinate > mean_LMWPE & mean_Succinate > mean_biofilm) | is.na(mean_Succinate)) %>%
    filter(!(is.na(mean_Succinate) & is.na(mean_LMWPE) & is.na(mean_biofilm))) %>%
    mutate(
        order = case_when(
            grepl("FadL", protein) ~ paste0(1, "_", protein),
            grepl("AlkN", protein) ~ paste0(2, "_", protein),
            grepl("AlkR", protein) ~ paste0(1, "_", protein),
            grepl("RubA", protein) ~ paste0(2, "_", protein),
            grepl("RubB", protein) ~ paste0(2, "_", protein),
            grepl("AlkB", protein) ~ paste0(3, "_", protein),
            grepl("AlkMa", protein) ~ paste0(3, "_", protein),
            grepl("AlmA", protein) ~ paste0(4, "_", protein),
            grepl("ADH", protein) ~ paste0(5, "_", protein),
            grepl("YiaY", protein) ~ paste0(5, "_", protein),
            grepl("ALDH", protein) ~ paste0(6, "_", protein),
            grepl("BVMO", protein) ~ paste0(1, "_", protein),
            grepl("CES", protein) ~ paste0(2, "_", protein),
            grepl("ACS_", protein) ~ paste0(3, "_", protein),
            grepl("ACSL_", protein) ~ paste0(1, "_", protein),
            grepl("ACAD", protein) ~ paste0(2, "_", protein),
            grepl("GCDH", protein) ~ paste0(3, "_", protein),
            grepl("FadJ", protein) ~ paste0(4, "_", protein),
            grepl("FadB", protein) ~ paste0(5, "_", protein),
            grepl("FadA", protein) ~ paste0(6, "_", protein),
            grepl("ACAT", protein) ~ paste0(7, "_", protein)
        ),
        pathway = case_when(
            pathway == "transport" ~ paste0(1, "_", pathway),
            pathway == "alkane" ~ paste0(2, "_", pathway),
            pathway == "ketone" ~ paste0(3, "_", pathway),
            pathway == "fatty acid" ~ paste0(4, "_", pathway)
        )
    ) %>%
    pivot_longer(cols = matches("^(Succinate|LMWPE|LMWPE_biofilm)_[1-3]$"), names_to = "Sample", values_to = "Intensity") %>%
    mutate(
        rep = str_extract(Sample, "[1-3]$"),
        condition = case_when(
            str_detect(Sample, "Succinate") ~ "Succinate",
            str_detect(Sample, "LMWPE_biofilm") ~ "LMWPE_biofilm",
            TRUE ~ "LMWPE"
        ),
        # Fix condition ordering
        condition = factor(condition, levels = c("Succinate", "LMWPE", "LMWPE_biofilm")),
        rep = factor(rep, levels = c("1", "2", "3"))
    ) %>%
    mutate(
        fill = case_when(
            Intensity < 15 ~ "14",
            Intensity >= 15 & Intensity < 17 ~ "15-17",
            Intensity >= 17 & Intensity < 19 ~ "17-19",
            Intensity >= 19 & Intensity < 21 ~ "19-21",
            Intensity >= 21 & Intensity < 23 ~ "21-23",
            Intensity >= 23 & Intensity < 25 ~ "23-25",
            Intensity >= 25 ~ "26",
            TRUE ~ NA
        )
    ) %>%
    drop_na(order)


panel_e <- ggplot(heat, aes(x = rep, y = order, fill = Intensity)) +
    geom_tile(color = "white", aes(height = 0.9, width = 0.9)) +
    scale_fill_gradientn(colors = colorRampPalette(c("#c7463e", "#f0755d", "#f5c464", "#f5d9a1", "#8ab891", "#5b795f", "#3a4d3c"))(10), na.value = "lightgray") +
    facet_grid(pathway ~ condition, scales = "free", space = "free") +
    theme +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12)
    ) +
    labs(fill = "log2(int)")

# Combine panels
library(cowplot)
fig2a_d <- plot_grid(
    panel_a, panel_b, panel_c, panel_d,
    label_size = 20,
    labels = "AUTO",
    ncol = 2
)


fig2 <- plot_grid(
    fig2a_d, panel_e,
    label_size = 20,
    labels = c("", "E"),
    ncol = 2,
    rel_widths = c(1, 0.7)
)
