

##################################################################
# Differential expression between high and low IDO transcripts and ATP
##################################################################
library("DESeq2")

# high IDO transcripts and high ATP (P12316 and P17083) and the other donors?

#################################################################
# All data
#################################################################



# count matrix
count_data_idotrans <- test_table %>%
  select(S01_P17040_unstim_lo:S34_P16088_stimul_lo) %>%
  as.matrix()
row.names(count_data_idotrans) <- test_table$gene_id


# define the experimental setup

# define the design formula (different)
design_idotrans <- "~ idotrans"

coldata_idotrans <- coldata

# create a DESeq dataset object from the count matrix and the coldata_idotrans
dds_idotrans <- DESeqDataSetFromMatrix(countData = count_data_idotrans,
                              colData = coldata_idotrans,
                              design = as.formula(design_idotrans))
# print dds_idotrans object to see the contents
# print(dds_idotrans) # dim: 25324

# for each gene, count the total number of reads for that gene
# in all samples and remove those that don't have at least 2 reads
dds_idotrans <- dds_idotrans[ rowSums(DESeq2::counts(dds_idotrans)) > 1, ]
print(dds_idotrans) # dim: 25324

# estimate size factors to normalize the counts and dispersion values
# compute GLM model based on the experimental design formula
dds_idotrans <- DESeq(dds_idotrans) #replacing outliers and refitting for 32 genes

# compute the contrast for the idotrans variable where "low"
# is the base
# this means +'ve FC are upreg in high relative to low
de_results_idotrans  <- results(dds_idotrans,
                                contrast = c("idotrans", "high", "low"))

###############################################################################
# Volcano data for Differential expression analysis for (P12316 and P17083)
# vs other donors: stim only
###############################################################################

de_results_idotrans_vol <- de_results_idotrans %>%
  data.frame(gene_id = row.names(de_results_idotrans)) %>%
  filter(!is.na(padj))
# dim(de_results_idotrans_vol) 16480
# length(unique(de_results_idotrans_vol$gene_id)) 16480
# AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)
# Annotation
annotations_idotrans <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                     keys = de_results_idotrans_vol$gene_id,
                                     columns = c("SYMBOL",
                                                 "ENTREZID",
                                                 "GENENAME",
                                                 "ENSEMBL",
                                                 "GENETYPE"),
                                     keytype = "ENSEMBL")
# dim(annotations_idotrans) 16568        we have 88 duplicates

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_idotrans$ENSEMBL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_idotrans <- annotations_idotrans[non_duplicates_idx, ]
# dim(annotations_idotrans) 16480


# merge annotations in to results
names(annotations_idotrans)[1] <- "gene_id"
de_results_idotrans_vol <- de_results_idotrans_vol %>%
  merge(annotations_idotrans, by = "gene_id")


# add a indicator for impt genes - those with FC > 0 and FDR <=0.05
# set FDR values less then 1e-50 to 0
de_results_idotrans_vol <- de_results_idotrans_vol %>%
  mutate(sig = padj <= 0.05,
         impt = abs(log2FoldChange) > 0 &
           padj <= 0.05,
         fdr = case_when(padj < 1e-50 ~ 0,
                         padj >= 1e-50 ~ padj))
write_csv(de_results_idotrans_vol,
          file = "data-processed/de_ido_resp_idotrans_results.csv")

de_results_idotrans_vol %>%
  group_by(sig) %>%
  summarise(n = length(GENETYPE)) %>%
  kable(row.names = FALSE) %>%
  kable_styling(font_size = 12)

de_results_idotrans_vol %>%
  group_by(GENETYPE, sig) %>%
  summarise(n = length(GENETYPE)) %>%
  kable(row.names = FALSE) %>%
  kable_styling(font_size = 12)


de_results_idotrans_vol %>%
  filter(sig) %>%
  select(SYMBOL, gene_id, log2FoldChange, padj, GENENAME) %>%
  kable(row.names = FALSE, digits = 4) %>%
  kable_styling(font_size = 12)


# genes that are significantly up regulated in the high DO transcripts pair
de_results_sig_up_idotrans_high <- de_results_idotrans_vol %>%
  filter(sig & log2FoldChange > 0)

write_csv(de_results_sig_up_idotrans_high,
          file = "data-processed/diff-expr-by-idotrans-sig-up-high.csv")

# genes that are significantly down regulated in the high DO transcripts pair
de_results_sig_down_idotrans_high <- de_results_idotrans_vol %>%
  filter(sig & log2FoldChange < 0)

write_csv(de_results_sig_down_idotrans_high,
          file = "data-processed/diff-expr-by-idotrans-sig-down-high.csv")


# volcano
de_results_idotrans_vol  %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(fdr),
             colour = impt)) +
  geom_point() +
  scale_colour_manual(values = c("#999999", "#D9017A")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  scale_y_continuous(name = "-log10[FDR])",
                     limits = c(0, 5)) +
  scale_x_continuous(name = "log2 Foldchange Stimulated/Unstimulated",
                     limits = c(-5, 5)) +
  annotate("text", x = -5,  y = 4.5,
           label = "Down regulated in P12316 and P17083 \nFDR < 0.05, log2FC < 0",
           hjust = 0) +
  annotate("text", x = 5,  y = 4.5,
           label = "Up regulated in P12316 and P17083 \nFDR < 0.05, log2FC > 0",
           hjust = 1) +
  annotate("text", x = -5,  y = 0,
           label = "NS",
           hjust = 0) +
  geom_text_repel(data = subset(de_results_idotrans_vol,
                                impt == TRUE),
                  aes(label = SYMBOL),
                  size = 3, max.overlaps = 100) +
  guides(colour = "none") +
  ggtitle("Volcano plot comparing high IDO transcripts and high\n ATP (P12316 and P17083) and the others",
          subtitle = "Pink = fdr < 0.05") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black",
                                    fill = NA))

ggsave(filename = "reports/figures/volcano-idotrans.png",
       units = "in",
       height = 8,
       width = 10,
       dpi = 300,
       device = "png")


##################################################################
#GO Terms for differentially expressed between low and high
##################################################################

library(gprofiler2)


# significantly up regulated genes
go_results_up <- gost(query = de_results_sig_up_idotrans_high$gene_id,
                      organism = "hsapiens",
                      ordered_query = FALSE,
                      multi_query = FALSE,
                      significant = FALSE,
                      exclude_iea = FALSE,
                      measure_underrepresentation = FALSE,
                      evcodes = FALSE,
                      user_threshold = 0.05,
                      correction_method = "g_SCS",
                      domain_scope = "annotated",
                      custom_bg = NULL,
                      numeric_ns = "",
                      sources = NULL,
                      as_short_link =  FALSE)

gostplot(go_results_up,
         capped = TRUE,
         interactive = TRUE)



# Available data sources and their abbreviations are:
#
#   Gene Ontology (GO or by branch GO:MF, GO:BP, GO:CC)
# KEGG (KEGG)
# Reactome (REAC)
# WikiPathways (WP)
# TRANSFAC (TF)
# miRTarBase (MIRNA)
# Human Protein Atlas (HPA)
# CORUM (CORUM)
# Human phenotype ontology (HP)



