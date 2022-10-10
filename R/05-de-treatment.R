####################        NEDDS EDITING         ##################



##################################################################
# Differential expression between high and low IDO transcripts
##################################################################
library("DESeq2")



###############################################################################
# Differential expression analysis for Treatment
###############################################################################


# count matrix
count_data <- test_table %>%
  select(S01_P17040_unstim_lo:S34_P16088_stimul_lo) %>%
  as.matrix()
row.names(count_data) <- test_table$gene_id


# coldata


# define the experimental setup

#define the design formula
design_treat <- "~ treatment"


# create a DESeq dataset object from the count matrix and the colData
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = as.formula(design_treat))
# print dds object to see the contents
print(dds) # dim: 25324

# estimate size factors to normalize the counts and dispersion values
# compute GLM model based on the experimental design formula
dds <- DESeq(dds)

# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 40 genes
# -- DESeq argument 'minReplicatesForReplace' = 7
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing


# compute the contrast for the 'treatment' variable where "unstim"
# is the base
# this means +'ve FC are upreg in stim relative to unstim
de_results = results(dds, contrast = c("treatment", "stimul", "unstim"))

# sort results by increasing p-value
# de_results <- de_results[order(de_results$pvalue),]

print(de_results)
###############################################################################
# Volcano plot for Differential expression analysis for Treatment
###############################################################################

de_results_vol <- de_results %>%
  data.frame(gene_id = row.names(de_results)) %>%
  filter(!is.na(padj))
# dim(de_results_vol) 17468
length(unique(de_results_vol$gene_id))
AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)
# Annotation
annotations <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                     keys = de_results_vol$gene_id,
                                     columns = c("SYMBOL",
                                                 "ENTREZID",
                                                 "GENENAME",
                                                 "ENSEMBL",
                                                 "GENETYPE"),
                                     keytype = "ENSEMBL")
# dim(annotations) 17558 we have 17558 - 17468 = 90 duplicates

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations$ENSEMBL) == FALSE)

# Return only the non-duplicated genes using indices
annotations <- annotations[non_duplicates_idx, ]
# dim(annotations) 17468

# merge annotations in to results
names(annotations)[1] <- "gene_id"
de_results_vol <- de_results_vol %>% merge(annotations, by = "gene_id")

# genetypes
de_results_vol %>% group_by(GENETYPE) %>% summarise(n = length(GENETYPE))
# GENETYPE           n
# 1 ncRNA           1236
# 2 protein-coding 11807
# 3 pseudo           485
# 4 snoRNA             2
# 5 unknown            2
# 6 NA              3936

# add a indicator for impt genes - those with FC >=3 and FDR <=0.01
# set FDR values less then 1e-50 to 0
de_results_vol <- de_results_vol %>%
  mutate(sig = padj <= 0.05,
         impt = abs(log2FoldChange) >= 3 &
           padj <= 0.01,
         fdr = case_when(padj < 1e-50 ~ 0,
                         padj >= 1e-50 ~ padj))

de_results_vol %>% group_by(GENETYPE, impt) %>% summarise(n = length(GENETYPE))
# GENETYPE       impt      n
# 1 ncRNA          FALSE  1201
# 2 ncRNA          TRUE     35
# 3 protein-coding FALSE 11615
# 4 protein-coding TRUE    192
# 5 pseudo         FALSE   475
# 6 pseudo         TRUE     10
# 7 snoRNA         FALSE     2
# 8 unknown        FALSE     2
# 9 NA             FALSE  3822
# 10 NA             TRUE    114

de_results_vol %>% group_by(GENETYPE, sig) %>% summarise(n = length(GENETYPE))
#   GENETYPE       sig       n
# 1 ncRNA          FALSE   986
# 2 ncRNA          TRUE    250
# 3 protein-coding FALSE 10253
# 4 protein-coding TRUE   1554
# 5 pseudo         FALSE   442
# 6 pseudo         TRUE     43
# 7 snoRNA         FALSE     2
# 8 unknown        FALSE     2
# 9 NA             FALSE  3194
# 10 NA             TRUE    742


# write to file
write_csv(de_results_vol,
          file = "data-processed/diff-expr-by-treatment.csv")

# genes that are significantly up regulated in the stimulation treatment
de_results_sig_up_stimul <- de_results_vol %>%
  filter(sig & log2FoldChange > 1)
# Note 252/965 sig up results have no annotation
# table(is.na(de_results_sig_up_stimul$SYMBOL))
# table(is.na(de_results_sig_up_stimul$ENTREZID))
# table(is.na(de_results_sig_up_stimul$GENENAME))
write_csv(de_results_sig_up_stimul,
          file = "data-processed/diff-expr-by-treatment-sig-up-stim.csv")

# genes that are significantly down regulated in the stimulation treatment
de_results_sig_down_stimul <- de_results_vol %>%
  filter(sig & log2FoldChange < -1)
# Note 208/734  sig down results have no annotation
# table(is.na(de_results_sig_down_stimul$SYMBOL))
# table(is.na(de_results_sig_down_stimul$ENTREZID))
# table(is.na(de_results_sig_down_stimul$GENENAME))
write_csv(de_results_sig_down_stimul,
          file = "data-processed/diff-expr-by-treatment-sig-down-stim.csv")


# those with FC >=3 and FDR <=0.01 UP
de_results_impt_up_stimul <- de_results_vol %>%
  filter(impt & log2FoldChange > 0)
# Note 84/273 impt up results have no annotation
# table(is.na(de_results_impt_up_stimul$SYMBOL))
# table(is.na(de_results_impt_up_stimul$ENTREZID))
# table(is.na(de_results_impt_up_stimul$GENENAME))
write_csv(de_results_impt_up_stimul,
          file = "data-processed/diff-expr-by-treatment-impt-up-stim.csv")

# genes that "important" i.e, those with FC >=3 and FDR <=0.01 DOWN
de_results_impt_down_stimul <- de_results_vol %>%
  filter(impt & log2FoldChange < 0)
# Note 30/78  impt down results have no annotation
# table(is.na(de_results_impt_down_stimul$SYMBOL))
# table(is.na(de_results_impt_down_stimul$ENTREZID))
# table(is.na(de_results_impt_down_stimul$GENENAME))
write_csv(de_results_impt_down_stimul,
          file = "data-processed/diff-expr-by-treatment-impt-down-stim.csv")


# volcano Protein coding
# 11807
de_results_vol %>% filter(GENETYPE == "protein-coding") %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(fdr),
             colour = impt)) +
  geom_point() +
  scale_colour_manual(values = c("#999999", "#D9017A")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_vline(xintercept = -3, linetype = "dashed") +
  scale_y_continuous(name = "-log10[FDR])",
                     limits = c(0, 60)) +
  scale_x_continuous(name = "log2 Foldchange Stimulated/Unstimulated",
                     limits = c(-10, 15)) +
  annotate("text", x = -10,  y = 35,
           label = "Unstimulated\nsignificantly higher\nFDR < 0.01, log2FC ≤ -3",
           hjust = 0) +
  annotate("text", x = 15,  y = 35,
           label = "Stimulated\nsignificantly higher\nFDR < 0.01, log2FC ≥ 3",
           hjust = 1) +
  annotate("text", x = -9.5,  y = 0,
           label = "NS",
           hjust = 0) +
  geom_text_repel(data = subset(de_results_vol,
                                impt == TRUE & GENETYPE == "protein-coding"),
                  aes(label = SYMBOL),
                  size = 2, max.overlaps = 100) +
  guides(colour = "none") +
  ggtitle("Volcano plot for protein coding genes (11807)",
          subtitle = "Pink = fdr < 0.01 and a log2 FC of at least 3 (192 genes)") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black",
                                    fill = NA))

ggsave(filename = "reports/figures/volcano-protein-coding.png",
       units = "in",
       height = 8,
       width = 10,
       dpi = 300,
       device = "png")


# volcano ncRNA
# 1236
de_results_vol %>% filter(GENETYPE == "ncRNA") %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(fdr),
             colour = impt)) +
  geom_point() +
  scale_colour_manual(values = c("#999999", "#D9017A")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_vline(xintercept = -3, linetype = "dashed") +
  scale_y_continuous(name = "-log10[FDR])",
                     limits = c(0, 60)) +
  scale_x_continuous(name = "log2 Foldchange Stimulated/Unstimulated",
                     limits = c(-10, 15)) +
  annotate("text", x = -10,  y = 35,
           label = "Unstimulated\nsignificantly higher\nFDR < 0.01, log2FC ≤ -3",
           hjust = 0) +
  annotate("text", x = 15,  y = 35,
           label = "Stimulated\nsignificantly higher\nFDR < 0.01, log2FC ≥ 3",
           hjust = 1) +
  annotate("text", x = -9.5,  y = 0,
           label = "NS",
           hjust = 0) +
  geom_text_repel(data = subset(de_results_vol,
                                impt == TRUE & GENETYPE == "ncRNA"),
                  aes(label = SYMBOL),
                  size = 2, max.overlaps = 100) +
  guides(colour = "none") +
  ggtitle("Volcano plot for ncRNA (1236)",
          subtitle = "Pink = fdr < 0.01 and a log2 FC of at least 3 (35)") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black",
                                    fill = NA))

ggsave(filename = "reports/figures/volcano-ncRNA.png",
       units = "in",
       height = 8,
       width = 10,
       dpi = 300,
       device = "png")


# volcano unannotated
# 3936
de_results_vol %>% filter(is.na(GENETYPE)) %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(fdr),
             colour = impt)) +
  geom_point() +
  scale_colour_manual(values = c("#999999", "#D9017A")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_vline(xintercept = -3, linetype = "dashed") +
  scale_y_continuous(name = "-log10[FDR])",
                     limits = c(0, 60)) +
  scale_x_continuous(name = "log2 Foldchange Stimulated/Unstimulated",
                     limits = c(-10, 15)) +
  annotate("text", x = -10,  y = 35,
           label = "Unstimulated\nsignificantly higher\nFDR < 0.01, log2FC ≤ -3",
           hjust = 0) +
  annotate("text", x = 15,  y = 35,
           label = "Stimulated\nsignificantly higher\nFDR < 0.01, log2FC ≥ 3",
           hjust = 1) +
  annotate("text", x = -9.5,  y = 0,
           label = "NS",
           hjust = 0) +
  geom_text_repel(data = subset(de_results_vol,
                                impt == TRUE & is.na(GENETYPE)),
                  aes(label = gene_id),
                  size = 2, max.overlaps = 100) +
  guides(colour = "none") +
  ggtitle("Volcano plot for unannotated (3936)",
          subtitle = "Pink = fdr < 0.01 and a log2 FC of at least 3 (114)") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black",
                                    fill = NA))

ggsave(filename = "reports/figures/volcano-unannotated.png",
       units = "in",
       height = 8,
       width = 10,
       dpi = 300,
       device = "png")



##################################################################
#GO Terms for differentially expressed between stim and unstim
##################################################################

library(gprofiler2)


# significantly up regulated genes
go_results_up <- gost(query = de_results_sig_up_stimul$gene_id,
                      organism = "hsapiens",
                      ordered_query = FALSE,
                      multi_query = FALSE,
                      significant = TRUE,
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


# significantly up regulated genes
go_results_down <- gost(query = de_results_sig_down_stimul$gene_id,
                        organism = "hsapiens",
                        ordered_query = FALSE,
                        multi_query = FALSE,
                        significant = TRUE,
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

gostplot(go_results_down,
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



