# this script use the given test_table contain the p values and
# fdr and logfc as given

library("DESeq2")
library("pheatmap")
library("viridis")
library("GGally")
library("ggrepel") # for figures
# library("heatmaply")


# filter rows that have no or nearly no information about
# the amount of gene expression
# zeros are are already removed from testTable.tab
# this will also filter those with sum of counts 10 or less
nrow(test_table) # 29970

# calculate the sum of counts across samples
# and filter to keep only those with a sum at least 10
test_table <- test_table %>%
  rowwise() %>%
  mutate(sum = sum(c_across(S01_P17040_unstim_lo:S34_P16088_stimul_lo),
                   na.rm = T)) %>% ungroup() %>%
  filter(sum >= 10)
nrow(test_table) # 18169
# has taken it down to 18169 rows
# 29970 - 18169 = 11801 have counts < 10

# p-values set to NA if
#   - all samples have zero counts, the baseMean  FC, p and adj p
#     will be zero
#   - If a row contains a sample with an extreme count outlier then
#     the p value and adjusted p value will be set to NA. outliers
#     are detected by Cook’s distance.
#   - If a row is automatically filtered for having a
#     low mean normalized count, then only the adjusted p value
#     will be NA.
# Remove and rows where FDR is NA
test_table <- test_table %>%
  filter(!is.na(`FDR.Unstimulated/Stimulated`))
# has taken it down to 18129 rows
# 18169 - 18129 = 40 have NA in FDR  for outliers of low means


# VST transform the raw counts for putting in PCA

# make the count matrix
# select just the counts
cols <- coldata$names
countdata <- test_table  %>%
  select(all_of(cols)) %>%
  round(1)
row.names(countdata) <- test_table$gene_id
head(countdata, 3)

# check that the order of the columns in countdata is the same as the order
# of the rows in coldata
names(countdata) == coldata$names

# build a DESeqDataSet from a count matrix and table of sample info

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ donor + treatment)


# VST transformation
# methods like PCA work best if data homoskedastic. but with counts, the
# variance increases with the mean
# we apply a *variance stabilizing transformation* (VST) for
# negative binomial data with a dispersion-mean trend
# [@Anders2010Differential],
dds_vst <- vst(dds, blind = FALSE)

# transpose matrix so genes are in columns, samples in rows
# this is need to calculaute the distance between samples and
# for PCA
dds_vst_t <- t(assay(dds_vst))

# Heatmap of sample-to-sample distances using the variance
# stabilizing transformed values.
# sample distances
sampleDists <- dist(dds_vst_t)
sampleDistMatrix <- as.matrix( sampleDists)
rownames(sampleDistMatrix) <- paste(dds_vst$treatment,
                                    dds_vst$donor)

colnames(sampleDistMatrix) <- NULL
p1 <- pheatmap(sampleDistMatrix,
               clustering_distance_rows = sampleDists,
               clustering_distance_cols = sampleDists,
               col = viridis(10))
ggsave("reports/figures/sample_distance_heatmap-given.png",
       plot = p1,
       device = "png",
       width = 6,
       height = 7,
       units = "in",
       dpi = 300)


# PCA
# carry out PCA with out scaling
pca <- dds_vst_t %>%
  prcomp(scale. = FALSE)
# examine amount of variance captured
summary(pca)[["importance"]][,1:10]
dat <-  data.frame(pca$x, sample = rownames(dds_vst_t))
dat <- dat %>%
  extract(sample,
          c("samp", "donor", "treatment", "idoresp"),
          "(S[0-9]{2})_(P[0-9]{5})_([a-z]{6})_([a-z]{2})")

ggplot(dat,
       aes(x = PC1, y = PC2, colour = donor, shape = treatment)) +
  geom_point(size = 3)

dat %>% filter(treatment == "stimul") %>%
  ggplot(aes(x = PC1, y = PC2, colour = idoresp)) +
  geom_point(size = 3)


p2 <- dat %>%
  select(PC1:PC6, donor, treatment) %>%
  ggpairs(aes(colour = donor, shape = treatment),
          upper = NULL,
          lower = list(continuous = wrap("points", size = 2)),
          columns = 1:6,
          diag = NULL) +
  theme_minimal()
ggsave("reports/figures/pca_pairs.png",
       plot = p2,
       device = "png",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)

# heatmap clustering

# Matrix of selected data with a log FC at least 4.5 and FDR 0.01
tmp_FDR_0.01 <- test_table %>%
  filter(`FDR.Unstimulated/Stimulated` <= 0.01) %>%
  filter(abs(`logFC.Unstimulated/Stimulated`) >= 4.5) %>%
  select(starts_with("TMP_"),
         gene_name,
         `logFC.Unstimulated/Stimulated`)
# there are 117 gene_names with a high FC and low FDR
length(unique(tmp_FDR_0.01$gene_name))

tmp_FDR_0.01_mat <- tmp_FDR_0.01 %>%
  select(-gene_name,
         -`logFC.Unstimulated/Stimulated`) %>%
  as.matrix()

row.names(tmp_FDR_0.01_mat) <- tmp_FDR_0.01$gene_name


n_treatment_clusters <- 2
n_gene_clusters <- 2


heatmaply::heatmaply(tmp_FDR_0.01_mat,
                     scale = "row",
                     grid_color = "white",
                     hide_colorbar = TRUE,
                     k_col = n_treatment_clusters,
                     k_row = n_gene_clusters,
                     label_names = c("Gene", "Sample", "TMP"),
                     fontsize_row = 6, fontsize_col = 10,
                     labCol = toupper(colnames(tmp_FDR_0.01_mat)),
                     labRow = rownames(tmp_FDR_0.01_mat),
                     heatmap_layers = theme(axis.line = element_blank()))



# list of genes up and down regulated
# note: Unstimulated/Stimulated thus if unstim > stim lofFC is +'ve
# if unstim < stim lofFC is -'ve

# sig up regulated in unstimulated
# FDR < 0.05 and FC > 0
unstim_sig_up <- test_table %>%
  filter(`logFC.Unstimulated/Stimulated` > 0) %>%
  filter(`FDR.Unstimulated/Stimulated` < 0.05) %>%
  arrange(`FDR.Unstimulated/Stimulated`,
          dplyr::desc(`logFC.Unstimulated/Stimulated`))
write(unstim_sig_up$gene_name, "reports/unstim_sig_up.txt")

# sig up regulated in stimulated
# FDR < 0.05 and FC < 0
stimul_sig_up <- test_table %>%
  filter(`logFC.Unstimulated/Stimulated` < 0) %>%
  filter(`FDR.Unstimulated/Stimulated` < 0.05) %>%
  arrange(`FDR.Unstimulated/Stimulated`,
          `logFC.Unstimulated/Stimulated`)
write(stimul_sig_up$gene_name, "reports/stimul_sig_up.txt")

# volcano plot with their values
# note: Unstimulated/Stimulated thus if unstim > stim lofFC is +'ve
# if unstim < stim lofFC is -'ve
# add a indicator for impt genes - those with FC >=3 and FDR <=0.01
# set FDR values less then 1e-100 to 0
test_table <- test_table %>%
  mutate(impt = abs(`logFC.Unstimulated/Stimulated`) >= 3 &
           `FDR.Unstimulated/Stimulated` <= 0.01,
         fdr = case_when(`FDR.Unstimulated/Stimulated` < 1e-100 ~ 0,
                         `FDR.Unstimulated/Stimulated` >= 1e-100 ~ `FDR.Unstimulated/Stimulated`))



test_table %>%
  select(gene_name,
         impt,
         `logFC.Unstimulated/Stimulated`,
         fdr) %>%
  ggplot(aes(x = `logFC.Unstimulated/Stimulated`,
             y = -log10(fdr),
             colour = impt)) +
  # scale_y_continuous(limits = c(0 , 200)) +
  geom_point() +
  scale_colour_manual(values = c("#999999", "#D9017A")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_vline(xintercept = -3, linetype = "dashed") +
  scale_y_continuous(name = "-log10[FDR])") +
  scale_x_continuous(name = "log2 Foldchange Unstimulated/stimulated",
                     limits = c(-10, 15)) +
  annotate("text", x = -7.5,  y = 80,
           label = "Stimulated\nsignificantly higher") +
  annotate("text", x = 12,  y = 80,
           label = "Unstimulated\nsignificantly higher") +
  annotate("text", x = -7.5,  y = 0,
           label = "NS") +
  geom_text_repel(data = subset(test_table,
                                impt == TRUE),
                  aes(label = gene_name),
                  size = 3, max.overlaps = 30) +
  guides(colour = "none") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black",
                                    fill = NA))








# go analysis
# computational genomics
