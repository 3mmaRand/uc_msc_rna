library("tidyverse")
library("pheatmap")
library("corrplot")
library("DESeq2")
library("viridis")
library("GGally")
library("ggrepel")
library("patchwork")
# library("org.Hs.eg.db")
# library("AnnotationDbi")
# library("heatmapply")


##################################################################
#
counts <- test_table %>%
  select(S01_P17040_unstim_lo:S34_P16088_stimul_lo,
         code_length) %>% as.matrix()

row.names(counts) <- test_table$gene_id

# # Computing CPM
#
# cpm <- apply(subset(counts, select = c(-code_length)), 2,
#              function(x) x/sum(as.numeric(x)) * 10^6)
#
# # Computing RPKM
# # create a vector of gene lengths
# geneLengths <- as.vector(subset(counts, select = c(code_length)))
#
# # compute rpkm
# rpkm <- apply(X = subset(counts, select = c(-code_length)),
#               MARGIN = 2,
#               FUN = function(x) {
#                 10^9 * x / geneLengths / sum(as.numeric(x))
#               })
# # Computing TPM
# # find gene length normalized values
# rpk <- apply( subset(counts, select = c(-code_length)), 2,
#               function(x) x/(geneLengths/1000))
# #normalize by the sample size using rpk values
# tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
# # LIVERPOOL HAVE USED code_length not end - start to calculate tmp
# # so i have too

##################################################################

# CPM, RPKM/FPKM, TPM do not account for the library composition,
# aka relative size of the compared transcriptomes
# but DESeq2 does


# ########################################################################
# # Exploratory analysis of the read count table
# ############################################################################
# # compute the variance of each gene across samples
# # note that genes with the most variance might not be those
# # with the great FC between treatments
# V <- apply(tpm, 1, var)
# # sort the results by variance in decreasing order
# # and select the top 100 genes
# selectedGenes <- names(V[order(V, decreasing = T)][1:100])
#
# anno_coldata <- coldata %>% select(-fcido, -names)
#
# # note heat map is of the genes with the great variance.
# # variance could be lower and still significant
# pheatmap(tpm[selectedGenes,], scale = 'row',
#          show_rownames = FALSE,
#          annotation_col = anno_coldata,
#          width = 6,
#          height = 9,
#          filename = "reports/figures/clust-heatmap.png",
#          main = "Top 100 most variable genes (by TPM)")
# # two of the unstim seem to cluster with the stimul in this
# # view
# # dev.off()
# ######################################################################
# # PCA on all genes TMP
# ######################################################################
#
# M <- t(tpm)
# M <- log2(M + 1)
# #compute PCA
#
# # carry out PCA
# pca <- M %>%
#   prcomp()
# # examine amount of variance captured
# pcvar <- summary(pca)[["importance"]][2,1:2]  %>% round(4) * 100
#
# # put sample name and pc score in a dataframe
# dat <-  data.frame(pca$x, names = rownames(M))
#
# # add the metadata from coldata for annotation of plot
# dat <- dat %>% merge(coldata, by = "names")
#
# pca_tmp <- ggplot(dat,
#                   aes(x = PC1, y = PC2, colour = donor, shape = treatment)) +
#   geom_point(size = 3) +
#   xlab(paste("PC1: ", pcvar[1],"%")) +
#   ylab(paste("PC2: ", pcvar[2],"%")) +
#
#   theme_classic()
#
#
# pca_tmp2 <- ggplot(dat,
#                   aes(x = PC1, y = PC2,
#                       colour = factor(ATPrank),
#                       shape = treatment)) +
#   geom_point(size = 3) +
#   xlab(paste("PC1: ", pcvar[1],"%")) +
#   ylab(paste("PC2: ", pcvar[2],"%")) +
#   theme_classic()
#
# pca <- pca_tmp + pca_tmp2 +
#   plot_annotation(title = "PCA on log2(TPM + 1): 25324 genes (counts > 1)")
#
# ggsave(pca,
#        filename = "reports/figures/pca-log2tmp-all.png",
#        height = 5,
#        width = 10)
#
#
#
# pca_pairs_tmp <- dat %>%
#   select(PC1:PC6, donor, treatment) %>%
#   ggpairs(aes(colour = donor, shape = treatment),
#           upper = NULL,
#           lower = list(continuous = wrap("points", size = 2)),
#           columns = 1:6,
#           diag = NULL,
#           legend = grab_legend(pca_tmp)) +
#   ggtitle("Pairwise plots of first 6 PCs",
#           subtitle = "PCA on log2(TPM + 1): 25324 genes (counts > 1)") +
#   theme_minimal()
#
# ggsave("reports/figures/pca_pairs_log2tmp.png",
#        plot = pca_pairs_tmp,
#        device = "png",
#        width = 8,
#        height = 6,
#        units = "in",
#        dpi = 300)
#
#
# pca_pairs_tmp2 <- dat %>%
#   select(PC1:PC6, ATPrank, treatment) %>%
#   ggpairs(aes(colour = factor(ATPrank), shape = treatment),
#           upper = NULL,
#           lower = list(continuous = wrap("points", size = 2)),
#           columns = 1:6,
#           diag = NULL,
#           legend = grab_legend(pca_tmp2)) +
#   ggtitle("Pairwise plots of first 6 PCs",
#           subtitle = "PCA on log2(TPM + 1): 25324 genes (counts > 1)") +
#   theme_minimal()
# ggsave("reports/figures/pca_pairs_log2tmp2.png",
#        plot = pca_pairs_tmp2,
#        device = "png",
#        width = 8,
#        height = 6,
#        units = "in",
#        dpi = 300)


# #############################################################################
# # correlation between samples
# # ############################################################################
# corr_matrix <- cor(tpm)
# anno_coldata <- coldata %>%
#   select(-names, -fcido)
# corrplot(corr_matrix, order = 'hclust',
#          addrect = 2, addCoef.col = 'white',
#          number.cex = 0.7)
# pheatmap(corr_matrix,
#          annotation_col = anno_coldata,
#        #  cutree_cols = 3,
#          width = 8,
#          height = 6,
#          filename = "reports/figures/correlation.png",
#          main = "Correlation between in TMP scores between samples")
# # dev.off()

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
class(de_results) # "DESeqResults"

# # Diagnostic plots
# DESeq2::plotMA(object = dds, alpha = 0.05)
# # most genes should be on the horizontal line as most genes
# # are not differentially expressed
#
# de_results %>% data.frame() %>%
#   ggplot(aes(x = pvalue)) +
#   geom_histogram(bins = 100)
# # there seems to be a few too many high p values

# PCA on normalised counts
# extract normalized counts from the DESeqDataSet object
countsNormalized <- DESeq2::counts(dds, normalized = TRUE)

# all genes

M <- t(countsNormalized)

pca <- M %>%
  prcomp()
# examine amount of variance captured
pcvar <- summary(pca)[["importance"]][2,1:2]  %>% round(4) * 100
# put sample name and pc score in a dataframe
dat <-  data.frame(pca$x, names = rownames(M))

# add the metadata from coldata for annotation of plot
dat <- dat %>% merge(coldata, by = "names")

pca_normc <- ggplot(dat,
                    aes(x = PC1, y = PC2,
                        colour = donor,
                        shape = treatment)) +
  geom_point(size = 3) +
  xlab(paste("PC1: ", pcvar[1],"%")) +
  ylab(paste("PC2: ", pcvar[2],"%")) +
  ggtitle("PCA: normalised counts") +
  theme_classic()

ggsave(pca_normc,
       filename = "reports/figures/pca-normalised-counts-all.png",
       height = 5,
       width = 5)



# select top 500 most variable genes
selectedGenes <- names(sort(apply(countsNormalized, 1, var),
                            decreasing = TRUE)[1:500])

M <- t(countsNormalized[selectedGenes,])

pca <- M %>%
  prcomp()
# examine amount of variance captured
pcvar <- summary(pca)[["importance"]][2,1:2]  %>% round(4) * 100
# put sample name and pc score in a dataframe
dat <-  data.frame(pca$x, names = rownames(M))

# add the metadata from coldata for annotation of plot
dat <- dat %>% merge(coldata, by = "names")

pca_normc <- ggplot(dat,
       aes(x = PC1, y = PC2, colour = donor, shape = treatment)) +
  geom_point(size = 3) +
  xlab(paste("PC1: ", pcvar[1],"%")) +
  ylab(paste("PC2: ", pcvar[2],"%")) +
  ggtitle("PCA: normalised counts top 500 most variable genes") +
  theme_classic()

ggsave(pca_normc,
       filename = "reports/figures/pca-normalised-counts-500.png",
       height = 5,
       width = 5)

#  Relative Log Expression (RLE) plot
# Like the MA plot is useful to see if normalization needed (Gandolfo
# and Speed 2018). Might need normalization more than above because
# of variation in library preparation, experiment conditions etc.
# The RLE plot is a quick diagnostic that can be done to raw or
# normalized counts
# Here RLE plots on the raw counts and normalized counts using
# the EDASeq package (Risso, Schwartz, Sherlock, et al. 2011)
# the boxplots should be centred around the line and be as tightly distrubted as
# possible. we can see how the normalisation has helped.
# there are other packages that can deal with any addition noise.
EDASeq::plotRLE(count_data,
        outline = FALSE,
        col = as.numeric(coldata$treatment),
        main = "Raw Counts")

EDASeq::plotRLE(DESeq2::counts(dds, normalized = TRUE),
        outline = FALSE,
        col = as.numeric(coldata$treatment),
        main = "Normalized Counts")

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

##################################################################
# Differential expression analysis for PCR validated response
# on unstimulated low vs high
##################################################################

# count matrix
count_data <- test_table %>%
  select(S01_P17040_unstim_lo:S34_P16088_stimul_lo) %>%
  as.matrix()
row.names(count_data) <- test_table$gene_id


# define the experimental setup

# define the design formula (different)
design_idoresp <- "~ idoresp"

# create a DESeq dataset object from the count matrix and the coldata
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = as.formula(design_idoresp))
# print dds object to see the contents
print(dds) # dim: 29970 18

# for each gene, count the total number of reads for that gene
# in all samples and remove those that don't have at least 1 read. (two surely?)
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]
print(dds) # dim: 25324 18 4646 removed

# estimate size factors to normalize the counts and dispersion values
# compute GLM model based on the experimental design formula
dds <- DESeq(dds)

# compute the contrast for the idoresp variable where "low"
# is the base
# this means +'ve FC are upreg in high relative to low
de_results = results(dds, contrast = c("idoresp", "high", "low"))

# sort results by increasing p-value
# de_results <- de_results[order(de_results$pvalue),]

print(de_results)
class(de_results) # "DESeqResults"

# not between 70000000000000000000000000000000000000000
# Diagnostic plots
DESeq2::plotMA(object = dds, alpha = 0.05)
# most genes should be on the horizontal line as most genes
# are nor differentially expressed. there are very few diff expressed here

de_results %>% data.frame() %>%
  ggplot(aes(x = pvalue)) +
  geom_histogram(bins = 100)
# there are many high p, suggests no effect

# PCA on normalised counts
# extract normalized counts from the DESeqDataSet object
countsNormalized <- DESeq2::counts(dds, normalized = TRUE)

# select top 500 most variable genes
selectedGenes <- names(sort(apply(countsNormalized, 1, var),
                            decreasing = TRUE)[1:500])

M <- t(countsNormalized[selectedGenes,])

pca <- M %>%
  prcomp()
# examine amount of variance captured
pcvar <- summary(pca)[["importance"]][2,1:2]  %>% round(4) * 100
dat <-  data.frame(pca$x, names = rownames(M))

# add the metadata from coldata for annotation of plot
dat <- dat %>% merge(coldata, by = "names")
ggplot(dat,
       aes(x = PC1, y = PC2, colour = idoresp)) +
  geom_point(size = 3) +
  xlab(paste("PC1: ", pcvar[1],"%")) +
  ylab(paste("PC2: ", pcvar[2],"%")) +
  geom_text_repel(aes(label = donor)) +
  ggtitle("PCA: normalised counts top 500 most variable genes") +
  theme_classic()

#  Relative Log Expression (RLE) plot
EDASeq::plotRLE(count_data,
        outline = FALSE,
        col = as.numeric(coldata$treatment),
        main = "Raw Counts")

EDASeq::plotRLE(DESeq2::counts(dds, normalized = TRUE),
        outline = FALSE,
        col = as.numeric(coldata$treatment),
        main = "Normalized Counts")

# not between 70000000000000000000000000000000000000000

#####below here
###############################################################################
# Volcano plot for Differential expression analysis for IDO response
###############################################################################

de_results_vol <- de_results %>%
  data.frame(gene_id = row.names(de_results)) %>%
  filter(!is.na(padj))
# dim(de_results_vol) 13049
# length(unique(de_results_vol$gene_id))
# AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)
# Annotation
annotations <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                     keys = de_results_vol$gene_id,
                                     columns = c("SYMBOL",
                                                 "ENTREZID",
                                                 "GENENAME",
                                                 "ENSEMBL",
                                                 "GENETYPE"),
                                     keytype = "ENSEMBL")
# dim(annotations) 13113      we have 13113- 13049  = 64 duplicates

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations$ENSEMBL) == FALSE)

# Return only the non-duplicated genes using indices
annotations <- annotations[non_duplicates_idx, ]
# dim(annotations) 13049

# merge annotations in to results
names(annotations)[1] <- "gene_id"
de_results_vol <- de_results_vol %>% merge(annotations, by = "gene_id")

# genetypes
de_results_vol %>% group_by(GENETYPE) %>% summarise(n = length(GENETYPE))
# GENETYPE           n
# ncRNA            911
# 2 protein-coding  9107
# 3 pseudo           244
# 4 unknown            1
# 5 NA              2786

# add a indicator for impt genes - those with FC >=0 and FDR <=0.05
# set FDR values less then 1e-50 to 0
de_results_vol <- de_results_vol %>%
  mutate(sig = padj <= 0.05,
         impt = abs(log2FoldChange) >= 0 &
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

#############################up to here
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
#GO Terms for differentially expressed between
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




























