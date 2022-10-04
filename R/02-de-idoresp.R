##################################################################
# Differential expression between high and low IDO responders
##################################################################
library("DESeq2")


#################################################################
# Using both the stimulated and unstimulated data
#################################################################

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
# print(dds) # dim: 29970 18

# for each gene, count the total number of reads for that gene
# in all samples and remove those that don't have at least 1 read. (two surely?)
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]
print(dds) # dim: 25324 18 4646 removed

# estimate size factors to normalize the counts and dispersion values
# compute GLM model based on the experimental design formula
dds <- DESeq(dds)
# fitting model and testing
# -- replacing outliers and refitting for 35 genes

# compute the contrast for the idoresp variable where "low"
# is the base
# this means +'ve FC are upreg in high relative to low
de_results  <- results(dds, contrast = c("idoresp", "high", "low"))

###############################################################################
# Volcano data for Differential expression analysis for IDO response
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


# add a indicator for impt genes - those with FC > 0 and FDR <=0.05
# set FDR values less then 1e-50 to 0
de_results_vol <- de_results_vol %>%
  mutate(sig = padj <= 0.05,
         impt = abs(log2FoldChange) > 0 &
           padj <= 0.05,
         fdr = case_when(padj < 1e-50 ~ 0,
                         padj >= 1e-50 ~ padj))
write_csv(de_results_vol,
          file = "data-processed/de_ido_resp_results.csv")

de_results_vol %>%
  group_by(sig) %>%
  summarise(n = length(GENETYPE)) %>%
  kable(row.names = FALSE) %>%
  kable_styling(font_size = 12)

# it is not worth doing a volcano plot for these results.



#################################################################
# Using just the stimulated data
#################################################################

# count matrix
count_data_stim <- test_table %>%
  select(S02_P17040_stimul_lo,
         S04_P17041_stimul_lo,
         S06_P17083_stimul_hi,
         S14_P14086_stimul_hi,
         S16_P12316_stimul_hi,
         S12_P17055_stimul_hi,
         S32_P16143_stimul_hi,
         S26_P15142_stimul_lo,
         S34_P16088_stimul_lo) %>%
  as.matrix()
row.names(count_data_stim) <- test_table$gene_id


# define the experimental setup

# define the design formula (different)
design_idoresp <- "~ idoresp"

# filter the coldata so there's only stimul

coldata_stim <- coldata %>%
  filter(treatment == "stimul") %>% droplevels()

# create a DESeq dataset object from the count matrix and the coldata_stim
dds_stim <- DESeqDataSetFromMatrix(countData = count_data_stim,
                              colData = coldata_stim,
                              design = as.formula(design_idoresp))
# print dds_stim object to see the contents
# print(dds_stim) # dim: 25324

# for each gene, count the total number of reads for that gene
# in all samples and remove those that don't have at least 2 reads
dds_stim <- dds_stim[ rowSums(DESeq2::counts(dds_stim)) > 1, ]
print(dds_stim) # dim: 21755  3569 removed

# estimate size factors to normalize the counts and dispersion values
# compute GLM model based on the experimental design formula
dds_stim <- DESeq(dds_stim)

# compute the contrast for the idoresp variable where "low"
# is the base
# this means +'ve FC are upreg in high relative to low
de_results_stim  <- results(dds_stim, contrast = c("idoresp", "high", "low"))

###############################################################################
# Volcano data for Differential expression analysis for IDO response: stim only
###############################################################################

de_results_stim_vol <- de_results_stim %>%
  data.frame(gene_id = row.names(de_results_stim)) %>%
  filter(!is.na(padj))
# dim(de_results_stim_vol) 1931
# length(unique(de_results_stim_vol$gene_id)) 1931
# AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)
# Annotation
annotations_stim <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                     keys = de_results_stim_vol$gene_id,
                                     columns = c("SYMBOL",
                                                 "ENTREZID",
                                                 "GENENAME",
                                                 "ENSEMBL",
                                                 "GENETYPE"),
                                     keytype = "ENSEMBL")
# dim(annotations_stim) 1941    we have 10 duplicates

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_stim$ENSEMBL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_stim <- annotations_stim[non_duplicates_idx, ]
# dim(annotations_stim) 1931

# merge annotations in to results
names(annotations_stim)[1] <- "gene_id"
de_results_stim_vol <- de_results_stim_vol %>% merge(annotations_stim, by = "gene_id")


# add a indicator for impt genes - those with FC > 0 and FDR <=0.05
# set FDR values less then 1e-50 to 0
de_results_stim_vol <- de_results_stim_vol %>%
  mutate(sig = padj <= 0.05,
         impt = abs(log2FoldChange) > 0 &
           padj <= 0.05,
         fdr = case_when(padj < 1e-50 ~ 0,
                         padj >= 1e-50 ~ padj))
write_csv(de_results_stim_vol,
          file = "data-processed/de_ido_resp_stim_results.csv")

de_results_stim_vol %>%
  group_by(sig) %>%
  summarise(n = length(GENETYPE)) %>%
  kable(row.names = FALSE) %>%
  kable_styling(font_size = 12)

de_results_stim_vol %>%
  filter(sig) %>%
  select(SYMBOL, gene_id, log2FoldChange, padj, GENENAME) %>%
  kable(row.names = FALSE, digits = 4) %>%
  kable_styling(font_size = 12)

# it is not worth doing a volcano plot for these results.


#################################################################
# Using just the unstim data
#################################################################

# count matrix
count_data_unstim <- test_table %>%
  select(S01_P17040_unstim_lo,
         S03_P17041_unstim_lo,
         S05_P17083_unstim_hi,
         S07_P17055_unstim_hi,
         S13_P14086_unstim_hi,
         S15_P12316_unstim_hi,
         S31_P16143_unstim_hi,
         S25_P15142_unstim_lo,
         S33_P16088_unstim_lo) %>%
  as.matrix()
row.names(count_data_unstim) <- test_table$gene_id


# define the experimental setup

# define the design formula (different)
design_idoresp <- "~ idoresp"

# filter the coldata so there's only stimul

coldata_unstim <- coldata %>%
  filter(treatment == "unstim") %>% droplevels()

# create a DESeq dataset object from the count matrix and the coldata_unstim
dds_unstim <- DESeqDataSetFromMatrix(countData = count_data_unstim,
                                     colData = coldata_unstim,
                                     design = as.formula(design_idoresp))
# print dds_unstim object to see the contents
# print(dds_unstim) # dim: 25324

# for each gene, count the total number of reads for that gene
# in all samples and remove those that don't have at least 2 reads
dds_unstim <- dds_unstim[ rowSums(DESeq2::counts(dds_unstim)) > 1, ]
print(dds_unstim) # dim: 21755  3569 removed

# estimate size factors to normalize the counts and dispersion values
# compute GLM model based on the experimental design formula
dds_unstim <- DESeq(dds_unstim)

# compute the contrast for the idoresp variable where "low"
# is the base
# this means +'ve FC are upreg in high relative to low
de_results_unstim  <- results(dds_unstim, contrast = c("idoresp", "high", "low"))

###############################################################################
# Volcano data for Differential expression analysis for IDO response: unstim only
###############################################################################

de_results_unstim_vol <- de_results_unstim %>%
  data.frame(gene_id = row.names(de_results_unstim)) %>%
  filter(!is.na(padj))
# dim(de_results_unstim_vol) 22505
# length(unique(de_results_unstim_vol$gene_id)) 22505
# AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)
# Annotation
annotations_unstim <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                            keys = de_results_unstim_vol$gene_id,
                                            columns = c("SYMBOL",
                                                        "ENTREZID",
                                                        "GENENAME",
                                                        "ENSEMBL",
                                                        "GENETYPE"),
                                            keytype = "ENSEMBL")
# dim(annotations_unstim) 1941    we have 122 duplicates

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_unstim$ENSEMBL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_unstim <- annotations_unstim[non_duplicates_idx, ]
# dim(annotations_unstim) 1931

# merge annotations in to results
names(annotations_unstim)[1] <- "gene_id"
de_results_unstim_vol <- de_results_unstim_vol %>% merge(annotations_unstim, by = "gene_id")


# add a indicator for impt genes - those with FC > 0 and FDR <=0.05
# set FDR values less then 1e-50 to 0
de_results_unstim_vol <- de_results_unstim_vol %>%
  mutate(sig = padj <= 0.05,
         impt = abs(log2FoldChange) > 0 &
           padj <= 0.05,
         fdr = case_when(padj < 1e-50 ~ 0,
                         padj >= 1e-50 ~ padj))
write_csv(de_results_unstim_vol,
          file = "data-processed/de_ido_resp_unstim_results.csv")

de_results_unstim_vol %>%
  group_by(sig) %>%
  summarise(n = length(GENETYPE)) %>%
  kable(row.names = FALSE) %>%
  kable_styling(font_size = 12)
# there are no sig results
de_results_unstim_vol %>%
  filter(sig) %>%
  select(SYMBOL, gene_id, log2FoldChange, padj, GENENAME) %>%
  kable(row.names = FALSE, digits = 4) %>%
  kable_styling(font_size = 12)

# it is not worth doing a volcano plot for these results.
