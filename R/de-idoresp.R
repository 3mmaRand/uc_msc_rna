# Differential expression between high and low IDO responders

library("DESeq2")

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
de_results = results(dds, contrast = c("idoresp", "high", "low"))

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
