# this script use the counts matrix

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
nrow(reads_count) # 61541

# calculate the sum of counts across samples
# and filter to keep only those with a sum at least 10
reads_count <- reads_count %>%
  rowwise() %>%
  mutate(sum = sum(c_across(S01_P17040_unstim:S34_P16088_stimul),
                   na.rm = T)) %>% ungroup() %>%
  filter(sum >= 10)
nrow(reads_count) # 18169
# has taken it down to 18169 rows
# 61541 - 18169 = 43372 have counts < 10

# VST transform the raw counts for putting in PCA

# make the count matrix
# select just the counts
cols <- coldata$names
countdata <- reads_count  %>%
  select(all_of(cols)) %>%
  round(1)
row.names(countdata) <- reads_count$gene_id
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
ggsave("reports/figures/sample_distance_heatmap-counts.png",
       plot = p1,
       device = "png",
       width = 6,
       height = 7,
       units = "in",
       dpi = 300)
