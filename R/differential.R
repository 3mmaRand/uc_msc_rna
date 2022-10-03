library("tidyverse") # data import, wrangling
library("DESeq2") # differential expression analysis
library("pheatmap")
library("RColorBrewer")

# create a count matrix
# create a metadata df
# check the unstim is level 1 of the factor
# use to make a DESeqDataSet

# files
test_table_f <- "data-raw/DE_Results/testTable.tab"
coldata_f    <- "background/coldata"

# COLUMN (META) DATA

# import column data
# information on each of the samples (the columns of the count matrix)
coldata <- read.table(coldata_f,
                      stringsAsFactors = T,
                      header = TRUE)
row.names(coldata) <- coldata$names
coldata
# check the order of the levels
levels(coldata$treatment)
# [1] "stimul" "unstim"

# make unstim the base level
coldata$treatment <- relevel(coldata$treatment, "unstim")


# COUNT DATA
# import count data for non-zero counts
test_table    <- read_delim(test_table_f)

# rename samples in test_table using "test-table-names.txt"
# so that the sample encodes sample number, patient and treatment
# S##_p#####_treatment so that S01_P17040_unstim is patient 17/040 unstimulated
# aka sample 1 and S26_P15142_stimul is patient 15/142 stimulated aka sample 26
names(test_table) <- readLines("background/test-table-names.txt")

# make the count matrix
# select just the counts
cols <- coldata$names
countdata <- test_table  %>%
  select(all_of(cols)) %>%
  as.data.frame() %>% # because tibbles don't have rownames but the example does
  round(1)

row.names(countdata) <- test_table$gene_id

head(countdata, 3)

# check that the order of the columns in countdata is the same as the order
# of the rows in coldata
names(countdata) == coldata$names

# build a DESeqDataSet from a count matrix and table of sample info


ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ donor + treatment)


head(assay(ddsMat))

# filter rows that have no or nearly no information about the amount of gene
# expression (i think no rows will be removed as thet empty rows , or only a single count across all samples, from readCounts.tab
# zeros are are already removed from testTable.tab)
nrow(ddsMat)
keep <- rowSums(counts(ddsMat)) > 1
table(keep)
ddsMat <- ddsMat[keep,]
nrow(ddsMat)
# has taken it down to 25324 rows - 29970 - 25324 = 4646 havw only one count

# variance stabilising transformation and rlog
# methods like PCA work best if data homoskedastic. with counts, the
# variance increases with the mean
# if we apply PCA to counts or normalised counts (for differences in
# sequencing depth) then the PCA will be dominated by the highest counts
# One solution is to take the log of the normalised count plus a
# pseudo count of 1. But this means the genes with the lowest counts
# now contribute a lot of noise because the log of small counts
# inflates their variance. this means low count genes with low signal-to-noise
# ration will over contribute to sample-sample differences
# DESeq2 has two transformation solutions to stabilise the variance across
# the mean:
#      1. *variance stabilizing transformation* (VST) for negative binomial
#      data with a dispersion-mean trend [@Anders2010Differential],
#      implemented in the *vst* function
#      2. the *regularized-logarithm transformation* or *rlog*
#      [@Love2014Moderated]
#

# VST or rlog transformed data can be used to compute distances between
# samples, making PCA plots or as input to downstream methods that need
# homoskedastic data
#
# VST - quick and less senstive to high count outliers
# rlog - good on small datasets and when sequencing depth varies
# by an order of magnitude

# use VST for medium to large datasets
# can use both and compare the meanSdPlot

# note the two transformations the DESeq2 will do are for tasks OTHER
# than differential expression.
# for differential expression, use the raw counts - DESeq() takes account
# of the dependence of the variance on the mean

# note, because these are made from the count matrix not a summarisedExperiemnt,
# they are not corrected for library size
# Sequencing depth correction is done automatically for the *vst* and *rlog*
# ONLY if done with the summarizedExperiment

# VST transformation
dds_vst <- vst(ddsMat, blind = FALSE)
head(assay(dds_vst))

# rlog transformation
dds_rlog <- vst(ddsMat, blind = FALSE)
head(assay(dds_rlog))

# we need to first estimate *size factors* to
# account for sequencing depth, and then specify `normalized=TRUE`.
ddsMat <- estimateSizeFactors(ddsMat)

# GLM-PCA, a generalization of PCA to exponential family likelihoods.
#  GLM-PCA operates on raw counts, avoiding the pitfalls of normalization


# PCA
# on counts
library("glmpca")
gpca <- glmpca(counts(ddsMat), L = 2)
gpca_dat <- gpca$factors
gpca_dat$treatment <- ddsMat$treatment
gpca_dat$donor <- ddsMat$donor

ggplot(gpca_dat,
       aes(x = dim1, y = dim2, color = donor)) +
  geom_point(size = 3, aes(shape = treatment)) +
  geom_line() +
  coord_fixed() +
  ggtitle("glmpca - Generalized PCA")
# there is definitely an effect of patient as well as treatment

# PCA
# on transformed counts
pcaData <- plotPCA(dds_vst,
                   intgroup = c( "treatment", "donor"),
                   returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))


ggplot(pcaData, aes(x = PC1, y = PC2, color = donor, shape = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")



# Heatmap of sample-to-sample distances using the variance
# stabilizing transformed values.
# sample distances
sampleDists <- dist(t(assay(dds_vst)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists)
rownames(sampleDistMatrix) <- paste(dds_vst$treatment,
                                    dds_vst$donor, sep = " - " )

colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)




pcaData <- plotPCA(dds_vst,
                   intgroup = c( "treatment", "donor"),
                   returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))


ggplot(pcaData, aes(x = PC1, y = PC2, color = donor, shape = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")





# differential expression

ddsMat <- DESeq(ddsMat)

results <- results(ddsMat, alpha = 0.05)
# automatically does
# results(ddsMat, contrast=c("treatment","stimul","unstim"))

mcols(results, use.names = TRUE)
summary(results)
test <- data.frame(results, gene_id = rownames(results))

write_delim(test, "test.tab")







