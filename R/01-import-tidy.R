
# imports the count data from testTable.tab
# renames the samples for ease and consistency
# filters out the genes with fewer than 2 counts.


library("tidyverse")

# filtering reads_count > 1 gives same results as test_table
# using test_table because it has the info on code length

test_table_f  <- "data-raw/DE_Results/testTable.tab"
test_table    <- read_delim(test_table_f)

# rename samples in test_table using "reads-count-names.txt"
# so that the sample encodes sample number, patient and treatment
# S##_p#####_treatment so that S01_P17040_unstim is patient 17/040 unstimulated
# aka sample 1 and S26_P15142_stimul is patient 15/142 stimulated aka sample 26
names(test_table) <- readLines("background/test-table-names.txt")


# calculate the sum of counts across samples
# and filter to keep only those with a sum greater than 1
# dim(test_table) 29970
test_table <- test_table %>%
  rowwise() %>%
  mutate(sum = sum(c_across(S01_P17040_unstim_lo:S34_P16088_stimul_lo),
                   na.rm = T)) %>% ungroup() %>%
  filter(sum > 1)
# dim(test_table)  25324  (29970 - 25324 = 4646 genes removed)

# COLUMN (META) DATA
# import column data
# information on each of the samples (the columns of the count matrix)
coldata_f    <- "background/coldata"
coldata <- read.table(coldata_f,
                      stringsAsFactors = T,
                      header = TRUE)
row.names(coldata) <- coldata$names

coldata <- coldata %>% mutate("ATP_norm" = ATP/ATP_PBMC)

# make unstim the base level and low the base level
coldata$treatment <- relevel(coldata$treatment, "unstim")
coldata$idoresp <- relevel(coldata$idoresp, "low")
coldata$idotrans <- relevel(coldata$idotrans, "low")
