library("tidyverse")



# files

# the count data for all IDs including unexpressed
reads_count_f   <- "data-raw/DE_Results/readsCount.tab"
# the IDs which are not expressed
zero_reads_f    <- "data-raw/DE_Results/allZeroreads.tab"
# differential expression analysis results for counts above 0
# where the comparison is between stimul and unstimul
# note: if is meaningful to examine response of donors
test_table_f    <- "data-raw/DE_Results/testTable.tab"
top_100_table_f <- "data-raw/DE_Results/top.100.Table.tab"
coldata_f    <- "background/coldata"

reads_count   <- read_delim(reads_count_f)
test_table    <- read_delim(test_table_f)
top_100_table <- read_delim(top_100_table_f)


# rename samples in test_table using "test-table-names.txt"
# so that the sample encodes sample number, patient and treatment
# S##_p#####_treatment so that S01_P17040_unstim is patient 17/040 unstimulated
# aka sample 1 and S26_P15142_stimul is patient 15/142 stimulated aka sample 26
names(reads_count) <- readLines("background/reads-count-names.txt")
names(test_table) <- readLines("background/test-table-names.txt")
names(top_100_table) <- readLines("background/test-table-names.txt")
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

levels(coldata$idoresp)
# [1] "high" "low"

# make unstim the base level and low the base level
coldata$treatment <- relevel(coldata$treatment, "unstim")
coldata$idoresp <- relevel(coldata$idoresp, "low")


