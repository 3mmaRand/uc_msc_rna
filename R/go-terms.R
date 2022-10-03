
library("topGO")
library("org.Hs.eg.db")

# mapping gene_id with go terms
gene2GOdf <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = test_table$gene_id,
                                   columns = c("GO", "ONTOLOGY", "GENENAME"),
                                   keytype = "ENSEMBL")[1:5]
# 'select()' returned 1:many mapping between keys and columns

length(unique(test_table$gene_id))
length(unique(test_table$gene_name))


# turn in to a list of named char vectors required by topGO
gene2GO <- split(gene2GOdf$GO, gene2GOdf$ENSEMBL)
