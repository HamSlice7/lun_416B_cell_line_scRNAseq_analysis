library(scRNAseq)
library(AnnotationHub)
library(scater)


##Loading in the data

#Loading in the SingleCellExperiment object
sce.416b <- LunSpikeInData(which = "416b")

#Convert block to a factor
sce.416b$block <- factor(sce.416b$block)

#AnnotationDb object - Gene and protein annotations for Mus musculus based on Ensembl version 97
ens.mm.v97 <- AnnotationHub()[["AH73905"]]

#Adding some meta data on the genes
rowData(sce.416b)$ENSMEBL <- rownames(sce.416b)
rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v97, keys = rownames(sce.416b), keytype = "GENEID", column = "SYMBOL")
rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v97, keys = rownames(sce.416b), keytype = "GENEID", column = "SEQNAME")

#Replace row names with the gene symbol if unique but if not unique or missing use ensembl id
rownames(sce.416b) <- uniquifyFeatureNames(rowData(sce.416b)$ENSMEBL, rowData(sce.416b)$SYMBOL)


##Quality Control