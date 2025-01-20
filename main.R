library(scRNAseq)
library(AnnotationHub)
library(scater)
library(scran)

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

#Saving unfiltered object
unfiltered <- sce.416b

#Filtering based on MT genes
mito <- which(rowData(sce.416b)$SEQNAME=="MT")
stats <- perCellQCMetrics(sce.416b, subsets=list(Mt=mito))
qc <- quickPerCellQC(stats, percent_subsets = c("subsets_Mt_percent", "altexps_ERCC_percent"), batch=sce.416b$block)
sce.416b <- sce.416b[,!qc$discard]

#Plotting the cells that are to be discarded or not and QC metrics
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$block <- factor(unfiltered$block)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
  plotColData(unfiltered, x="block", y="sum", 
              colour_by="discard") + scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, x="block", y="detected", 
              colour_by="discard") + scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, x="block", y="subsets_Mt_percent", 
              colour_by="discard") + ggtitle("Mito percent"),
  plotColData(unfiltered, x="block", y="altexps_ERCC_percent", 
              colour_by="discard") + ggtitle("ERCC percent"),
  nrow=2,
  ncol=2
)

#plotting MT percent vs total count and spike-in percent
gridExtra::grid.arrange(
  plotColData(unfiltered, x="sum", y="subsets_Mt_percent", 
              colour_by="discard") + scale_x_log10(),
  plotColData(unfiltered, x="altexps_ERCC_percent", y="subsets_Mt_percent",
              colour_by="discard"),
  ncol=2
)

#Examine the number of cells for each reason
colSums(as.matrix(qc))

##Normalization

#calculating size factor for each cell
sce.416b <- computeSumFactors(sce.416b)

#Compute log-transformed normalized expression values
sce.416b <- logNormCounts(sce.416b)

summary(sizeFactors(sce.416b))