library(scRNAseq)
library(AnnotationHub)
library(scater)
library(scran)
library(limma)

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

#plotting deconvolution factor vs library size factor
plot(librarySizeFactors(sce.416b), sizeFactors(sce.416b), pch = 16, 
     xlab = "Library size factors", ylab = "Deconvolution factors", 
     col = c("black", "red")[grepl("induced", sce.416b$phenotype)+1],
     log="xy")

##Variance modelling. - block on the plate of origin to minimize plate effects
dec.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC", block=sce.416b$block)
#Take the top 10% of genes with the largest biological component (variance)
chosen.hvgs <- getTopHVGs(dec.416b, prop = 0.1)


par(mfrow=c(1,2))
blocked.stats <- dec.416b$per.block
for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  plot(current$mean, current$total, main=i, pch = 16, cex=0.5,
       xlab = "Mean of log-expression", ylab = "Variance of log-expression")
  curfit <- metadata(current)
  points(curfit$mean, curfit$var, col = "red", pch=16)
  curve(curfit$trend(x), col = "dodgerblue", add=TRUE, lwd=2)
}


##Batch correction
#Composition of cells is expected to be the same across the two plates, hence the use of removeBatchEffect() rather than more complex methods
assay(sce.416b, "corrected") <- removeBatchEffect(logcounts(sce.416b), design = model.matrix(~sce.416b$phenotype), batch = sce.416b$block)


##Dimensionality reduction

#We do not expect a great deal of heterogeneity in this dataset, so we only request 10 PCs. 
#We use an exact SVD to avoid warnings from irlba about handling small datasets.
sce.416b <- runPCA(sce.416b, ncomponents = 10, subset_row = chosen.hvgs,
                   exprs_values="corrected", BSPARAM=BiocSingular::ExactParam())

set.seed(1010)
sce.416b <- runTSNE(sce.416b, dimred="PCA", perplexity = 10)


##Clustering
