# DESeq based SVA workflow, modelled on https://biodatascience.github
# By Mesude Bicak and Benjamin Readhead
######################################################################

library(DESeq2)
library(variancePartition)
library('doParallel')
library(sva)
library(pbapply)
library(plyr)

options(stringsAsFactors = FALSE)

allCountsFile <- "borchelt.final"

#Read in read counts
data <- read.delim(allCountsFile)
rowNames <- as.character((data$GeneID))
data <- data.matrix(data)
row.names(data) <- rowNames
data <- data[,-1]
dim(data)

#Assign sample groups
coldata <- data.frame(samples = colnames(data), group = factor(c(rep("WT",5), rep("BNN",5), rep("BNY",4), rep("BYN",5), rep("BYY",5))))

dds <- DESeqDataSetFromMatrix(data, coldata, ~group)
dds <- estimateSizeFactors(dds)
norm.cts <- counts(dds, normalized=TRUE)

#Filter samples with zero counts
norm.cts <- norm.cts[rowSums(norm.cts) > 0,]
norm.cts <- norm.cts[which(rowSums(norm.cts == 0) <= 3),] 

#Normalised data
write.table(norm.cts, file="normalised.borchelt.final", row.names=T, sep="\t", quote=F)

#Surrogate Variable Analysis (SVA) will be used to estimate surrogate variables to adjust for unknown sources of noise
mm <- model.matrix(~ group, colData(dds))
mm0 <- model.matrix(~ 1, colData(dds))
fit <- svaseq(norm.cts, mod=mm, mod0=mm0)#leave n.sv blank, and function will automatically estimate
#fit$n.sv will give # of surrogate variables

colnames(fit$sv) <- paste("SV", 1:fit$n.sv, sep = "")
colData(dds) <- cbind(colData(dds), fit$sv)

design(dds) <- as.formula(paste("~ BCI + Exercise + BCI:Exercise",  paste(paste("SV", 1:fit$n.sv, sep = ""),collapse = " + "), sep = " + "))

dds <- DESeq(dds)

#Comparisons to do:
contrastList <- list(
  BorcheltControl_vs_WT = c("BNN","WT"),
  BorcheltBCI_vs_WT = c("BNY","WT"),
  BorcheltExercise_vs_WT = c("BYN","WT"),
  BorcheltBCIExercise_vs_WT = c("BYY","WT"),
  BorcheltBCIExercise_vs_BorcheltControl = c("BYY","BNN"),
  BorcheltBCI_vs_BorcheltControl = c("BNY","BNN"),
  BorcheltExercise_vs_BorcheltControl = c("BYN","BNN"),
  BorcheltBCI_vs_BorcheltExercise = c("BNY","BYN"),
  BorcheltBCIExercise_vs_BorcheltExercise = c("BYY","BYN"),
  BorcheltBCIExercise_vs_BorcheltBCI = c("BYY","BNY")
)

#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering
keep <- which(rowSums(counts(dds) == 0) <= 3) 
dds <- dds[keep,]
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


DEList <-
  pblapply(1:length(contrastList), function(x) {
    res <- results(dds, contrast=c("group",contrastList[[x]][1],contrastList[[x]][2]))
    res <- res[complete.cases(res),]
    dds_DE <- data.frame(res[order(res$pvalue),])
    dds_DE <- data.frame(symbol = rownames(dds_DE), dds_DE, row.names = NULL)
    dds_DE
})

names(DEList) <- names(contrastList)

DE_df <- ldply(DEList, rbind, .id = "Comparison")
table(DE_df$Comparison[which(DE_df$padj <= 0.05)])

write.table(table(DE_df$Comparison[which(DE_df$padj <= 0.05)]), file="post-sva-borchelt2019-TableOfDEGs-allcomparisons", row.names=T, sep="\t", quote=F)
write.table(DE_df, file="post-sva-borchelt2019-toptable-allcomparisons", row.names=T, sep="\t", quote=F)


