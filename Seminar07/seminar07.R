library(limma)
library(edgeR)

dat <- read.table("../data/bottomly_count_table.tsv", header = TRUE, row.names = 1)
des <- read.table("../data/bottomly_phenodata.tsv", header = TRUE, row.names = 1)
str(dat)

identical(rownames(des), colnames(dat))

with(des, table(strain))

group <- factor(c(rep("1", 10), rep("2", 11)))
group

dge.glm <- DGEList(counts = dat, group = group)
str(dge.glm)

names(dge.glm)
dge.glm[["samples"]]
nrow(dge.glm[[1]])
ncol(dge.glm[[1]])

design <- model.matrix(~group)
design

#common line
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)

#fits the trend
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp)

#squeezes the trend 
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)

plotBCV(dge.glm.tag.disp)

fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

tt.glm <- topTags(lrt, n = Inf)
class(tt.glm)

nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])

eHits <- tt.glm$table[tt.glm$table$FDR < 1e-4, ]

interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]

summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))

# plotting the tagwise log fold changes against log-cpm
tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")


#Mini exercise

newDat <- dat[rowSums(dat)>0, ]

selRows  <- apply(newDat, 1, function(row) all(row[1:10]!=0) | all(row[11:21]!=0))

str(selRows)

newestDat <- newDat[selRows, ]

dge.glm <- DGEList(counts = newestDat, group = group)
str(dge.glm)

names(dge.glm)
dge.glm[["samples"]]
nrow(dge.glm[[1]])
ncol(dge.glm[[1]])

design <- model.matrix(~group)
design

#common line
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)

#fits the trend
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp)

#squeezes the trend 
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)

plotBCV(dge.glm.tag.disp)

fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

tt.glm <- topTags(lrt, n = Inf)
class(tt.glm)

nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])

interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]

summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))

# plotting the tagwise log fold changes against log-cpm
tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")

################################################################ DESeq

library(DESeq)

# reading in the same count table data and grouping information
deSeqDat <- newCountDataSet(dat, group)
head(counts(deSeqDat))

deSeqDat <- estimateSizeFactors(deSeqDat)
sizeFactors(deSeqDat)

deSeqDat <- estimateDispersions(deSeqDat)
# plotting the estimated dispersions against the mean normalized counts
plotDispEsts(deSeqDat)

## this takes a minute or so for JB
results <- nbinomTest(deSeqDat, levels(group)[1], levels(group)[2])
str(results)


dHits <- results[results$padj < 1e-4, ]

plotMA(results)

################################################################ Voom

library(limma)
norm.factor <- calcNormFactors(dat)
dat.voomed <- voom(dat, design, plot = TRUE, lib.size = colSums(dat) * norm.factor)

dat.voomed

fit <- lmFit(dat.voomed, design)
fit <- eBayes(fit)
topTable(fit)

############################################################### Take home
library(gplots)

#edgeR (calculated above)
eHits <- eHits

#DESeq (calculated above)
dHits <- na.exclude(dHits)
rownames(dHits) <- dHits$id

#voom + limma
vHits <- topTable(fit, p.value=1e-04, number=nrow(dat), coef="group2")

input = list(rownames(vHits), rownames(dHits), rownames(eHits))

str(input)

venn(input)
