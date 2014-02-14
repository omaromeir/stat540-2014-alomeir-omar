library(lattice)
library(limma)
prDat <- read.table("../data/GSE4051_data.tsv")
prDes <- readRDS("../data/GSE4051_design.rds")
source("../Seminar05/seminar05.R")

m <- 1000
n <- 3
x <- matrix(rnorm(m * n), nrow = m)

obsVars <- apply(x, 1, var)
summary(obsVars)
mean(obsVars < 1/3)
densityplot(~obsVars, n=200)

#take home exercise: Make the above simulation more realistic with two (or more) groups, different data-generating means and group differences, different data-generating gene-wise variances, etc.

#2 different groups here
m <- 1000
n <- 3
g <- 2

x <- matrix(rnorm(m * n, mean = rep(c(1, 2), each=(nrow(x)/g)), sd = rep(c(sqrt(1), sqrt(2)), each=(nrow(x)/g))), nrow = m)
gNames  <- paste0("Group ", seq(1, g))
gMat <- cbind(rep(gNames , each=(nrow(x)/g)), x)
colnames(gMat) <- c("Group", "dat1", "dat2", "dat3")

obsVars <- apply(gMat[,grep("dat", colnames(gMat))], 1, var)
summary(obsVars)
mean(obsVars < 1/3)
densityplot(~obsVars, n=200)

#Fit a linear model

wtDes <- subset(prDes, gType == "wt")
str(wtDes)
wtDat <- subset(prDat, select = prDes$gType == "wt")
str(wtDat, max.level = 0)

wtDesMat <- model.matrix(~devStage, wtDes)
str(wtDesMat)

wtFit <- lmFit(wtDat, wtDesMat)
wtEbFit <- eBayes(wtFit)

topTable(wtEbFit)
colnames(coef(wtEbFit))

(dsHits <- topTable(wtEbFit, coef = grep("devStage", colnames(coef(wtEbFit)))))

stripplot(gExp ~ devStage | gene, subset= gType=="wt", prepareData(rownames(dsHits)[c(3, 6, 9)]), jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE)

#Does it look plausible to you that -- using only wild type data -- these probes show the most compelling evidence for expression change over development? yep

#Exercise: lm() on one or all 3 of these probes and check if the F stats and p-values are similar. Don't expect exact equality because you must remember that limma has moderated the estimated error variance.
mFit <- lm(formula = gExp ~ devStage, data=prepareData(rownames(dsHits)[3]), subset = gType=="wt")
summary(mFit)
#The numbers are similar, though not exactly the same

#Be the boss of top table
newHits <- topTable(wtEbFit, p.value=1e-05, number=nrow(wtDat), coef = grep("devStage", colnames(coef(wtEbFit))))
str(newHits)
newHits[63, c("F", "adj.P.Val", 'devStageP6')]

p2Hits <- topTable(wtEbFit, number=nrow(wtDat), coef = "devStageP2", sort.by="none")
str(p2Hits)

p10Hits <- topTable(wtEbFit, number=nrow(wtDat), coef = "devStageP10", sort.by="none")
str(p10Hits)

xyplot(p10Hits$t ~ p2Hits$t, scales= list(limit = c(-20, 20)), asp = 1, panel=function(...) {panel.smoothScatter(...); panel.abline(0,1,col=3)})

densityplot(~ p10Hits$adj.P.Val + p2Hits$adj.P.Val, auto.key=TRUE, type = 'p')

addmargins(table(p2Hits$adj.P.Val < 1e-03, p10Hits$adj.P.Val < 1e-03, dnn=c("p2", "p10")))

p10HitsBY <- topTable(wtEbFit, number=nrow(wtDat), coef = "devStageP10", adjust.method="BY", sort.by="none")
pVals <- data.frame(Raw = p10Hits$P.Value, BH = p10Hits$adj.P.Val, BY = p10HitsBY$adj.P.Val)
pairs(pVals)

#Relationship is linear between raw and BH
tail(pVals)

#Inference for some contrasts
colnames(wtDesMat)
(cont.matrix <- makeContrasts(P10VsP6 = devStageP10 - devStageP6, fourweeksVsP10 = devStage4_weeks - devStageP10, levels = wtDesMat))
wtFitCont <- contrasts.fit(wtFit, cont.matrix)
wtEbFitCont <- eBayes(wtFitCont)
latestHits <- topTable(wtEbFitCont)
stripplot(gExp ~ devStage | gene, subset= gType=="wt", asp=1, prepareData(rownames(latestHits)[1:4]), jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE)

cutoff <- 1e-04
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)

(hits1 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] < 0)])
stripplot(gExp ~ devStage | gene, subset= gType=="wt", asp=1, prepareData(hits1), jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE)
(hits2 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] < 0)])
stripplot(gExp ~ devStage | gene, subset= gType=="wt", asp=1, prepareData(hits2), jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE)
intersect(hits1, hits2)
(hits3 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] > 0)])
stripplot(gExp ~ devStage | gene, subset= gType=="wt", asp=1, prepareData(hits3[1:8]), jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE)
intersect(hits1, hits3)
intersect(hits2, hits3)

cutoff <- 0.01
nHits <- 8
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)

hits1 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] < 0)]
hits2 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] < 0)]
hits3 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] > 0)]
hits4 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] > 0)]
vennDiagram(wtResCont)
hits6 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] > 0 & wtResCont[, "fourweeksVsP10"] < 0)]
stripplot(gExp ~ devStage | gene, subset= gType=="wt", asp=1, prepareData(hits6[1:8]), jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE)

#take home exercise
(cont.matrix <- makeContrasts(P2VsI = devStageP2 - Intercept, P6VsP2 = devStageP6 - devStageP2, P10VsP6 = devStageP10 - devStageP6, fourweeksVsP10 = devStage4_weeks - devStageP10, levels = wtDesMat))
wtFitCont <- contrasts.fit(wtFit, cont.matrix)
wtEbFitCont <- eBayes(wtFitCont)
toptable(wtEbFitCont)
cutoff <- 1e-04
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)

greatestHits <- rownames(prDat)[which(wtResCont[, "P2VsI"] < 0 & wtResCont[, "P6VsP2"] > 0 & wtResCont[, "P10VsP6"] == 0 & wtResCont[, "fourweeksVsP10"] == 0)]
stripplot(gExp ~ devStage | gene, subset= gType=="wt", asp=1, prepareData(greatestHits), jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE)
