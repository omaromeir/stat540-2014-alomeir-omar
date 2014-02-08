library(lattice)
kDat <- readRDS("GSE4051_MINI.rds")
str(kDat)
table(kDat$devStage)
xyplot(eggBomb + poisonFang ~ crabHammer, kDat, groups=gType, auto.key = TRUE, outer=TRUE, grid=TRUE)

nDat <-with(kDat, data.frame(sidChar, sidNum, devStage, gType, crabHammer, probeset = factor(rep(c("eggBomb", "poisonFang"), each = nrow(kDat))), geneExp = c(eggBomb, poisonFang)))
str(nDat)

xyplot(geneExp ~ crabHammer | probeset, nDat, grid = TRUE, groups = devStage, auto.key = TRUE)

oDat <- with(kDat, data.frame(sidChar, sidNum, devStage, gType, probeset = factor(rep(c("crabHammer", "eggBomb", "poisonFang"), each = nrow(kDat))), geneExp = c(crabHammer, eggBomb, poisonFang)))
str(oDat)

stripplot(~ geneExp | probeset , oDat, layout = c(nlevels(oDat$probeset), 1), groups = gType, auto.key= TRUE)

stripplot(geneExp ~ devStage | probeset, oDat, groups = gType, auto.key = TRUE, type = c('p','a'))

densityplot(~ geneExp | gType, oDat)

densityplot(~ geneExp, oDat, groups = gType)

jBw <- 0.2
jn <- 400
densityplot(~ geneExp, oDat,
            groups = devStage, auto.key = TRUE,
            bw = jBw, n = jn,
            main = paste("bw =", jBw, ", n =", jn))

bwplot(geneExp ~ devStage, oDat)

bwplot(geneExp ~ devStage | gType, oDat,panel = panel.violin)

# Heat maps 
prDat <- read.table("GSE4051_data.tsv")
str(prDat, max.level = 0)

prDes <- readRDS("GSE4051_design.rds")
str(prDes)

set.seed(1)
(yo <- sample(1:nrow(prDat), size = 50))

hDat <- prDat[yo, ]
str(hDat)

hDat <- as.matrix(t(hDat))
rownames(hDat) <- with(prDes, paste(devStage, gType, sidChar, sep="_"))
str(hDat)

heatmap(hDat)
heatmap(hDat, Rowv = NA, Colv = NA, scale="none", margins = c(5, 8))

heatmap(hDat, Rowv = NA, Colv = NA, col = cm.colors(256), scale="none", margins = c(5, 8))

