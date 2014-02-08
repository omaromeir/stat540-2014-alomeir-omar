library(lattice)
prDat = read.table("../data/GSE4051_data.tsv")
str(prDat, max.level=0)
prDes = readRDS("../data/GSE4051_design.rds")
str(prDes)

set.seed(987)
(theGene <- sample(1:nrow(prDat), 1))
pDat  <- data.frame(prDes, gExp = unlist(prDat[theGene, ]))
str(pDat)

aggregate(gExp ~ gType, pDat, FUN = mean)
stripplot(gType ~ gExp, pDat)
ttRes  <- t.test(gExp ~ gType, pDat)
str(ttRes)
ttRes$statistic
ttRes$p.value

# You try

myGene <- 14158
pDat2  <- data.frame(prDes, gExp = unlist(prDat[myGene, ]))
str(pDat2)
aggregate(gExp ~ gType, pDat2, FUN = mean)
stripplot(gType ~ gExp, pDat2)
tt  <- t.test(gExp ~ gType, pDat2)
str(tt)
tt$statistic
tt$p.value

tt2  <- t.test(gExp ~ gType, pDat2, var.equal=TRUE)
str(tt2)

wt <- wilcox.test(gExp ~ gType, pDat2)
str(wt)
wt$statistic
wt$p.value

ks <- ks.test(pDat2$gExp[pDat2$gType == "wt"], pDat2$gExp[pDat2$gType == "NrlKO"])
str(ks)
ks$statistic
ks$p.value

# Can you pull test statistics and/or p-values from the different approaches into an common object, like a readable table? Are you getting the same message from the various approaches?

rbind(c(tt$statistic, tt$p.value, wt$statistic, wt$p.value, ks$statistic, ks$p.value))
#According to p values all approaches show that there is no difference in means.


#apply() stuff from here on out
kDat <- readRDS("../data/GSE4051_MINI.rds")
kMat <- as.matrix(kDat[c('crabHammer', 'eggBomb', 'poisonFang')])
str(kMat)
median(kMat[, "eggBomb"])
apply(kMat, 2, median)

apply(kMat, 1, min)
colnames(kMat)[apply(kMat, 1, which.min)]

#aggregate() stuff
aggregate(eggBomb ~ devStage, kDat, FUN = mean)
aggregate(eggBomb ~ gType * devStage, kDat, FUN = mean)

#two sample test
keepGenes <- c("1431708_a_at", "1424336_at", "1454696_at","1416119_at", "1432141_x_at", "1429226_at" )
miniDat <- subset(prDat, rownames(prDat) %in% keepGenes)
miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
                      gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                    levels = keepGenes))
miniDat <- suppressWarnings(data.frame(prDes, miniDat))
str(miniDat)
stripplot(gType ~ gExp | gene, miniDat, scales = list(x = list(relation = "free")), group = gType, auto.key = TRUE)

someDat <- droplevels(subset(miniDat, gene == keepGenes[1]))
t.test(gExp ~ gType, someDat)

#plyr stuff
library(plyr)

d_ply(miniDat, ~ gene, function(x) t.test(gExp ~ gType, x), .print=TRUE)
ttRes <- dlply(miniDat, ~ gene, function(x) t.test(gExp ~ gType, x))
names(ttRes)

ttRes <- ddply(miniDat, ~ gene, function(z) {
  zz <- t.test(gExp ~ gType, z)
  round(c(tStat = zz$statistic, pVal = zz$p.value), 4)
})
ttRes

