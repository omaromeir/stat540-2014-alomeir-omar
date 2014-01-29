library(lattice)
prDat = read.table("GSE4051_data.tsv")
str(prDat, max.level=0)
prDes = readRDS("GSE4051_design.rds")
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
