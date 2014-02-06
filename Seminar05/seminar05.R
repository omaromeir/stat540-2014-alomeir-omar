library(lattice)
prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level = 0)

prDes  <- readRDS("../data/GSE4051_design.rds")
str(prDes)

#g: the selected genes
prepareData  <- function(g){
  #prDes row + gExp of the gene from prDat + gene from the input
  pDat <- data.frame()
  for (i in 1:length(g)){
    pDat1 <- data.frame(prDes, gExp = as.vector(as.matrix(prDat[g[i], ])), gene = g[i])
    pDat <- rbind(pDat, pDat1)
  } 
  pDat
}

(magicGenes = c("1419655_at", "1438815_at"))
jDat  <- prepareData(magicGenes)
str(jDat)
head(jDat)
tail(jDat)

stripplot(gExp ~ devStage | gene, jDat, group = gType, jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE)

#the data
makeStripplot  <- function(d){
  stripplot(gExp ~ devStage | gene, d, group = gType, jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE)
}

