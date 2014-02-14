library(lattice)
prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level = 0)

prDes  <- readRDS("../data/GSE4051_design.rds")
str(prDes)

#g: the selected genes
prepareData  <- function(g, d = prDes){
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

#d: the data
makeStripplot  <- function(d){
  stripplot(gExp ~ devStage | gene, d, group = gType, jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE)
}

makeStripplot(jDat)

#t test
tDat <- prepareData("1456341_a_at")
str(tDat)

t.test(gExp ~ devStage, subset(tDat, devStage == "P2" | devStage == "4_weeks" ))

#Linear model
makeStripplot(mDat <- prepareData("1438786_a_at"))
str(mDat)

mFit <- lm(formula = gExp ~ devStage, data=mDat, subset = gType=="wt")

summary(mFit)

#Questions: does the intercept look plausible given the plot? How about the devStageP2 effect, etc.?
#Yes, the P2 and P10 effect is clear using E16 as an intercept.

#Iference for a contrast

contMat = matrix(c(0, 1, 0, -1, 0), nrow=1)

(obsDiff <- contMat %*% coef(mFit))

(sampMeans <- aggregate(gExp ~ devStage, mDat, FUN = mean,
                        subset = gType == "wt"))
with(sampMeans, gExp[devStage == "P2"] - gExp[devStage == "P10"])

sqrt(diag(vcov(mFit)))
summary(mFit)$coefficients[ , "Std. Error"]

(estSe <- contMat %*% vcov(mFit) %*% t(contMat))
(testStat <- obsDiff/estSe)

2 * pt(abs(testStat), df = df.residual(mFit), lower.tail = FALSE)

#Two categorical covariates

makeStripplot(oDat <- prepareData("1448690_at"))

oFitBig <- lm(formula = gExp ~ gType * devStage, data=oDat)

oFitSmall <- lm(formula = gExp ~ gType + devStage, data=oDat)

anova(oFitSmall, oFitBig)

#Is the intercept plausible? How about the various effects? Do the ones with small p-values, e.g. meeting a conventional cut-off of 0.05, look 'real' to you?
#Intercept seems reasonable, the effects with p values smaller than 0.05 are clear in the plot.

#Looking at a more interesting gene

makeStripplot(iDat <- prepareData("1429225_at"))

iFitBig <- lm(formula = gExp ~ gType * devStage, data=iDat)

iFitSmall <- lm(formula = gExp ~ gType + devStage, data=iDat)

anova(iFitSmall, iFitBig)

#Further experimentation

#adding age as a quantitative variable
library(car)

mDat$age <- recode(mDat$devStage, "'E16'=-2; 'P2'=2; ' P6'=6; 'P10'=10; '4_weeks'=28", as.factor.result = FALSE)

# xyplot(gExp ~ age, mDat, type=c("p", "smooth"))

#Quadratic model
# qFit = lm(formula=gExp~age + I(age^2), data=mDat)
# summary(qFit)
