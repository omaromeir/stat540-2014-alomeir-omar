prDat <- read.table("GSE4051_MINI.txt", header=TRUE, row.names=1)

# How many rows are there? Hint: nrow(), dim().

nrow(prDat)

# How many columns or variables are there? Hint: ncol(), length(), dim().

ncol(prDat)

dim(prDat) #returns both number of rows and number of columns

# Inspect the first few observations or the last few or a random sample. Hint: head(), tail(), x[i, j] combined with sample().

head(prDat)
tail(prDat)
prDat[sample(nrow(prDat), size = 6), ] #sample takes a random sample from the number of rows

# What does row correspond to -- different genes or different mice?

#Different genes

# What are the variable names? Hint: names(), dimnames().

names(prDat) #sample, devStage, gType, crabHammer, eggBomb, poisonFang

# What "flavor" is each variable, i.e. numeric, character, factor? Hint: str().

str(prDat) # int, Factor, Factor, num, num, num

# For sample, do a sanity check that each integer between 1 and the number of rows in the dataset occurs exactly once. Hint: a:b, seq(), seq_len(), sort(), table(), ==, all(), all.equal(), identical().

identical(seq(1, nrow(prDat)), sort(prDat[["sample"]]))

# For each factor variable, what are the levels? Hint: levels(), str().

str(prDat) # devStage: 5. gType: 2

# How many observations do we have for each level of devStage? For gType? Hint: summary(), table().

summary(prDat) # devStage: 8, 7, 8, 8, 8. gType: 19, 20

# Perform a cross-tabulation of devStage and gType. Hint: table().

table(c(prDat["devStage"], prDat["gType"]))

# If you had to take a wild guess, what do you think the intended experimental design was? What actually happened in real life?

# Analyzing data of photo receptor cells in mice. Various developmental stages, two genotypes for 3 gene expressions.

# For each quantitative variable, what are the extremes? How about average or median? Hint: min(), max(), range(), summary(), fivenum(), mean(), median(), quantile().

summary(prDat$sample)
summary(prDat$crabHammer)
summary(prDat$eggBomb)
summary(prDat$poisonFang)

# Create a new data.frame called weeDat only containing observations for which expression of poisonFang is above 7.5.

weeDat  <- subset(prDat, subset = poisonFang > 7.5)

# For how many observations poisonFang > 7.5? How do they break down by genotype and developmental stage?

nrow(weeDat) # 9. In terms of gType, 5 are wt and 4 are NrlKO. In terms of devStage, 4 are P2, 1 is P6, 4 are P10

# Print the observations with row names "Sample_16" and "Sample_38" to screen, showing only the 3 gene expression variables.

(subset(prDat, subset = sample == 38, select = c("crabHammer", "eggBomb", "poisonFang")))
(subset(prDat, subset = sample == 16, select = c("crabHammer", "eggBomb", "poisonFang")))

# Which samples have expression of eggBomb less than the 0.10 quantile?

subset(prDat, subset = eggBomb < quantile(prDat$eggBomb, probs = seq(0, 1, 0.1))[2], select = "sample") # samples: 25, 14, 3, 35
