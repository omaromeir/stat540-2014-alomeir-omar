Omar AlOmeir

Seminar 01
========================================================



```r
prDat <- read.table("GSE4051_MINI.txt", header = TRUE, row.names = 1)
```


How many rows are there? Hint: nrow(), dim().


```r
nrow(prDat)
```

```
## [1] 39
```


How many columns or variables are there? Hint: ncol(), length(), dim().


```r
ncol(prDat)
```

```
## [1] 6
```

```r

dim(prDat)  #returns both number of rows and number of columns
```

```
## [1] 39  6
```


Inspect the first few observations or the last few or a random sample. Hint: head(), tail(), x[i, j] combined with sample().


```r
head(prDat)
```

```
##           sample devStage gType crabHammer eggBomb poisonFang
## Sample_20     20      E16    wt     10.220   7.462      7.370
## Sample_21     21      E16    wt     10.020   6.890      7.177
## Sample_22     22      E16    wt      9.642   6.720      7.350
## Sample_23     23      E16    wt      9.652   6.529      7.040
## Sample_16     16      E16 NrlKO      8.583   6.470      7.494
## Sample_17     17      E16 NrlKO     10.140   7.065      7.005
```

```r
tail(prDat)
```

```
##           sample devStage gType crabHammer eggBomb poisonFang
## Sample_38     38  4_weeks    wt      9.767   6.608      7.329
## Sample_39     39  4_weeks    wt     10.200   7.003      7.320
## Sample_11     11  4_weeks NrlKO      9.677   7.204      6.981
## Sample_12     12  4_weeks NrlKO      9.129   7.165      7.350
## Sample_2       2  4_weeks NrlKO      9.744   7.107      7.075
## Sample_9       9  4_weeks NrlKO      9.822   6.558      7.043
```

```r
prDat[sample(nrow(prDat), size = 6), ]  #sample takes a random sample from the number of rows
```

```
##           sample devStage gType crabHammer eggBomb poisonFang
## Sample_17     17      E16 NrlKO     10.140   7.065      7.005
## Sample_12     12  4_weeks NrlKO      9.129   7.165      7.350
## Sample_36     36  4_weeks    wt      9.960   7.866      6.993
## Sample_28     28       P6    wt      8.214   6.530      7.428
## Sample_37     37  4_weeks    wt      9.667   6.992      7.324
## Sample_7       7       P6 NrlKO      8.803   6.188      7.754
```


What does row correspond to -- different genes or different mice?
Different genes

What are the variable names? Hint: names(), dimnames().


```r
names(prDat)  #sample, devStage, gType, crabHammer, eggBomb, poisonFang
```

```
## [1] "sample"     "devStage"   "gType"      "crabHammer" "eggBomb"   
## [6] "poisonFang"
```


What "flavor" is each variable, i.e. numeric, character, factor? Hint: str().


```r
str(prDat)  # int, Factor, Factor, num, num, num
```

```
## 'data.frame':	39 obs. of  6 variables:
##  $ sample    : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```


For sample, do a sanity check that each integer between 1 and the number of rows in the dataset occurs exactly once. Hint: a:b, seq(), seq_len(), sort(), table(), ==, all(), all.equal(), identical().


```r
identical(seq(1, nrow(prDat)), sort(prDat[["sample"]]))
```

```
## [1] TRUE
```


For each factor variable, what are the levels? Hint: levels(), str().


```r
str(prDat)  # devStage: 5. gType: 2
```

```
## 'data.frame':	39 obs. of  6 variables:
##  $ sample    : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```


How many observations do we have for each level of devStage? For gType? Hint: summary(), table().


```r
summary(prDat)  # devStage: 8, 7, 8, 8, 8. gType: 19, 20
```

```
##      sample        devStage   gType      crabHammer       eggBomb    
##  Min.   : 1.0   4_weeks:8   NrlKO:19   Min.   : 8.21   Min.   :6.14  
##  1st Qu.:10.5   E16    :7   wt   :20   1st Qu.: 8.94   1st Qu.:6.28  
##  Median :20.0   P10    :8              Median : 9.61   Median :6.76  
##  Mean   :20.0   P2     :8              Mean   : 9.43   Mean   :6.79  
##  3rd Qu.:29.5   P6     :8              3rd Qu.: 9.83   3rd Qu.:7.09  
##  Max.   :39.0                          Max.   :10.34   Max.   :8.17  
##    poisonFang  
##  Min.   :6.74  
##  1st Qu.:7.19  
##  Median :7.35  
##  Mean   :7.38  
##  3rd Qu.:7.48  
##  Max.   :8.58
```


Perform a cross-tabulation of devStage and gType. Hint: table().


```r
table(c(prDat["devStage"], prDat["gType"]))
```

```
##          gType
## devStage  NrlKO wt
##   4_weeks     4  4
##   E16         3  4
##   P10         4  4
##   P2          4  4
##   P6          4  4
```


If you had to take a wild guess, what do you think the intended experimental design was? What actually happened in real life?

Analyzing data of photo receptor cells in mice. Various developmental stages, two genotypes for 3 gene expressions.

For each quantitative variable, what are the extremes? How about average or median? Hint: min(), max(), range(), summary(), fivenum(), mean(), median(), quantile().


```r
summary(prDat$sample)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     1.0    10.5    20.0    20.0    29.5    39.0
```

```r
summary(prDat$crabHammer)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    8.21    8.94    9.61    9.43    9.83   10.30
```

```r
summary(prDat$eggBomb)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    6.14    6.28    6.76    6.79    7.09    8.17
```

```r
summary(prDat$poisonFang)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    6.74    7.19    7.35    7.38    7.48    8.58
```


Create a new data.frame called weeDat only containing observations for which expression of poisonFang is above 7.5.


```r
weeDat <- subset(prDat, subset = poisonFang > 7.5)
```


For how many observations poisonFang > 7.5? How do they break down by genotype and developmental stage?


```r
nrow(weeDat)  # 9. In terms of gType, 5 are wt and 4 are NrlKO. In terms of devStage, 4 are P2, 1 is P6, 4 are P10
```

```
## [1] 9
```


Print the observations with row names "Sample_16" and "Sample_38" to screen, showing only the 3 gene expression variables.


```r
(subset(prDat, subset = sample == 38, select = c("crabHammer", "eggBomb", "poisonFang")))
```

```
##           crabHammer eggBomb poisonFang
## Sample_38      9.767   6.608      7.329
```

```r
(subset(prDat, subset = sample == 16, select = c("crabHammer", "eggBomb", "poisonFang")))
```

```
##           crabHammer eggBomb poisonFang
## Sample_16      8.583    6.47      7.494
```


Which samples have expression of eggBomb less than the 0.10 quantile?


```r
subset(prDat, subset = eggBomb < quantile(prDat$eggBomb, probs = seq(0, 1, 0.1))[2], 
    select = "sample")  # samples: 25, 14, 3, 35
```

```
##           sample
## Sample_25     25
## Sample_14     14
## Sample_3       3
## Sample_35     35
```


