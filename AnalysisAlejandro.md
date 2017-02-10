---
output: html_document
---
<!--
To compile run:
library(knitr); knit2html( "AnalysisAlejandro.Rmd" )
To extract the R-code:
knitr::purl("AnalysisAlejandro.Rmd")
-->

Exon usage data: FDR analysis of tissue dependent  exon usage between species pairs
======================

 **LAST UPDATE AT**


```
## [1] "Sat Aug 23 14:53:03 2014"
```

Preparations
--------------

We first set global chunk options and load the
necessary packages.


```r
library(knitr)
library(reshape2)
library(ConcaveFDR)
library(ggplot2)
library(grid)
library(locfdr)
library(fdrtool)
library(mixfdr)
library(BiocParallel)
library(plyr)
library(grid)
library(RColorBrewer)
library(psych)

opts_chunk$set(fig.width=14, fig.height=10, cache=T, 
error=FALSE, cache.lazy=FALSE)
options(digits = 2, scipen = 10)

multicoreParam <- MulticoreParam(workers = round(parallel::detectCores() * 0.75),
                                 catch.errors = FALSE)
multicoreParam
```

```
## class: MulticoreParam; bpisup: TRUE; bpworkers: 60; catch.errors: FALSE
## setSeed: TRUE; recursive: TRUE; cleanup: TRUE; cleanupSignal: 15;
##   verbose: FALSE
```

```r
register(multicoreParam)
```

# Get the input data

All the data can be downloaded from [http://www-huber.embl.de].

```r
input.data = list.files("/g/huber/www-huber/pub/DEUprimates_supplement/SupplementaryFile3/objects")
input.data
```

```
## [1] "allEffects.RData"        "back-coding.RData"      
## [3] "back-exon.RData"         "CDUT-class.RData"       
## [5] "classesDefinition.RData" "ctsr.RData"             
## [7] "explainedVariance.RData" "pcaData.RData"          
## [9] "REUC_analysis.rda"
```

```r
## load all input data into the global workspace
invisible(lapply( paste0("/g/huber/www-huber/pub/DEUprimates_supplement/SupplementaryFile3/objects/",input.data), load, envir = .GlobalEnv ))
```


# Get correlations

We first inspect the data. What we analyze is the covariances/correlations 
across the 5  tissues of relative exon usages for an exon i for a given 
species pair. We first inspect the species pairs and the covariances / correlations.


```r
head( spPairs )
```

```
##         [,1]  [,2] 
## ggo:hsa "ggo" "hsa"
## ggo:mml "ggo" "mml"
## ggo:ppa "ggo" "ppa"
## ggo:ppy "ggo" "ppy"
## ggo:ptr "ggo" "ptr"
## hsa:mml "hsa" "mml"
```

```r
head( sppCovs, 20)
```

```
##                      ggo:hsa  ggo:mml  ggo:ppa  ggo:ppy  ggo:ptr   hsa:mml
## ENSG00000000419:E005  0.0072 -0.01675  0.00245  0.01296  0.00752 -0.013780
## ENSG00000000419:E010 -0.0047 -0.00099  0.01169  0.01206  0.01620 -0.001765
## ENSG00000000419:E011  0.0374 -0.00542  0.02842  0.01475  0.01353 -0.004804
## ENSG00000000419:E013  0.0107 -0.01828 -0.00395  0.01932 -0.00135 -0.005301
## ENSG00000000419:E015  0.0069  0.01346  0.01365  0.01876  0.00991  0.003715
## ENSG00000000419:E016  0.0408 -0.03285  0.02408  0.02178 -0.00404 -0.033239
## ENSG00000000457:E005  0.0077  0.01138  0.00649  0.00358  0.00280  0.042168
## ENSG00000000457:E006 -0.0109 -0.00356  0.02263 -0.01307  0.02320 -0.005847
## ENSG00000000457:E007 -0.0179 -0.00892  0.01169  0.01933  0.02158 -0.000049
## ENSG00000000457:E008  0.0051 -0.00069  0.00339  0.00537  0.00357 -0.017499
## ENSG00000000457:E009  0.0152  0.00203  0.01724  0.01512 -0.01111 -0.002204
## ENSG00000000457:E010  0.0054 -0.00761  0.00875  0.00545 -0.00177 -0.007603
## ENSG00000000457:E011  0.0034  0.00013 -0.00135 -0.00175  0.00273 -0.008569
## ENSG00000000457:E012  0.0090 -0.00173  0.01191 -0.02280 -0.00199  0.014533
## ENSG00000000457:E013  0.0013  0.00528  0.00077  0.00248  0.00087  0.002297
## ENSG00000000457:E014  0.0164  0.00147  0.01221  0.01276 -0.01316  0.001811
## ENSG00000000457:E015  0.0032  0.01126  0.00631  0.03412  0.01541  0.003870
## ENSG00000000457:E016  0.0212  0.00561 -0.02169  0.01150  0.01065  0.018098
## ENSG00000000457:E017  0.0076  0.00224  0.01194  0.00045  0.00237 -0.004977
## ENSG00000000457:E021      NA       NA       NA       NA       NA        NA
##                       hsa:ppa  hsa:ppy hsa:ptr  mml:ppa  mml:ppy  mml:ptr
## ENSG00000000419:E005  0.00202  0.00947 -0.0041  0.00420 -0.00032 -0.01577
## ENSG00000000419:E010  0.00409 -0.00589  0.0058 -0.00260 -0.00338 -0.00549
## ENSG00000000419:E011  0.03624 -0.00045  0.0098 -0.00194 -0.00694 -0.01047
## ENSG00000000419:E013 -0.00075  0.00977  0.0084  0.00515 -0.01334  0.00249
## ENSG00000000419:E015  0.00739  0.00372  0.0044  0.00851  0.01360  0.00636
## ENSG00000000419:E016  0.01925  0.02045  0.0102 -0.02064 -0.01176 -0.00623
## ENSG00000000457:E005  0.03232  0.02707  0.0214  0.03770  0.06971  0.06208
## ENSG00000000457:E006 -0.00518 -0.00563 -0.0066 -0.01309 -0.00631 -0.00523
## ENSG00000000457:E007 -0.01334 -0.00643 -0.0172  0.00696 -0.00760  0.00844
## ENSG00000000457:E008  0.04794  0.06142  0.0024 -0.01902 -0.00861 -0.00022
## ENSG00000000457:E009  0.03261  0.02050  0.0035  0.00038 -0.00374 -0.00145
## ENSG00000000457:E010  0.01091  0.00582  0.0038 -0.02750  0.00025 -0.01164
## ENSG00000000457:E011 -0.00705 -0.00592 -0.0122  0.00372 -0.00474  0.01164
## ENSG00000000457:E012  0.01253 -0.02692  0.0066  0.00570 -0.01257  0.00660
## ENSG00000000457:E013  0.01294  0.01528 -0.0050  0.01120  0.01237  0.00450
## ENSG00000000457:E014  0.01619  0.02508 -0.0371 -0.00725 -0.00415 -0.00433
## ENSG00000000457:E015 -0.00142  0.00170 -0.0128  0.00106  0.00976  0.00692
## ENSG00000000457:E016 -0.01305  0.00692 -0.0036  0.00163  0.00713  0.00440
## ENSG00000000457:E017  0.01181  0.00715  0.0045 -0.00405 -0.01035 -0.01904
## ENSG00000000457:E021       NA       NA      NA       NA       NA       NA
##                       ppa:ppy  ppa:ptr  ppy:ptr
## ENSG00000000419:E005  0.00402 -0.00966 -0.00337
## ENSG00000000419:E010  0.00468  0.01658  0.01015
## ENSG00000000419:E011  0.01119  0.02626  0.01718
## ENSG00000000419:E013 -0.00100  0.00017  0.00014
## ENSG00000000419:E015  0.01643  0.01384  0.01165
## ENSG00000000419:E016  0.01521  0.00231  0.00867
## ENSG00000000457:E005  0.02870  0.02210  0.04752
## ENSG00000000457:E006  0.00799  0.01765 -0.00396
## ENSG00000000457:E007 -0.00166  0.02302  0.00711
## ENSG00000000457:E008  0.06464 -0.01112 -0.02336
## ENSG00000000457:E009  0.02536  0.00453 -0.01816
## ENSG00000000457:E010 -0.00072  0.01587 -0.00933
## ENSG00000000457:E011  0.00590  0.01012  0.00322
## ENSG00000000457:E012 -0.01843  0.00223 -0.00517
## ENSG00000000457:E013  0.01766 -0.00355  0.00291
## ENSG00000000457:E014  0.02167 -0.00994 -0.01815
## ENSG00000000457:E015 -0.00164  0.02672  0.00339
## ENSG00000000457:E016  0.01241  0.00932  0.00847
## ENSG00000000457:E017  0.00021  0.02255  0.00277
## ENSG00000000457:E021       NA       NA       NA
```

```r
str( sppCovs )
```

```
##  num [1:119344, 1:15] 0.00724 -0.00471 0.03738 0.01067 0.00687 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:119344] "ENSG00000000419:E005" "ENSG00000000419:E010" "ENSG00000000419:E011" "ENSG00000000419:E013" ...
##   ..$ : chr [1:15] "ggo:hsa" "ggo:mml" "ggo:ppa" "ggo:ppy" ...
```

```r
head( sppCors, 20)
```

```
##                      ggo:hsa ggo:mml ggo:ppa ggo:ppy ggo:ptr hsa:mml
## ENSG00000000419:E005    0.47  -0.454   0.231   0.714   0.318 -0.5097
## ENSG00000000419:E010   -0.26  -0.121   0.528   0.595   0.617 -0.3531
## ENSG00000000419:E011    0.72  -0.169   0.401   0.387   0.333 -0.2137
## ENSG00000000419:E013    0.56  -0.384  -0.677   0.538  -0.044 -0.2776
## ENSG00000000419:E015    0.61   0.860   0.645   0.909   0.803  0.3823
## ENSG00000000419:E016    0.93  -0.593   0.700   0.813  -0.093 -0.6736
## ENSG00000000457:E005    0.79   0.618   0.651   0.283   0.253  0.6827
## ENSG00000000457:E006   -0.34  -0.177   0.660  -0.435   0.960 -0.3518
## ENSG00000000457:E007   -0.86  -0.473   0.329   0.875   0.556 -0.0048
## ENSG00000000457:E008    0.75  -0.176   0.496   0.516   0.612 -0.5994
## ENSG00000000457:E009    0.67   0.248   0.687   0.828  -0.543 -0.1402
## ENSG00000000457:E010    0.60  -0.668   0.635   0.728  -0.167 -0.4228
## ENSG00000000457:E011    0.38   0.016  -0.184  -0.209   0.263 -0.6249
## ENSG00000000457:E012    0.37  -0.080   0.760  -0.727  -0.228  0.7078
## ENSG00000000457:E013    0.15   0.672   0.071   0.239   0.166  0.1412
## ENSG00000000457:E014    0.81   0.186   0.673   0.763  -0.522  0.1287
## ENSG00000000457:E015    0.14   0.648   0.257   0.937   0.281  0.4132
## ENSG00000000457:E016    0.47   0.207  -0.567   0.317   0.314  0.8252
## ENSG00000000457:E017    0.69   0.151   0.637   0.043   0.120 -0.3279
## ENSG00000000457:E021      NA      NA      NA      NA      NA      NA
##                      hsa:ppa hsa:ppy hsa:ptr mml:ppa mml:ppy mml:ptr
## ENSG00000000419:E005    0.26   0.712  -0.237   0.224 -0.0098  -0.377
## ENSG00000000419:E010    0.30  -0.474   0.363  -0.427 -0.6057  -0.759
## ENSG00000000419:E011    0.73  -0.017   0.345  -0.064 -0.4224  -0.598
## ENSG00000000419:E013   -0.32   0.680   0.685   0.874 -0.3689   0.080
## ENSG00000000419:E015    0.56   0.290   0.576   0.470  0.7698   0.602
## ENSG00000000419:E016    0.63   0.857   0.262  -0.535 -0.3914  -0.128
## ENSG00000000457:E005    0.97   0.640   0.577   0.598  0.8737   0.888
## ENSG00000000457:E006   -0.18  -0.227  -0.328  -0.733 -0.4031  -0.416
## ENSG00000000457:E007   -0.68  -0.532  -0.809   0.394 -0.6926   0.438
## ENSG00000000457:E008    0.95   0.795   0.055  -0.651 -0.1934  -0.009
## ENSG00000000457:E009    0.68   0.585   0.090   0.022 -0.2944  -0.102
## ENSG00000000457:E010    0.50   0.493   0.227  -0.992  0.0169  -0.548
## ENSG00000000457:E011   -0.57  -0.422  -0.701   0.336 -0.3758   0.744
## ENSG00000000457:E012    0.84  -0.906   0.801   0.427 -0.4711   0.890
## ENSG00000000457:E013    0.58   0.712  -0.458   0.559  0.6463   0.467
## ENSG00000000457:E014    0.50   0.841  -0.826  -0.577 -0.3586  -0.248
## ENSG00000000457:E015   -0.11   0.086  -0.434   0.103  0.6444   0.303
## ENSG00000000457:E016   -0.42   0.236  -0.132   0.087  0.4002   0.265
## ENSG00000000457:E017    0.62   0.675   0.223  -0.157 -0.7249  -0.701
## ENSG00000000457:E021      NA      NA      NA      NA      NA      NA
##                      ppa:ppy ppa:ptr ppy:ptr
## ENSG00000000419:E005   0.436  -0.804  -0.164
## ENSG00000000419:E010   0.309   0.844   0.565
## ENSG00000000419:E011   0.309   0.680   0.826
## ENSG00000000419:E013  -0.225   0.046   0.006
## ENSG00000000419:E015   0.688   0.969   0.837
## ENSG00000000419:E016   0.815   0.076   0.368
## ENSG00000000457:E005   0.665   0.585   0.993
## ENSG00000000457:E006   0.299   0.822  -0.211
## ENSG00000000457:E007  -0.080   0.634   0.315
## ENSG00000000457:E008   0.836  -0.257  -0.354
## ENSG00000000457:E009   0.651   0.104  -0.572
## ENSG00000000457:E010  -0.040   0.618  -0.669
## ENSG00000000457:E011   0.520   0.720   0.201
## ENSG00000000457:E012  -0.955   0.415  -0.482
## ENSG00000000457:E013   0.668  -0.267   0.229
## ENSG00000000457:E014   0.813  -0.248  -0.491
## ENSG00000000457:E015  -0.077   0.828   0.071
## ENSG00000000457:E016   0.494   0.398   0.381
## ENSG00000000457:E017   0.012   0.659   0.146
## ENSG00000000457:E021      NA      NA      NA
```

```r
str( sppCors )
```

```
##  num [1:119344, 1:15] 0.473 -0.26 0.716 0.563 0.605 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:119344] "ENSG00000000419:E005" "ENSG00000000419:E010" "ENSG00000000419:E011" "ENSG00000000419:E013" ...
##   ..$ : chr [1:15] "ggo:hsa" "ggo:mml" "ggo:ppa" "ggo:ppy" ...
```

# get significant covariances according to the paper 


```r
head( padjCov, 20)
```

```
##                      ggo:hsa ggo:mml ggo:ppa ggo:ppy ggo:ptr hsa:mml
## ENSG00000000419:E005   0.283    0.89   0.433   0.366    0.54    0.92
## ENSG00000000419:E010   0.638    0.67   0.204   0.385    0.35    0.84
## ENSG00000000419:E011   0.023    0.77   0.040   0.327    0.41    0.87
## ENSG00000000419:E013   0.214    0.90   0.639   0.239    0.75    0.88
## ENSG00000000419:E015   0.291    0.35   0.170   0.251    0.48    0.77
## ENSG00000000419:E016   0.018    0.95   0.063   0.199    0.80    0.95
## ENSG00000000457:E005   0.272    0.39   0.315   0.609    0.65    0.13
## ENSG00000000457:E006   0.759    0.73   0.073   0.884    0.23    0.88
## ENSG00000000457:E007   0.840    0.82   0.204   0.239    0.26    0.82
## ENSG00000000457:E008   0.336    0.67   0.403   0.562    0.63    0.93
## ENSG00000000457:E009   0.145    0.59   0.123   0.319    0.87    0.85
## ENSG00000000457:E010   0.329    0.80   0.258   0.560    0.76    0.89
## ENSG00000000457:E011   0.380    0.65   0.565   0.748    0.65    0.90
## ENSG00000000457:E012   0.242    0.69   0.201   0.929    0.76    0.58
## ENSG00000000457:E013   0.447    0.51   0.492   0.639    0.70    0.78
## ENSG00000000457:E014   0.130    0.61   0.195   0.371    0.89    0.79
## ENSG00000000457:E015   0.388    0.39   0.320   0.075    0.37    0.76
## ENSG00000000457:E016   0.084    0.51   0.876   0.400    0.47    0.51
## ENSG00000000457:E017   0.276    0.59   0.200   0.688    0.66    0.87
## ENSG00000000457:E021      NA      NA      NA      NA      NA      NA
##                      hsa:ppa hsa:ppy hsa:ptr mml:ppa mml:ppy mml:ptr
## ENSG00000000419:E005   0.360   0.391    0.73    0.75   0.690   0.899
## ENSG00000000419:E010   0.300   0.772    0.48    0.85   0.757   0.828
## ENSG00000000419:E011   0.018   0.663    0.39    0.84   0.812   0.869
## ENSG00000000419:E013   0.452   0.384    0.42    0.73   0.875   0.710
## ENSG00000000419:E015   0.225   0.541    0.52    0.67   0.366   0.621
## ENSG00000000419:E016   0.075   0.182    0.38    0.93   0.863   0.836
## ENSG00000000457:E005   0.024   0.111    0.17    0.17   0.021   0.026
## ENSG00000000457:E006   0.585   0.768    0.77    0.91   0.803   0.827
## ENSG00000000457:E007   0.741   0.780    0.87    0.70   0.821   0.576
## ENSG00000000457:E008   0.009   0.015    0.58    0.93   0.833   0.764
## ENSG00000000457:E009   0.024   0.182    0.55    0.81   0.763   0.784
## ENSG00000000457:E010   0.163   0.483    0.54    0.95   0.676   0.878
## ENSG00000000457:E011   0.629   0.772    0.83    0.75   0.779   0.509
## ENSG00000000457:E012   0.141   0.924    0.46    0.72   0.870   0.618
## ENSG00000000457:E013   0.137   0.268    0.74    0.62   0.395   0.662
## ENSG00000000457:E014   0.101   0.131    0.93    0.89   0.771   0.817
## ENSG00000000457:E015   0.475   0.600    0.84    0.80   0.449   0.611
## ENSG00000000457:E016   0.737   0.455    0.72    0.79   0.508   0.665
## ENSG00000000457:E017   0.150   0.448    0.52    0.86   0.850   0.911
## ENSG00000000457:E021      NA      NA      NA      NA      NA      NA
##                      ppa:ppy ppa:ptr ppy:ptr
## ENSG00000000419:E005   0.472   0.729   0.827
## ENSG00000000419:E010   0.450   0.135   0.542
## ENSG00000000419:E011   0.286   0.055   0.373
## ENSG00000000419:E013   0.630   0.492   0.772
## ENSG00000000419:E015   0.190   0.170   0.503
## ENSG00000000419:E016   0.213   0.423   0.582
## ENSG00000000457:E005   0.066   0.082   0.038
## ENSG00000000457:E006   0.361   0.123   0.833
## ENSG00000000457:E007   0.649   0.074   0.620
## ENSG00000000457:E008   0.014   0.751   0.928
## ENSG00000000457:E009   0.089   0.357   0.914
## ENSG00000000457:E010   0.622   0.143   0.878
## ENSG00000000457:E011   0.414   0.231   0.709
## ENSG00000000457:E012   0.881   0.425   0.845
## ENSG00000000457:E013   0.169   0.602   0.716
## ENSG00000000457:E014   0.121   0.733   0.914
## ENSG00000000457:E015   0.649   0.053   0.704
## ENSG00000000457:E016   0.262   0.247   0.586
## ENSG00000000457:E017   0.591   0.077   0.718
## ENSG00000000457:E021      NA      NA      NA
```

```r
str( padjCov )
```

```
##  num [1:119344, 1:15] 0.2828 0.6382 0.0234 0.2138 0.2911 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:119344] "ENSG00000000419:E005" "ENSG00000000419:E010" "ENSG00000000419:E011" "ENSG00000000419:E013" ...
##   ..$ : chr [1:15] "ggo:hsa" "ggo:mml" "ggo:ppa" "ggo:ppy" ...
```

```r
table(na.omit(padjCov) < 0.1)
```

```
## 
##   FALSE    TRUE 
## 1245613  128297
```

```r
# prop.table(table(na.omit(padjCov) < 0.1))
```


# run the FDR analyses for correlations

Does not work well!

```r
test.stats <- as.vector(na.omit(sppCors))
qplot(test.stats)
```

```
## stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust this.
```

![plot of chunk run FDR analyses corrs](figure/run FDR analyses corrs1.png) 

```r
sub.sample <- round(length(test.stats)*0.7)
sub.sample
```

```
## [1] 961758
```

```r
res.fdrtool <- fdrtool(sample(test.stats, sub.sample), statistic = "correlation")
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
## Step 5... prepare for plotting
```

![plot of chunk run FDR analyses corrs](figure/run FDR analyses corrs2.png) 

```r
sum(res.fdrtool$lfdr < 0.2) / sub.sample
```

```
## [1] 0.00035
```

```r
sum(res.fdrtool$lfdr < 0.5) / sub.sample
```

```
## [1] 0.0059
```

```r
sum(res.fdrtool$qval < 0.2) / sub.sample
```

```
## [1] 0.00079
```

```r
#res.ConcaveFDR <- ConcaveFDR(sample(test.stats, sub.sample), statistic = "correlation")
#res.ConcaveFDR.non.sub <- ConcaveFDR(test.stats, statistic = "correlation")
```


We can also run the FDR analysis on the z-transformed correlations.

```r
### z-value transform 
test.stats.z <- fisherz(test.stats)
test.stats.z <- test.stats.z - median(test.stats.z)
summary(test.stats.z)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    -Inf       0       0               0     Inf
```

```r
test.stats.z <- test.stats.z[abs(test.stats.z) < 4]
summary(test.stats.z)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    -4.0    -0.5     0.0     0.0     0.5     4.0
```

```r
mad(test.stats.z)
```

```
## [1] 0.7
```

```r
sub.sample <- round(length(test.stats.z)*0.1)
sub.sample
```

```
## [1] 137338
```

```r
res.fdrtool <- fdrtool(sample(test.stats.z, sub.sample), statistic = "normal")
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
## Step 5... prepare for plotting
```

![plot of chunk run FDR analyses z-transformed corrs](figure/run FDR analyses z-transformed corrs1.png) 

```r
sum(res.fdrtool$lfdr < 0.2) / sub.sample
```

```
## [1] 0.0083
```

```r
sum(res.fdrtool$lfdr < 0.5) / sub.sample
```

```
## [1] 0.021
```

```r
sum(res.fdrtool$qval < 0.2) / sub.sample
```

```
## [1] 0.016
```

```r
### use transformed values

concaveCorr <- ConcaveFDR(sample(test.stats.z, sub.sample) , statistic = "normal", 
                            cutoff.method = "smoothing")
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
## Step 5... prepare for plotting
```

![plot of chunk run FDR analyses z-transformed corrs](figure/run FDR analyses z-transformed corrs2.png) 

```r
#ConcaveFDR:::smoothing.cutoff(test.stats.z, "normal")

sum(concaveCorr$lfdr.log < 0.2)  / sub.sample
```

```
## [1] 0.014
```

```r
sum(concaveCorr$lfdr.log < 0.5)  / sub.sample
```

```
## [1] 0.043
```

```r
sum(concaveCorr$qval.log < 0.2)  / sub.sample
```

```
## [1] 0.029
```

```r
sum(concaveCorr$qval.log < 0.1)  / sub.sample
```

```
## [1] 0.015
```


# run the FDR analyses for z-transformed covariances

This agrees with the paper results of 10% conserved exons.

## fdrtool

```r
### try on covariances on the natural log scale
test.stats.cov <- as.vector(na.omit(sppCovs))
test.stats.cov <- fisherz(test.stats.cov)
```

```
## Warning: NaNs produced
```

```r
test.stats.cov <- na.omit(test.stats.cov)
test.stats.cov <- test.stats.cov - median(test.stats.cov)
test.stats.cov <- as.vector(test.stats.cov)
##save(test.stats.cov, file ="AlejandroCovs.rda")
#load("AlejandroCovs.rda")

sd(test.stats.cov)
```

```
## [1] 0.044
```

```r
summary(test.stats.cov)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    -0.6     0.0     0.0     0.0     0.0     4.2
```

```r
test.stats.cov <- test.stats.cov[abs(test.stats.cov) < 4]
summary(test.stats.z)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    -4.0    -0.5     0.0     0.0     0.5     4.0
```

```r
mad(test.stats.z)
```

```
## [1] 0.7
```

```r
sub.sample <- round(length(test.stats.cov)*0.1)
sub.sample
```

```
## [1] 137369
```

```r
FDR.cov <- fdrtool(sample(test.stats.cov, sub.sample), 
                               statistic = "normal")
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
## Step 5... prepare for plotting
```

![plot of chunk run FDR analyses covs](figure/run FDR analyses covs.png) 

```r
sum(FDR.cov$lfdr < 0.2) / sub.sample
```

```
## [1] 0.078
```

```r
sum(FDR.cov$lfdr < 0.5) / sub.sample
```

```
## [1] 0.1
```

```r
sum(FDR.cov$qval < 0.2) / sub.sample
```

```
## [1] 0.13
```

```r
sum(FDR.cov$qval < 0.1)/ sub.sample
```

```
## [1] 0.1
```


## ConcaveFDR

```r
res.ConcaveFDR.cov  <-ConcaveFDR(sample(test.stats.cov, sub.sample),
                               statistic = "normal", cutoff = "smoothing")
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
## Step 5... prepare for plotting
```

![plot of chunk Concave FDR on covariances](figure/Concave FDR on covariances.png) 

```r
sum(res.ConcaveFDR.cov$lfdr.log < 0.2) / sub.sample
```

```
## [1] 0.15
```

```r
sum(res.ConcaveFDR.cov$lfdr.log < 0.5) / sub.sample
```

```
## [1] 0.22
```

```r
sum(res.ConcaveFDR.cov$qval.log < 0.2) / sub.sample
```

```
## [1] 0.24
```

```r
sum(res.ConcaveFDR.cov$qval.log < 0.1) / sub.sample
```

```
## [1] 0.18
```

**SESSION INFO**
 


```
## R Under development (unstable) (2014-08-04 r66305)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] grid      splines   stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] psych_1.4.8.11       RColorBrewer_1.0-5   plyr_1.8.1          
##  [4] BiocParallel_0.6.1   mixfdr_1.0           locfdr_1.1-7        
##  [7] ConcaveFDR_0.0.1     R.utils_1.32.4       R.oo_1.18.0         
## [10] R.methodsS3_1.6.1    logcondens_2.1.1     fdrtool_1.2.12      
## [13] ggplot2_1.0.0        pastecs_1.3-18       boot_1.3-11         
## [16] fda_2.4.3            Matrix_1.1-4         reshape2_1.4        
## [19] knitr_1.6            BiocInstaller_1.14.2
## 
## loaded via a namespace (and not attached):
##  [1] BatchJobs_1.3       BBmisc_1.7          BiocGenerics_0.10.0
##  [4] brew_1.0-6          checkmate_1.3       codetools_0.2-9    
##  [7] colorspace_1.2-4    DBI_0.2-7           digest_0.6.4       
## [10] evaluate_0.5.5      fail_1.2            foreach_1.4.2      
## [13] formatR_0.10        gtable_0.1.2        iterators_1.0.7    
## [16] labeling_0.2        lattice_0.20-29     MASS_7.3-34        
## [19] munsell_0.4.2       parallel_3.2.0      proto_0.3-10       
## [22] Rcpp_0.11.2         RSQLite_0.11.4      scales_0.2.4       
## [25] sendmailR_1.1-2     stringr_0.6.2       tools_3.2.0
```

