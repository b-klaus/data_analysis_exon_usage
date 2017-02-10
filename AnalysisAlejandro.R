## ---- echo=FALSE---------------------------------------------------------
print(date())

## ----setup---------------------------------------------------------------
library(knitr)
library(reshape2)
library(concaveFDR)
library(ggplot2)
library(grid)
library(locfdr)
library(fdrtool)
library(BiocParallel)
library(grid)
library(RColorBrewer)
library(psych)
library(DEXSeq)
library(httr)
library(xml2)
library(purrr)
library(stringr)

opts_chunk$set(fig.width=14, fig.height=10, cache=TRUE, 
error=FALSE, cache.lazy=FALSE)
options(digits = 2, scipen = 10)

multicoreParam <- MulticoreParam(workers = round(parallel::detectCores() * 0.75) )
multicoreParam
register(multicoreParam)


## ----get data Alejandro, cache=FALSE-------------------------------------
base_url <- "http://www-huber.embl.de/pub/DEUprimates_supplement/SupplementaryFile3/objects"
r <- GET(base_url)
x <- content(r, type = "text/html")
# html_structure(x)
potential_files  <- xml_text(xml_find_all(x, "//a"))
potential_files

idx_actual_files <-  str_detect(potential_files, pattern = ".[Rr][Dd]a")
files_to_download <- potential_files[idx_actual_files]

map(paste0(base_url, "/", files_to_download), ~ load(url(.x), envir = .GlobalEnv) )

## ----inspect data Alejandro----------------------------------------------
head( spPairs )

head( sppCovs, 20)
str( sppCovs )
head( sppCors, 20)
str( sppCors )




## ----get siginifcant correlations----------------------------------------
head( padjCov, 20)
str( padjCov )
table(na.omit(padjCov) < 0.1)
# prop.table(table(na.omit(padjCov) < 0.1))


## ----run FDR analyses corrs----------------------------------------------
test.stats <- as.vector(na.omit(sppCors))
qplot(test.stats)


sub.sample <- round(length(test.stats)*0.7)
sub.sample

res.fdrtool <- fdrtool(sample(test.stats, sub.sample), statistic = "correlation")
sum(res.fdrtool$lfdr < 0.2) / sub.sample
sum(res.fdrtool$lfdr < 0.5) / sub.sample
sum(res.fdrtool$qval < 0.2) / sub.sample


#res.concaveFDR <- concaveFDR(sample(test.stats, sub.sample), statistic = "correlation")
#res.concaveFDR.non.sub <- concaveFDR(test.stats, statistic = "correlation")

## ----run FDR analyses z-transformed corrs--------------------------------
### z-value transform 
test.stats.z <- fisherz(test.stats)
test.stats.z <- test.stats.z - median(test.stats.z)
summary(test.stats.z)

test.stats.z <- test.stats.z[abs(test.stats.z) < 4]
summary(test.stats.z)
mad(test.stats.z)


sub.sample <- round(length(test.stats.z)*0.1)
sub.sample

res.fdrtool <- fdrtool(sample(test.stats.z, sub.sample), statistic = "normal")
sum(res.fdrtool$lfdr < 0.2) / sub.sample
sum(res.fdrtool$lfdr < 0.5) / sub.sample
sum(res.fdrtool$qval < 0.2) / sub.sample


### use transformed values

concaveCorr <- concaveFDR(sample(test.stats.z, sub.sample) , statistic = "normal", 
                            cutoff.method = "smoothing")

#concaveFDR:::smoothing.cutoff(test.stats.z, "normal")

sum(concaveCorr$lfdr.log < 0.2)  / sub.sample
sum(concaveCorr$lfdr.log < 0.5)  / sub.sample
sum(concaveCorr$qval.log < 0.2)  / sub.sample
sum(concaveCorr$qval.log < 0.1)  / sub.sample

## ----run FDR analyses covs-----------------------------------------------
### try on covariances on the natural log scale
test.stats.cov <- as.vector(na.omit(sppCovs))
test.stats.cov <- fisherz(test.stats.cov)
test.stats.cov <- na.omit(test.stats.cov)
test.stats.cov <- test.stats.cov - median(test.stats.cov)
test.stats.cov <- as.vector(test.stats.cov)
##save(test.stats.cov, file ="AlejandroCovs.rda")
#load("AlejandroCovs.rda")

sd(test.stats.cov)
summary(test.stats.cov)

test.stats.cov <- test.stats.cov[abs(test.stats.cov) < 4]
summary(test.stats.z)
mad(test.stats.z)

sub.sample <- round(length(test.stats.cov)*0.1)
sub.sample

FDR.cov <- fdrtool(sample(test.stats.cov, sub.sample), 
                               statistic = "normal")
sum(FDR.cov$lfdr < 0.2) / sub.sample
sum(FDR.cov$lfdr < 0.5) / sub.sample
sum(FDR.cov$qval < 0.2) / sub.sample
sum(FDR.cov$qval < 0.1)/ sub.sample


## ----Concave FDR on covariances------------------------------------------

res.concaveFDR.cov  <-concaveFDR(sample(test.stats.cov, sub.sample),
                               statistic = "normal", cutoff = "smoothing")


sum(res.concaveFDR.cov$lfdr.log < 0.2) / sub.sample
sum(res.concaveFDR.cov$lfdr.log < 0.5) / sub.sample
sum(res.concaveFDR.cov$qval.log < 0.2) / sub.sample
sum(res.concaveFDR.cov$qval.log < 0.1) / sub.sample



## ----SessionInfo, echo=FALSE---------------------------------------------

sessionInfo()


