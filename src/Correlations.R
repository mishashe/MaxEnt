

#==============================================================================================================================================================
# KCL has been using GRCh19
setwd("/home/m.sheinman/Development/precision-Clonality")
require(doMC)
registerDoMC(36)
library(data.table)
library(readr)
library(GenomicRanges)
library(ggplot2)
library(vcfR)
library(gtools)
library(stringr)
library(discover)
library(dplyr)
library(devtools)
library(stats)
library(devtools)
library(roxygen2)
library(matrixStats)
library(MASS)
library(pracma)
library(maxLik)
source("src/models/ClonalityPackageCorrelations/Ising.R")
N <- 10
nSamples <- 10000
h <- rnorm(N)-0
J <- matrix(rnorm(N*N),nrow=N, ncol=N)/1
J <- (J+t(J))/2
diag(J) <- 0

s <- sampleIsing(h, J, nSamples, 100)

m <- colMeans(mat)
c <- cov(mat);  diag(c) <- 0
c[is.na(c)] <- 0


source("src/models/ClonalityPackageCorrelations/Ising.R")
fit <- fitIsing(s, method=c("NoCor", "nMFT", "TAP")[2])
h_nMFT <- fit[[1]]
J_nMFT <- fit[[2]]

fit <- fitLikelihood(s)
h_LH <- fit[[1]]
J_LH <- fit[[2]]


pdf("src/models/ClonalityPackageCorrelations/plots/h.pdf")
plot(h,h_LH,col="blue", main=paste0(cor(h,h_LH,)," ",cor(h,h_nMFT,)))
points(h,h_nMFT,pch=19,col="red")
lines(h,h)
dev.off()

a <- Reshape(J, prod(dim(J)), 1)
b_LH <- Reshape(J_LH, prod(dim(J_LH)), 1)
b_nMFT <- Reshape(J_nMFT, prod(dim(J_nMFT)), 1)
pdf("src/models/ClonalityPackageCorrelations/plots/J.pdf")
plot(a,b_LH,col="blue", main=paste0(cor(a,b_LH,)," ",cor(a,b_nMFT,)))
points(a,b_nMFT,pch=19,col="red")
lines(a,a)
dev.off()








