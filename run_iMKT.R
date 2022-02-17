#!/home/rs2474/R/bin/Rscript

## run iMKT analysis using R package

library(iMKT)
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("Three arguments must be supplied (daf file, divergence file).n", call.=FALSE)
} else {
  # default output file
  path = '.'
  dafData = args[1]
  divData = args[2]
}

setwd(path)
# getwd()

my_dafData = read.delim(dafData)
my_divData = read.delim(divData)

# Perform standard MKT
#cat("Standard MKT: \n")
#standardMKT(my_dafData, my_divData)

## Perform FWW correction
cat("FWW correction MKT: \n")
# FWW(myDafData, myDivergenceData)
FWW(my_dafData, my_divData, listCutoffs=c(0, 0.05, 0.1))

## Perform asymptotic MKT, failed
# cat("Asymptotic MKT: \n")
# asymptoticMKT(myDafData, myDivergenceData, xlow=0, xhigh=0.9)

## Perform asymptotic MKT
# cat("iMKT: \n")
# iMKT(myDafData, myDivergenceData, xlow=0, xhigh=0.9)

## Perform complete MKT methodologies
## Failed... for FWW correction
# cat("Complete MKTs (standard, FWW, DGRP, asymptotic, iMKT): \n")
# completeMKT(myDafData, myDivergenceData, xlow=0, xhigh=0.9)
