#This script was used to perform rma normalization

# install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")  # Get data from NCBI Gene Expression Omnibus (GEO)
BiocManager::install("affy", force=TRUE)
BiocManager::install("hgu133a.db")
BiocManager::install("hgu133acdf")

library(affy)
library(hgu133a.db)
library(hgu133acdf)
library(GEOquery)

#data source https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5462
untar("GSE5462_RAW.tar", exdir = "GSE5462")
cels = list.files("GSE5462/", pattern = "CEL")
sapply(paste("GSE5462", cels, sep = "/"), gunzip)
cels = list.files("GSE5462/", pattern = "CEL")
setwd("GSE5462")
raw.data = ReadAffy(verbose = FALSE, filenames = cels, cdfname = "hgu133acdf")
# perform RMA normalization (log2)
data.rma.norm = rma(raw.data,bgversion=1)
rma = exprs(data.rma.norm)