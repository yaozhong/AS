# Author: Yao-zhong Zhang
# Date: 10/06/2014
# Implementation of weighted reads estimation

# Pipeline
# step 1. loading bam file and reference file
# step 2. foreach gene in the reference, estimate isoform expression 
# step 3. calculate basic information for estimation.
# step 4. get the EM estimate results
# step 5. save the estimation results

#######################################
#step 1
#######################################
setwd("./AS/weighted/")
library(GenomicFeatures)
library(Rsamtools)
library(parallel)

ref_file <- '/Users/Yaozhong/Research/2014_JP_HGC/Database/UCSC/human/genes.gtf'
bam_file1 <- "/Users/Yaozhong/Research/2014_JP_HGC/R.workspace/AS/data/K562/rep1.bam"
bam_file2 <- "/Users/Yaozhong/Research/2014_JP_HGC/R.workspace/AS/data/K562/rep2.bam"

source("./rna_util.R")

cat("Loading REF and BAM... ")
t_loading <- proc.time()

refData <- loadRef(ref_file)
bamData_1 <- loadBam(bam_file1)
bamData_2 <- loadBam(bam_file2)

t_load <- proc.time() - t_loading

cat("Done!")
cat(paste("\n-> using time", round(t_load[3]/60,2) , " min. \n"))

#######################################
#step 2
#foreach gene in the reference, estimate isoform expression 
#######################################
t_processing <- proc.time()

#human_chr_namelist <- c(paste("chr",seq(1,22), sep=""), "chrX", "chrY")
human_chr_namelist <- c("chr21")

res <- list()
res.w <- list()

for(chrName in human_chr_namelist){
  cat(paste("\n* Trnascript expression on ", chrName, " ..."))
  txEx <- processByChr(refData=refData, bam=bamData_1, chrName, readLength=75)
  res[[chrName]] <- txEx$"EM"
  res.w[[chrName]] <- txEx$"weighted"
}

t_proc <- proc.time() - t_processing
cat("Done!")
cat(paste("\n-> using time", round(t_proc[3]/60,2) , " min. \n"))


reads.chr.count <- unlist(lapply(res, function(x) { x$chrReadCount}))
reads.total <- sum(reads.chr.count)

dataDump <- function(res){
  
  for(chr in names(res)){
    res.chr <- res[[chr]]
    fpkm <- unlist(lapply(unname(res$chr$fpkmC), function(x, all) { fpkm <- x/all * 10^9 }, all=reads.total))
    
  }
  
}


