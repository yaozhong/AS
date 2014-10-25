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


main()

main <- function(){

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

#human_chr_namelist <- c(paste("chr",seq(1,22), sep=""), "chrX", "chrY")
human_chr_namelist <- c("chr21")

res <- list()
res.w <- list()
bam <- bamData_1

for(chrName in human_chr_namelist){
  
  cat(paste("\n* Trnascript expression on ", chrName, " ..."))
  
  # process the reference file 
  seqlevels(refData, force=TRUE) <- c(chrName)
  annoDB.chr <- procGenome(refData, chrName, mc.cores=4)
  bam0 <- rmShortInserts(bam, isizeMin=100)
  bam.chr <- getBamReadsByChr(bam0, chrName)
  
  # pre-process for chromesome
  t_processing <- proc.time()
  txEx <- processByChr(annoDB.chr=annoDB.chr, bam.chr=bam.chr, chrName, readLength=75, mc=1)
  t_proc <- proc.time() - t_processing
  cat(paste("\n-> using time", round(t_proc[3]/60,2) , " min. \n"))
  # reset
  refData <- restoreSeqlevels(refData) 
  
  
  #res[[chrName]] <- txEx$"EM"
  #res.w[[chrName]] <- txEx$"weighted"
}




#reads.chr.count <- unlist(lapply(res, function(x) { x$chrReadCount}))
#reads.total <- sum(reads.chr.count)

}

dataDump <- function(res){
  
  for(chr in names(res)){
    res.chr <- res[[chr]]
    fpkm <- unlist(lapply(unname(res$chr$fpkmC), function(x, all) { fpkm <- x/all * 10^9 }, all=reads.total))
    
  }
  
}


