####################################################
# Author: Yao-zhong Zhang
# Date: 10/31/2014
# Refined version for the weighted EM algorithm for RNA-seq estimation
####################################################

# Part 1: Data Split and prepare data chunks for allocated machines
# Part 2: Raw Data Processing and prepare R-T matrix for EM algorithms
# Part 3: EM iteration
# Part 4: Results and analysis module

## <Part 1>

# Env Setting
setwd("/Users/Yaozhong/Research/2014_JP_HGC/R.workspace/AS/weighted/")
library(GenomicFeatures)
library(casper)
library(Rsamtools)
library(parallel)

source("./data_util.R")
source("./read_processing.R")


main <- function(){

	ref_file  <- "/Users/Yaozhong/Research/2014_JP_HGC/Database/UCSC/human/genes.gtf"
	bam_file1 <- "/Users/Yaozhong/Research/2014_JP_HGC/R.workspace/AS/data/K562/rep1.bam"
	bam_file2 <- "/Users/Yaozhong/Research/2014_JP_HGC/R.workspace/AS/data/K562/rep2.bam"
	bam_files <- c(bam_file1)
	human_chr_namelist <- c(paste("chr",seq(1,22), sep=""), "chrX", "chrY")
  #human_chr_namelist <- c("chr21", "chr22")

  # loading the data.
	data <- data_load_split(bam_files, ref_file, human_chr_namelist)
  
  # inference processing
  output <- list()
  for(chrName in human_chr_namelist){
    cat("\n Processing ", chrName, "...");
    ress <- processByChr(tsDB=data$REF[[chrName]], bam.chr=data$BAMS[[1]][[chrName]], readLength=75, mc=6);
    output[[chrName]] <- ress
  }
  
  return(output)
  
} 


# loading reference Genome and Bam files
data_load_split <- function(bamFiles, ref_file, human_chr_namelist, mc=6){

	cat("> Loading [Genome Reference] and [Aligned BAM Files] ... ")

	cat("\n >> Reference Genome loading ... ")
  refData <- loadRef(ref_file)
  cat("DONE!")
	

	cat("\n >> Reference Genome spliting into chromesomes ... ")
	refData.byChr <- list()
	for(chrName in human_chr_namelist){
		  seqlevels(refData, force=TRUE) <- c(chrName)
  		annoDB.chr <- procGenome(refData, chrName, mc.cores=mc)
  		refData <- restoreSeqlevels(refData)
  		refData.byChr[[chrName]] <- annoDB.chr
	}
  cat("DONE!")


  bamDataList.byChr <- list()
	for(idx in 1:length(bamFiles)){
	  bamFile <- bamFiles[idx]  
		bamData <- loadBam(bamFile)
		cat("\n >> Bam file [", bamFile, "] loaded!")

		bam.byChr <- list()
		bam0 <- rmShortInserts(bamData, isizeMin=100)

		cat("\n >> Bam file [", bamFile, "] spliting into chromesome ...")
		for(chrName in human_chr_namelist){

  			bam.chr <- getBamReadsByChr(bam0, chrName)
  			bam.byChr[[chrName]] <- bam.chr
		}

		bamDataList.byChr[[idx]] <- bam.byChr
    
		cat("DONE!\n")
		
	}

	data <- list("REF"=refData.byChr, "BAMS"=bamDataList.byChr)
}




