## basic raw data (bam file and reference processing files)

loadRef <- function(ref){
  
  # loading the reference
  genDB <- makeTranscriptDbFromGFF(ref,format="gtf")
  genDB  
  
}

loadBam <- function(bamfile){
  
  # loading Bam file
  what <- scanBamWhat(); 
  what <- what[!(what %in% c('seq','qual'))]
  flag <- scanBamFlag(isPaired=T, hasUnmappedMate=F)
  param <- ScanBamParam(flag=flag, what=what, tag='XS')
  bam <- scanBam(bamfile, param=param)
  
  bam[[1]]
}

getBamReadsByChr <- function(bam, chr){
  
  sel <- which(bam$rname == chr)
  bam <- lapply(bam, function(z) z[sel])

}
