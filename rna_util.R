# RNA-seq processing utility

library(casper)
library(parallel)
source("./estimator.R")

# fragment 
FRAG.LEN.MEAN <- 187
FRAG.LEN.SD <- 50


loadingData <- function(ref, bamfile, genomeName="hg19"){
  
  # loading the reference
  print("Loading reference annotation ...")
  genDB <- makeTranscriptDbFromGFF(ref,format="gtf")
  #humanDB <- procGenome(genDB, genomeName, mc.cores=6)
  
  # loading Bam file
  print("Loading alignment BAM file ...")
  what <- scanBamWhat(); 
  what <- what[!(what %in% c('seq','qual'))]
  flag <- scanBamFlag(isPaired=T, hasUnmappedMate=F)
  param <- ScanBamParam(flag=flag, what=what, tag='XS')
  bam <- scanBam(bamfile, param=param)
  
  input <- list(humanDB, bam[[1]])
}

loadRef <- function(ref){
  
  # loading the reference
  cat(paste("\n* Loading Reference Genome Data [", ref, "]..."))
  genDB <- makeTranscriptDbFromGFF(ref,format="gtf")
  
  cat("Done")
  genDB  
  
}

loadBam <- function(bamfile){
  
  # loading Bam file
  cat(paste("\n* Loading BAM file [", bamfile, "]..."))
  what <- scanBamWhat(); 
  what <- what[!(what %in% c('seq','qual'))]
  flag <- scanBamFlag(isPaired=T, hasUnmappedMate=F)
  param <- ScanBamParam(flag=flag, what=what, tag='XS')
  bam <- scanBam(bamfile, param=param)
  
  cat("done")
  bam[[1]]
}

##### casper further step of annotation processing.
#humanDB <- procGenome(genDB, genomeName, mc.cores=6)

#### one chromesone processing
processByChr <- function(annoDB.chr, bam.chr, chrName, readLength, mc){
  
  txExList <- list()
  fpkmCountList <- list()
  txExList.w <- list()
  fpkmCountList.w <- list()
  r.cnt.chr <- 0
  
  # setting chromesome filter
#   seqlevels(refData, force=TRUE) <- c(chrName)
#   annoDB.chr <- procGenome(refData, chrName, mc.cores=6)
#   
#   bam0 <- rmShortInserts(bam, isizeMin=100)
#   bam.chr <- getBamReadsByChr(bam0, chrName)
  
  # be careful with pbam case
  #pbam.chr <- procBam(bam.chr)  # reads is plit accroding to CIGAR codes
  #frags <- getFragsFromPBam(pbam.chr)
  
  frags <- getFragPairFromBam(bam.chr, readLength)
  frags.left <- frags[["1"]]
  frags.right <- frags[["2"]]
  
  
  #for(idx in 1:length(annoDB.chr@transcripts)){
  
  txp <- function(idx){
    
    cat(paste("\n* processing", idx, "/", length(annoDB.chr@transcripts)," transcript"))
    
    ts <- annoDB.chr@transcripts[[idx]]
    ts.name <- names(ts)
    ts.num <- length(ts.name)
    
    # ignore the inference for the single isoform case
#     if(ts.num == 1){
#       #next
#       return(NULL)
#     }
    
    ts.len <- getTxLength(annoDB.chr, idx)
    
    exonGRs <- annoDB.chr@islands[[idx]]
    etMat <- getExTxMat(annoDB.chr, idx)
    
    # calculate exon weight
    ex.weight <- 1/apply(etMat,1,sum)
    ex.weight <- ex.weight/sum(ex.weight)

    exonGRs.start <- start(exonGRs)
    exonGRs.end <- end(exonGRs)
    
    tx.start <- min(exonGRs.start)
    tx.end   <- max(exonGRs.end)
    
    # filter fragments in the transcript
    sel.idx.left <-  which(start(frags.left) >= tx.start & (end(frags.left) <= tx.end) )
    sel.idx.right <-  which(start(frags.right) >= tx.start & (end(frags.right) <= tx.end) )
    sel.idx <- intersect(sel.idx.left, sel.idx.right)
    frags.tx.left <- frags.left[sel.idx]
    frags.tx.right <- frags.right[sel.idx]
    
    ## find overlap of frags and exons
    over.left <- findOverlaps(frags.tx.left, exonGRs)
    over.right <- findOverlaps(frags.tx.right, exonGRs)
    #cnt <- countOverlaps(frags.tx, exonGRs)
    
    # if no reads aligned to this reigon, no expresssion directly assign zero
    # next step is to estimate FPKM-value
    if(length(over.left) == 0 & length(over.right) == 0){
      #next
      return(NULL)
    }
     
    rid.set <- intersect(unique(queryHits(over.left)), unique(queryHits(over.right)))
    rid.num <- length(rid.set)

    ### reture value for single exon case
    if(ts.num == 1){
      fpkm.count <- rid.num/ts.len[1]
      tx.res <- list("EM"=c(1), "WEM"=c(1), "fpkmC"=fpkm.count, "fpkmC.wem"= fpkm.count, "readC"= rid.num, "tsName" = ts.name)
      return(tx.res)
    }

    ###prepare the process for the multi-exon case
    # maintanced results for the transcript
    rtMat <- matrix(nrow= rid.num, ncol=ts.num)
    ftMat.len <- matrix(nrow= rid.num, ncol=ts.num)
    r.weights <- numeric(rid.num)
    r.exPathList <- list()
    ri <- 1
    ridnames <- numeric(rid.num)
    
    
    
    #for(rid in rid.set){
    r.p <- function(rid){
       
       sel.idx <- which(queryHits(over.left) == rid)
       e.idx.left <- subjectHits(over.left)[sel.idx]
       e.ids.left <- names(exonGRs)[e.idx.left]
       
       sel.idx <- which(queryHits(over.right) == rid)
       e.idx.right <- subjectHits(over.right)[sel.idx]
       e.ids.right <- names(exonGRs)[e.idx.right]
       
       ## pay attention for the results, this 
       e.ids <- union(e.ids.left, e.ids.right)
       
       ## NOTE DONOT directly subset GRanges OBJECT!!! Very Very slow 
       cfrag.start <- start(frags.tx.left)[rid]
       cfrag.end <- end(frags.tx.right)[rid]
         
      
       # save the exon path information
       r.exPath <- e.ids
       r.exPathList[[as.character(rid)]] <<- r.exPath 
       
      ## generate compatible matrix  
      for(j in  1:ts.num ){
        res <- as.numeric(r.exPath) %in% ts[[j]]
        if(all(res)) { rtMat[ri, j] <<- 1 }
        else{ rtMat[ri, j] <<- 0 }
      } 
       
      # avoid invalid paried-end reads, next will over write current line
      if(sum(rtMat[ri,]) == 0){
        #next 
        return(NULL)
      }

      #ex.first.idx <- e.idx.left[1]
      #ex.last.idx  <- e.idx.right[length(e.idx.right)]
      #note the exonGRs is in reverse order
      ex.first.idx <- e.idx.left[length(e.idx.left)]
      ex.last.idx  <- e.idx.right[1] 
     
      start.ex.fg <- max(cfrag.start, exonGRs.start[ex.first.idx])
      end.ex.fg   <- min(cfrag.end, exonGRs.end[ex.last.idx])

      if(ex.first.idx != ex.last.idx){
        
        # check if there is any ambigiuous exons in between
        s.idx <- min(ex.first.idx, ex.last.idx)
        e.idx <- max(ex.first.idx, ex.last.idx)
        
        ## scan any exons in between and extend
        for(ei in (s.idx+1):(e.idx-1)){ 
          if(names(exonGRs)[ei] %in% e.ids == FALSE) e.ids <- c(e.ids, names(exonGRs)[ei])
        }
        
        # multi exon cases
        for(tj in 1:ts.num){
        
          if(rtMat[ri, tj] == 0){
            ftMat.len[ri,tj] <<- 0
          }else{
            
            cfrag.len <- (exonGRs.end[ex.first.idx] - start.ex.fg + 1) + (end.ex.fg - exonGRs.start[ex.last.idx] + 1)
            
            #if(abs(ex.first.idx - ex.last.idx) > 1){
            if(e.idx - s.idx > 1){
                for(ei in (s.idx+1):(e.idx-1)){
                  cfrag.len <- cfrag.len + width(exonGRs)[ei]*etMat[ei, tj]
                }
            }
            ftMat.len[ri,tj] <<- cfrag.len
          }
        }
      }else{
        # single exon case
        cfrag.len <- end.ex.fg - start.ex.fg + 1
        ftMat.len[ri,] <<- rep(cfrag.len, ts.num) * rtMat[ri,]
      }
      
      # read filter 2:  consider the p
      #fp <- dnorm(ftMat.len[ri,], FRAG.LEN.MEAN, FRAG.LEN.SD)
      #cond <- fp * rtMat[ri,]
      #if(sum(cond) == 0){
      #  return(NULL)
      #}
      
      r.weight <- max(ex.weight[e.ids])
      r.weights[ri] <<- r.weight
      
      ridnames[ri] <<- rid
      ri <<- ri + 1   
     } # loop of read process
    
    # note: shared varible can not be used with mclapply
    invisible(sapply(rid.set, r.p))
    # sapply(rid.set, r.p)

    # the reads is too few to be used set threshold here.
    if(ri == 1) return(NULL) #next
    
    # 

    #Trunck the unused
    ftMat.len <- matrix(ftMat.len[1:(ri-1),], nrow= ri-1)
    rtMat <- matrix(rtMat[1:(ri-1),], nrow = ri-1)
    r.weights <- r.weights[1:(ri-1)]
    ridnames <- ridnames[1:(ri-1)]
    
    colnames(rtMat) <- names(ts)
    
    #r.weights <- r.weights/sum(r.weights) 

    rownames(ftMat.len) <- ridnames
    rownames(rtMat) <- ridnames
    #colnames(ftMat.len) <- ts.name
    #colnames(rtMat) <- ts.name
    
    # do Parameters Estimation
    theta.w <- EM_estimator(r.weights, 1000, rtMat, ftMat.len, ts.len, TRUE)
    theta   <- EM_estimator(r.weights, 1000, rtMat, ftMat.len, ts.len)
    
    # calcualate FPKM value
    rtMat.th <- t(apply(rtMat, 1, function(x, t) {x*t}, t=theta) )
    fpkm.count <- apply(rtMat.th, 2, sum)
    fpkm.count <- fpkm.count/ts.len 

    rtMat.th.w <- t(apply(rtMat, 1, function(x, t) {x*t}, t=theta.w) )
    fpkm.count.w <- apply(rtMat.th.w, 2, sum)
    fpkm.count.w <- fpkm.count.w/ts.len 
     
    
#     for(k in 1:length(ts.len)){
#       txExList[[ts.name[k]]] <<- theta[k]
#       fpkmCountList[[ts.name[k]]] <<- fpkm.count[k]
#       
#       txExList.w[[ts.name[k]]] <<- theta.w[k]
#       fpkmCountList.w[[ts.name[k]]] <<- fpkm.count.w[k]    
#     }
    
    # real active reads
    r.cnt.tx <- nrow(rtMat.th)
    # r.cnt.chr <- r.cnt.chr + r.cnt.tx
    

    ## to do things
    #return necessary results for each transcript
    
    tx.res <- list("EM"=theta, "WEM"=theta.w, "fpkmC"=fpkm.count, "fpkmC.wem"= fpkm.count.w, "readC"= r.cnt.tx, "tsName" = ts.name)
    
  }# end loop for transcript processing
  
  txs.res <- mclapply(1:length(annoDB.chr@transcripts), txp, mc.cores = mc)
  
  # restore the default
  # refData <- restoreSeqlevels(refData)  #note this is necessary before next run
  
  ## res <- list("p"=txExList, "fpkmC"=fpkmCountList, "chrReadCount"=r.cnt.chr)
  ## res.w <- list("p"=txExList.w, "fpkmC"=fpkmCountList.w, "chrReadCount"=r.cnt.chr)
  ## out <- list("EM"=res, "weighted"=res.w)

  ## accemble the results for the transcript
}


#########################################
# filter the bam by chromesome name
#########################################
getBamReadsByChr <- function(bam, chr){
  
  sel <- which(bam$rname == chr)
  bam <- lapply(bam, function(z) z[sel])
}


## process exon related information
getExTxMat <- function(annoDB.chr, idx){
  
  eidVec <- names(annoDB.chr@islands[[idx]])
  txs <- annoDB.chr@transcripts[[idx]]
  
  num.e <- length(eidVec)
  num.t <- length(txs)
  
  etMat <- matrix(nrow=num.e, ncol=num.t)
  
  for(i in seq(1, num.e)){
    
    eid <- eidVec[i]
    for(j in seq(1, num.t)){
      
      if(eid %in% txs[[j]]){
        etMat[i,j] = 1
      }else{ etMat[i,j] = 0 }
    }   
  }
  
  rownames(etMat) <- eidVec
  colnames(etMat) <- names(txs)
  
  etMat
}

## get the paired-end reads fragments information
getFragsFromBam <- function(bam,readLength){
  
  # right reads process
  d <- bam$mpos - bam$pos
  sel.idx <- d<0
  n <- bam$qname[sel.idx] # name of reads
  sp <- bam$rname[sel.idx]
  names(sp) <- n
  
  en <- bam$pos[sel.idx] + readLength - 1
  names(en) <- n
  
  # left reads process
  sel.idx <- d>0
  st <- bam$pos[sel.idx]
  names(st) <- bam$qname[sel.idx]
  sel.idx <- match(n, names(st))
  
  st <- st[sel.idx[!is.na(sel.idx)]]
  en <- en[names(st)]
  sp <- sp[names(st)]
  
  sel.idx <- st<en
  st <- st[sel.idx]
  en <- en[sel.idx]
  sp <- sp[sel.idx]
  
  frags <- GRanges(IRanges(start=st, end=en), seqnames=sp)
}

getFragPairFromBam <- function(bam,readLength){
  
  # right reads process
  d <- bam$mpos - bam$pos
  sel.idx <- d<0
  n <- bam$qname[sel.idx] # name of reads
  sp <- bam$rname[sel.idx]
  names(sp) <- n
  
  en <- bam$pos[sel.idx] + readLength - 1
  names(en) <- n
  
  # left reads process
  sel.idx <- d>0
  st <- bam$pos[sel.idx]
  names(st) <- bam$qname[sel.idx]
  sel.idx <- match(n, names(st))
  
  st <- st[sel.idx[!is.na(sel.idx)]]
  en <- en[names(st)]
  sp <- sp[names(st)]
  
  sel.idx <- st<en
  st <- st[sel.idx]
  en <- en[sel.idx]
  sp <- sp[sel.idx]
  
  # frags <- GRanges(IRanges(start=st, end=en), seqnames=sp)
  frags.left <- GRanges(IRanges(start=st, end=st+readLength), seqnames=sp)
  frags.right <- GRanges(IRanges(start=en-readLength, end=en), seqnames=sp)
  
  output <- list("1"=frags.left, "2"=frags.right)
}

# to get fragment from processed Bam file
# by product of the exon path (main comparision with casper)
getFragsFromPBam <- function(pbam){
  
  # in pbam one reads could be splited into several parts, just keep the first one that is different from previous one.
  sel.idx <- c(TRUE, pbam@pbam$names[-1] != pbam@pbam$names[-length(pbam@pbam)] | (pbam@pbam$rid[-1] != pbam@pbam$rid[-length(pbam@pbam)]) ) 
  pbam@pbam <- pbam@pbam[sel.idx,]

  # select the left one, left one is indexed with 1, right one is index with 2
  sel.idx <- pbam@pbam$rid==1
  # seqname is the chromesome name
  frags <- GRanges(IRanges(start(pbam@pbam)[sel.idx], end(pbam@pbam)[pbam@pbam$rid==2]), seqnames=seqnames(pbam@pbam)[sel.idx])
  
}

###############################################
# get the transcript length for the idx-th transcripts
###############################################
getTxLength <- function(annoDB.chr,idx){
  
  ts <- annoDB.chr@transcripts[[idx]]
  
  ts.len <- numeric(length(ts))
  
  for(i in 1:length(ts)) {
    exonGRs <- annoDB.chr@exonsNI[ ts[[i]] ]
    ts.len[i] <- sum(width(exonGRs))
  }
  
  ts.len
}


## Develop LOG
# 2014/10/28 revised with parallel processing power
