####################################################
# Author: Yao-zhong Zhang
# Date: 10/31/2014
# Refined version for the weighted EM algorithm for RNA-seq estimation
####################################################

# Part 2: Raw Data Processing and prepare R-T matrix for EM algorithms

source("./EM.R")

# fragment length distribution 
FRAG.LEN.MEAN <- 187
FRAG.LEN.SD <- 50

# general basic process unit
processByChr <- function(ref.chr, bam.chr, chrName, readLength, mc){

  # Process the paried read alignemtn
  frags <- getFragPairFromBam(bam.chr, readLength)
  frags.left <- frags[["1"]]
  frags.right <- frags[["2"]]
  tsDB <- ref.chr

  z <- function(idx){
  	res <- processByGene(tsDB, idx, readLength)
  }

#   t_processing <- proc.time()
#   for( idx in 1:length(tsDB@transcripts)){
#   	cat("\nidx-", idx)
#   	res <- processByGene(tsDB, idx, readLength)
# 
#   }
#   t_proc <- proc.time() - t_processing
#   cat(paste("\n-> using time", round(t_proc[3]/60,2) , " min. \n"))

  t_processing <- proc.time()
  res <- mclapply(1:length(tsDB@transcripts), z, mc.cores=6)
  t_proc <- proc.time() - t_processing
  cat(paste("\n-> using time", round(t_proc[3]/60,2) , " min. \n"))
  

}


processByGene <- function(tsDB, idx, readLength=75){

	  ts <- tsDB@transcripts[[idx]]
	  ts.name <- names(ts)
    ts.num <- length(ts.name)

    exonGRs <- tsDB@islands[[idx]]
    etMat <- getExTxMat(tsDB, idx)

    ex.weight <- 1/apply(etMat,1,sum)
    ex.weight <- ex.weight/sum(ex.weight)
    
	  exonGRs.start <- start(exonGRs)
	  exonGRs.end <- end(exonGRs)
	  exonGRs.width <- width(exonGRs)
	  names(exonGRs.start) <- names(exonGRs)
	  names(exonGRs.end) <- names(exonGRs)
	  names(exonGRs.width) <- names(exonGRs)
    
	  ts.len <- sapply(ts, function(x) {sum(exonGRs.width[as.character(x)])})
	  tx.start <- min(exonGRs.start)
	  tx.end   <- max(exonGRs.end)

    # filter fragments in the transcript, torrent partial matching  ** revised one
	  sel.idx.left <-  which(start(frags.left) >= tx.start -readLength & (end(frags.left) <= tx.end + readLength) )
	  sel.idx.right <-  which(start(frags.right) >= tx.start -readLength & (end(frags.right) <= tx.end + readLength) )	    
    sel.idx <- intersect(sel.idx.left, sel.idx.right)
    frags.tx.left <- frags.left[sel.idx]
    frags.tx.right <- frags.right[sel.idx]

    ## find overlap of frags and exons
    over.left <- findOverlaps(frags.tx.left, exonGRs)
    over.right <- findOverlaps(frags.tx.right, exonGRs)

    # none of the alignment
    if(length(over.left) == 0 & length(over.right) == 0){
      #next
      return(NULL)
    }

    # get the exon region aligned paired-end readds
    rid.set <- intersect(unique(queryHits(over.left)), unique(queryHits(over.right)))
    rid.num <- length(rid.set)

    if(ts.num == 1){
      fpkm.count <- rid.num/ts.len[1]
      tx.res <- list("EM"=c(1), "WEM"=c(1), "fpkmC"=fpkm.count, "fpkmC.wem"= fpkm.count, "readC"= rid.num)
      return(tx.res)
    }

    ###############################################################
    # read process
    getRead2Exons <- function(rid){

       sel.idx <- which(queryHits(over.left) == rid)
       e.idx.left <- subjectHits(over.left)[sel.idx]
       e.ids.left <- names(exonGRs)[e.idx.left]
       
       sel.idx <- which(queryHits(over.right) == rid)
       e.idx.right <- subjectHits(over.right)[sel.idx]
       e.ids.right <- names(exonGRs)[e.idx.right]

       e.ids <- sort(unique(c(e.ids.left, e.ids.right)))
    }

    # matrix generation
    r2exs <- lapply(rid.set, getRead2Exons)
    rtMat <- getReadTxMat(r2exs, ts)
    ###############################################################

    
    ###############################################################
    # key problems here: how to calcualte fragment length based on reads 

    ## use local read index for this
    cfrag.start <- start(frags.tx.left)[rid.set]
    cfrag.end <- end(frags.tx.right)[rid.set]


    getFragLenInTx <- function(r.idx){

    	e.ids <- r2exs[[r.idx]]

    	start.ex.fg <- max(cfrag.start[r.idx], exonGRs.start[e.ids[1]])
        end.ex.fg   <- min(cfrag.end[r.idx], exonGRs.end[e.ids[length(e.ids)]])

        cutB <- (start.ex.fg - exonGRs.start[e.ids[1]]) + (exonGRs.end[e.ids[length(e.ids)]] - end.ex.fg)
        cfrag.len <- sum(exonGRs.width[e.ids]) - cutB

        fit <- rtMat[r.idx,]*cfrag.len

        # consider the case
        p.ids <- intersect(as.character(e.ids[1]:e.ids[length(e.ids)]), names(exonGRs))

        if(length(p.ids) != length(e.ids)){
        	a.ids <- setdiff(p.ids, e.ids)

        	a.sum <- etMat[a.ids,]* exonGRs.width[a.ids]
        	if(length(a.ids) > 1){
        	  a.sum <- apply(etMat[a.ids,]* exonGRs.width[a.ids], 2, sum)
       		}
        	fit <- fit + a.sum*rtMat[r.idx,]
        }
        fit
    }

    ftMat <- t(matrix(sapply(1:rid.num, getFragLenInTx), nrow=ts.num))

    ###############################################################


    ###############################################################
    getReadWeights <- function(r.idx){

    	e.ids <- r2exs[[r.idx]]
    	p.ids <- intersect(as.character(e.ids[1]:e.ids[length(e.ids)]), names(exonGRs))

    	r.weight <- max(ex.weight[p.ids])
    }

    r.weights <- sapply(1:rid.num, getReadWeights)

    ##############################################################

    # EM-algorithm
	  ts.len <- unname(ts.len)
    theta.w <- EM_estimator(r.weights, 1000, rtMat, ftMat, ts.len, TRUE)
    theta   <- EM_estimator(r.weights, 1000, rtMat, ftMat, ts.len)


    # calcualate FPKM value
    rtMat.th <- t(apply(rtMat, 1, function(x, t) {x*t}, t=theta) )
    fpkm.count <- apply(rtMat.th, 2, sum)
    fpkm.count <- fpkm.count/ts.len 

    rtMat.th.w <- t(apply(rtMat, 1, function(x, t) {x*t}, t=theta.w) )
    fpkm.count.w <- apply(rtMat.th.w, 2, sum)
    fpkm.count.w <- fpkm.count.w/ts.len 

    tx.res <- list("EM"=theta, "WEM"=theta.w, "fpkmC"=fpkm.count, "fpkmC.wem"= fpkm.count.w, "readC"= length(r.weights))

}


getFragLenTxMat <- function(r2exs, txs){

	ex.first.idx <- e.idx.left[1]
    ex.last.idx  <- e.idx.right[length(e.idx.right)]
}


getReadTxMat <- function(r2exs, txs){
  res <- lapply(r2exs, function(exs) { rtCompt <- unlist(lapply(txs, function(t, ep) { all(as.numeric(ep) %in% t)+0} , ep=exs))} )
  rtMat <- t(matrix(unlist(res), nrow=length(txs)))
}


getExTxMat <- function(tsDB, idx){
  # note to consider consistency of order
  eidVec <- sort(as.numeric(names(tsDB@islands[[idx]])))
  txs <- tsDB@transcripts[[idx]]

  num.e <- length(eidVec)
  num.t <- length(txs)

  etMat <- matrix(unlist(lapply(txs, function(x) {eidVec %in% x})), nrow=num.e, ncol=num.t)+0
  
  rownames(etMat) <- eidVec
  colnames(etMat) <- names(txs)
  
  etMat
}


getFragPairFromBam <- function(bam, readLength){
  
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
  

  output <- list("1"=frags.left, "2"=frags.right, "fl"=bam$isize[sel.idx])
}
