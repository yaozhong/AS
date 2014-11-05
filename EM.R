# Author: Yao-zhong Zhang
# Date: 10/06/2014
# Description: Parameters Estimation

EM_estimator <- function(r.weights, maxIter, rtMat, ftMat.len, ts.len, RD=FALSE, thre=1e-4, verbose=FALSE){
  
  t.num <- length(ts.len)
  theta.0 <- rep(1/t.num, t.num)
  
  theta.c <- theta.0
  theta.p <- theta.0
  
  #r.weights <- r.weights/sum(r.weights)
  
  for(i in seq(1,maxIter)){
    
    if(verbose){
      cat(paste("\n EM iteration:[", i, "] ... "))
    }
    
    # E-step ** Note P(fragment is needed to be considerd)
    pt.mat <- t(apply(ftMat.len, 1, function(x) { dnorm(x,FRAG.LEN.MEAN,FRAG.LEN.SD)*(theta.c)/(ts.len-x+1) }))
    pt.mat <- pt.mat * rtMat
    
    pt.mat <- t(apply(pt.mat, 1, function(x) {x/sum(x)}))
    
    # M-step
    if(RD){
      theta <- apply(pt.mat, 2, function(x){ sum(r.weights*x, na.rm=TRUE)})
    }else{ 
      theta <- apply(pt.mat, 2, function(x) {sum(x,na.rm=TRUE)})
    }
    
    
    theta <- theta/sum(theta)
    if(any(is.na(theta))) {return(NULL)}
    
    # Update
    theta.p <- theta.c
    theta.c <- theta
  
    if(verbose){
    cat("Done. ")
    cat(theta)
    }
    
    if ( all( abs(theta.p - theta.c) <= thre) ) { break }
    
  }
  if(verbose){
  cat("\n")
  }
  
  theta.c
}





