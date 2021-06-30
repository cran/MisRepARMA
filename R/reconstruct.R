reconstruct <- function(data, q, w)
{
  # cv    <- matrix(rep(1,length(data)),ncol=1)
  # sd_x  <- sqrt(var(data)/(w+(q^2)*(1-w)))
  # invisible(capture.output(res <- try(regmixEM(data, cv, addintercept=FALSE, k=2, lambda=c(w, (1-w)), 
  #                                             beta=matrix(c(mean(data), mean(data)/(w+q*(1-w))), ncol=2, nrow=1),
  #                                             sigma=c(sd(data), sd_x), epsilon=1e-16, maxit=50000), silent=TRUE)))
  # k1    <- which(res$beta==min(res$beta))
  # k2    <- which(res$beta==max(res$beta))
  # aux   <- ifelse(res$posterior[, k2]>0.5, 0, 1)
  # 
  # x_rec   <- ifelse(aux==0, data, data/q)
  sd_x  <- sqrt(var(data)/(w+(q^2)*(1-w)))
  invisible(capture.output(res <- try(normalmixEM(data, lambda=c(w, (1-w)), 
                                              mu=c(mean(data), mean(data)/(w+q*(1-w))),
                                              sigma=c(sd(data), sd_x), epsilon=1e-16, maxit=50000), silent=TRUE)))
  k1    <- which(res$mu==min(res$mu))
  k2    <- which(res$mu==max(res$mu))
  aux   <- ifelse(res$posterior[, k2]>0.5, 0, 1)
  
  x_rec   <- ifelse(aux==0, data, data/q)
  return(x_rec)
}
