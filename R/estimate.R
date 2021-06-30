estimate <- function(data, tol, p_AR, q_MA, covars=NULL, misReport="U")
{
  addintercept <- FALSE
  cv    <- matrix(rep(1,length(data)), ncol=1)
  num_rows <- 1
  if (!is.null(covars))
  {
    cv <- covars
    addintercept <- TRUE
    num_rows <- (ncol(cv)+1)
  }
  w0    <- 0.5
  q0    <- 0.5
  sd_x  <- sqrt(var(data)/(w0+(q0^2)*(1-w0)))
  if (!is.null(covars))
  {
  invisible(capture.output(res1 <- try(regmixEM(data, cv, addintercept, k=2, lambda=c(w0, (1-w0)), 
                                                beta=matrix(c(mean(data), mean(data)/(w0+q0*(1-w0))), ncol=2, nrow=num_rows),
                                                sigma=c(sd(data), sd_x), epsilon=1e-16, maxit=50000), silent=TRUE)))
  }else{
    invisible(capture.output(res1 <- try(normalmixEM(data, lambda=c(w0, (1-w0)), 
                                                  mu=c(mean(data), mean(data)/(w0+q0*(1-w0))),
                                                  sigma=c(sd(data), sd_x), epsilon=1e-16, maxit=50000), silent=TRUE))) 
  }
  if (class(res1)=="try-error") return(NA)
  if(!is.null(covars))
  {
    if(misReport=="U")
    {
      q0   <- min(res1$beta[1,])/max(res1$beta[1,])
      w0   <- res1$lambda[which(res1$beta[1,]==min(res1$beta[1,]))]
      k1   <- which(res1$beta[1,]==min(res1$beta[1,]))
      k2   <- which(res1$beta[1,]==max(res1$beta[1,]))
    }
    if(misReport=="O")
    {
      q0   <- max(res1$beta[1,])/min(res1$beta[1,])
      w0   <- res1$lambda[which(res1$beta[1,]==max(res1$beta[1,]))]
      k1   <- which(res1$beta[1,]==max(res1$beta[1,]))
      k2   <- which(res1$beta[1,]==min(res1$beta[1,]))
    }
  }
  if(is.null(covars))
  {
    if(misReport=="U")
    {
      q0   <- min(res1$mu)/max(res1$mu)
      w0   <- res1$lambda[which(res1$mu==min(res1$mu))]
      k1   <- which(res1$mu==min(res1$mu))
      k2   <- which(res1$mu==max(res1$mu))
    }
    if(misReport=="O")
    {
      q0   <- max(res1$mu)/min(res1$mu)
      w0   <- res1$lambda[which(res1$mu==max(res1$mu))]
      k1   <- which(res1$mu==max(res1$mu))
      k2   <- which(res1$mu==min(res1$mu))
    }
  }
  ### Step 1: Reconstrueixo x amb els nous valors de q i w
  q   <- q0
  w   <- w0
  eps <- tol + 1
  eps_ant <- tol + 2
  res <- res1
  i=1
  while (eps >= tol & eps < eps_ant)
  {
    eps_ant <- eps
    aux   <- ifelse(res$posterior[, k2]>0.5, 0, 1)
    x1    <- ifelse(aux==0, data, NA)   ### Sèrie dels valors que no estan infrareportats
    x2    <- ifelse(aux==1, data, NA)   ### Sèrie dels valors que estan infrareportats
    q_ant <- q
    w_ant <- w
    invisible(capture.output(mod_infra   <- try(arima(x2, order=c(p_AR, 0, q_MA), xreg=covars, method="ML"), silent=TRUE)))
    invisible(capture.output(mod_noinfra <- try(arima(x1, order=c(p_AR, 0, q_MA), xreg=covars, method="ML"), silent=TRUE)))
    i=1
    while ((class(mod_infra)=="try-error" | class(mod_noinfra)=="try-error") & i < 10)
    {
      sd_x <- sqrt(var(data)/(w+(q^2)*(1-w)))
      if(!is.null(covars))
      {
      invisible(capture.output(res <- try(regmixEM(data, cv, addintercept, k=2, lambda=c(w, (1-w)), 
                                                   beta=matrix(c(mean(data), mean(data)/(w+q*(1-w))), ncol=2, nrow=num_rows),
                                                   sigma=c(sd(data), sd_x), epsilon=1e-16,maxit=50000), silent=TRUE)))
      }else{
        invisible(capture.output(res1 <- try(normalmixEM(data, lambda=c(w, (1-w)), 
                                                         mu=c(mean(data), mean(data)/(w+q*(1-w))),
                                                         sigma=c(sd(data), sd_x), epsilon=1e-16, maxit=50000), silent=TRUE)))
      }
      if (!is.null(covars))
      {
        if(misReport=="U")
        {
          q <- min(res$beta[1,])/max(res$beta[1,])
          w <- res$lambda[which(res$beta[1,]==min(res$beta[1,]))]
          k1    <- which(res$beta[1,]==min(res$beta[1,]))
          k2    <- which(res$beta[1,]==max(res$beta[1,]))
        }
        if(misReport=="O")
        {
          q <- max(res$beta[1,])/min(res$beta[1,])
          w <- res$lambda[which(res$beta[1,]==max(res$beta[1,]))]
          k1    <- which(res$beta[1,]==max(res$beta[1,]))
          k2    <- which(res$beta[1,]==min(res$beta[1,]))
        }
      }
      if (is.null(covars))
      {
          if(misReport=="U")
          {
            q <- min(res$mu)/max(res$mu)
            w <- res$lambda[which(res$mu==min(res$mu))]
            k1    <- which(res$mu==min(res$mu))
            k2    <- which(res$mu==max(res$mu))
          }
          if(misReport=="O")
          {
            q <- max(res$mu)/min(res$mu)
            w <- res$lambda[which(res$mu==max(res$mu))]
            k1    <- which(res$mu==max(res$mu))
            k2    <- which(res$mu==min(res$mu))
          }
      }
      aux   <- ifelse(res$posterior[, k2]>0.5, 0, 1)
      x1    <- ifelse(aux==0, data, NA)   ### Sèrie dels valors que no estan infrareportats
      x2    <- ifelse(aux==1, data, NA)   ### Sèrie dels valors que estan infrareportats
      invisible(capture.output(mod_infra   <- try(arima(x2, order=c(p_AR, 0, q_MA), xreg=covars), silent=TRUE)))
      invisible(capture.output(mod_noinfra <- try(arima(x1, order=c(p_AR, 0, q_MA), xreg=covars), silent=TRUE)))
      i<-i+1
    }
    if (class(mod_infra)!="try-error" & class(mod_noinfra)!="try-error")
    {
      q     <- as.numeric(coef(mod_infra)[names(coef(mod_infra))=="intercept"]/coef(mod_noinfra)[names(coef(mod_noinfra))=="intercept"])
      
      if(!is.null(covars))
      {
        invisible(capture.output(res <- try(regmixEM(data, cv, addintercept, k=2, lambda=c(w, (1-w)), 
                                                     beta=matrix(c(coef(mod_infra)[names(coef(mod_infra))=="intercept"],
                                                   coef(mod_noinfra)[names(coef(mod_noinfra))=="intercept"]), ncol=2, nrow=num_rows),
                                                     sigma=c(sd(x2, na.rm=T), sd(x1, na.rm=T)), epsilon=1e-16,maxit=50000), silent=TRUE)))
      }else{
        invisible(capture.output(res <- try(normalmixEM(data, lambda=c(w, (1-w)), 
                                                         mu=c(coef(mod_infra)[names(coef(mod_infra))=="intercept"],
                                                              coef(mod_noinfra)[names(coef(mod_noinfra))=="intercept"]),
                                                         sigma=c(sd(x2, na.rm=T), sd(x1, na.rm=T)), epsilon=1e-16, maxit=50000), silent=TRUE)))
      }
      
      
      # invisible(capture.output(res <- normalmixEM(data, k=2, mean.constr=c(coef(mod_infra)[names(coef(mod_infra))=="intercept"],
      #                                                                      coef(mod_noinfra)[names(coef(mod_noinfra))=="intercept"]), 
      #                                             sd.constr=c(sd(x2, na.rm=T), sd(x1, na.rm=T)), 
      #                                             maxit=50000, epsilon=1e-16, maxrestarts=1000)))
      #w <- sum(!is.na(x2))/length(data)
      if (!is.null(covars))
      {
        if(misReport=="U")
        {
          k1    <- which(res$beta[1,]==min(res$beta[1,]))
          k2    <- which(res$beta[1,]==max(res$beta[1,]))
        }
        if(misReport=="O")
        {
          k1    <- which(res$beta[1,]==max(res$beta[1,]))
          k2    <- which(res$beta[1,]==min(res$beta[1,]))
        }
      }
      if (is.null(covars))
      {
        if(misReport=="U")
        {
          k1    <- which(res$mu==min(res$mu))
          k2    <- which(res$mu==max(res$mu))
        }
        if(misReport=="O")
        {
          k1    <- which(res$mu==max(res$mu))
          k2    <- which(res$mu==min(res$mu))
        }
      }
      w     <- res$lambda[k1]
      eps   <- sqrt((q_ant-q)^2+(w_ant-w)^2)
    }else{
      return(NA)
    }
  }
  ### Step 2: Reconstrueixo x i estimo alpha, theta i sigma2
  x_rec   <- ifelse(aux==0, data, data/q)
  mod_rec <- arima(x_rec, order=c(p_AR, 0, q_MA), xreg=covars, method="ML")
  if (p_AR==0)
  {
    alphas <- 0
  }else{
    alphas <- as.numeric(coef(mod_rec)[1:p_AR])
  }
  if (q_MA==0)
  {
    thetas <- 0
  }else{
    thetas <- as.numeric(coef(mod_rec)[(p_AR+1):(p_AR+q_MA)])
  }
  sigma2  <- mod_rec$sigma2
  const   <- as.numeric(coef(mod_rec)[names(coef(mod_rec))=="intercept"])
  estimates <- c(const, alphas, thetas, sigma2, q, w, AIC(mod_rec))
  n_alphas <- vector()
  n_thetas <- vector()
  for (i in 1:p_AR)
  {
    n_alphas[i] <- paste0("alpha", i)
  }
  for (i in 1:q_MA)
  {
    n_thetas[i] <- paste0("theta", i)
  }
  names(estimates) <- c("const", n_alphas, n_thetas, "var", "q", "w", "AIC")
  return(estimates)
}
