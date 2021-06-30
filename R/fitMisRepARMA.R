### Estimació amb interval de confiança (bootstrap paramètric) per tots els paràmetres
fitMisRepARMA <- function(y, tol, B, alpha, p_AR, q_MA, covars=NULL, misReport="U")
{
  params_post <- NA
  class(params_post) <- "try-error"
  j <- 1
  while(class(params_post)=="try-error" & j < 10)
  {
    invisible(capture.output(params_post <- try(tsbootstrap(y, nb=B, statistic=estimate, 
              type="stationary", tol=tol, p_AR=p_AR, q_MA=q_MA, covars=covars, misReport=misReport), silent=TRUE)))
    j <- j + 1
  }
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
  colnames(params_post$statistic) <- c("const", n_alphas, n_thetas, "var", "q", "w", "AIC")
  return(params_post) 
}
