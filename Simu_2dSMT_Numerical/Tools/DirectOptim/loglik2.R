# for compare

loglik2 <- function(const, corr, corrmat, corrmodel, data, covariate,
                   dimat, fixed, fname, grid, ident, model, namescorr,  namescov, namesnuis,
                   param, setup) {
  loglik <- -1e+08
  param <- c(param, fixed)
  paramcorr <- param[namescorr]
  nuisance <- param[namesnuis]
  betac <- param[namescov]
  stdata <- data - matrix( covariate %*% betac,ncol=dim(data)[2]) #nuisance["mean"]
  cc = .C(corrmat, corr = corr, as.integer(corrmodel), 
          as.double(nuisance), as.double(paramcorr), DUP = TRUE, 
          NAOK = TRUE)$corr
  corr <- cc
  if (corr[1] == -2) 
    return(loglik)
  cova <- corr * nuisance["sill"]
  #nuisance["nugget"] <- 1-nuisance["sill"] 
  loglik <- sum(apply(stdata, 1, fname, const = const, 
                      cova = cova, dimat = dimat, ident = ident, nuisance = nuisance, 
                      setup = setup))
  #print(betac)
  #print(loglik)
  return(loglik)
}
