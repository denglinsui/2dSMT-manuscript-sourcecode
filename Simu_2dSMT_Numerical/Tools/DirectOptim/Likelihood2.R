Likelihood2 <- function (corrmodel, data,covariate, fixed, flagcor, flagnuis, grid, lower, 
          model, namescorr, namesnuis, namesparam, namescovariate, numcoord, numpairs, 
          numparamcor, numrep, numtime, optimizer, param, setup, spacetime, 
          varest, taper, type, upper){
  LogNormDenRestr <- function(const, cova, ident, dimat, nuisance, 
                              setup, stdata) {
    llik <- -1e+08
    varcov <- ident #* (nuisance["nugget"] + nuisance["sill"])
    varcov[lower.tri(varcov)] <- cova
    varcov <- t(varcov)
    varcov[lower.tri(varcov)] <- cova
    cholvarcov <- try(chol(varcov), silent = TRUE)
    if (!is.matrix(cholvarcov)) 
      return(llik)
    detvarcov <- sum(log(diag(cholvarcov)))
    ivarcov <- chol2inv(cholvarcov)
    sumvarcov <- sum(ivarcov)
    p <- ivarcov - array(rowSums(ivarcov), c(dimat, 1)) %*% 
      colSums(ivarcov)/sumvarcov
    llik <- -0.5 * (const + 2 * detvarcov + log(sumvarcov) + 
                      crossprod(t(crossprod(stdata, p)), stdata))
    return(llik)
  }
  LogNormDenTap <- function(const, cova, ident, dimat, nuisance, 
                            setup, stdata) {
    lliktap <- -1e+08
    cova[cova == (nuisance["sill"])] <- nuisance["sill"] + 
      nuisance["nugget"]
    varcovtap <- new("spam", entries = cova * setup$taps, 
                     colindices = setup$ja, rowpointers = setup$ia, dimension = as.integer(rep(dimat, 
                                                                                               2)))
    cholvctap <- try(spam::chol.spam(varcovtap), silent = TRUE)
    if (class(cholvctap) == "try-error") 
      return(lliktap)
    logdet <- c(spam::determinant(cholvctap)$modulus)
    inv <- spam::solve.spam(cholvctap)
    slot(varcovtap, "entries") <- inv[setup$idx] * setup$taps
    lliktap = -0.5 * (const + 2 * logdet + drop(t(stdata) %*% 
                                                  varcovtap %*% stdata))
    return(lliktap)
  }
  LogNormDenStand <- function(const, cova, ident, dimat, nuisance, 
                              setup, stdata) {
    llik <- -1e+08
    varcov <- ident #* (nuisance["nugget"] + nuisance["sill"])
    varcov[lower.tri(varcov)] <- cova
    varcov <- t(varcov)
    varcov[lower.tri(varcov)] <- cova
    cholvarcov <- try(chol(varcov), silent = TRUE)
    if (!is.matrix(cholvarcov)) 
      return(llik)
    detvarcov <- sum(log(diag(cholvarcov)))
    llik <- -0.5 * (const + 2 * detvarcov + sum(stdata * 
                                                  backsolve(cholvarcov, forwardsolve(cholvarcov, stdata, 
                                                                                     transpose = TRUE, upper.tri = TRUE))))
    return(llik)
  }
  
  loglik <- function(const, corr, corrmat, corrmodel, data, covariate,
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
  loglik_tap_comp <- function(const, corr, corrmat, corrmodel, 
                              data, dimat, fixed, fname, grid, ident, model, namescorr, 
                              namesnuis, param, setup) {
    loglik <- -1e+08
    param <- c(param, fixed)
    paramcorr <- param[namescorr]
    nuisance <- param[namesnuis]
    stdata <- data - nuisance["mean"]
    cc = .C(corrmat, corr = corr, as.integer(corrmodel), 
            as.double(nuisance), as.double(paramcorr), DUP = TRUE, 
            NAOK = TRUE)$corr
    corr <- cc
    if (corr[1] == -2) 
      return(loglik)
    cova <- corr * nuisance["sill"]
    loglik <- LogNormDenTap(const = const, cova = cova, dimat = dimat, 
                            ident = ident, nuisance = nuisance, setup = setup, 
                            stdata = stdata)
    sumloglik = sumloglik + loglik
    return(sumloglik)
  }
  dimat <- numcoord * numtime
  numpairstot <- dimat * (dimat - 1)/2
  const <- dimat * log(2 * pi)
  corr <- double(numpairstot)
  corrmat <- "CorrelationMat"
  ident <- diag(dimat)
  if (spacetime) {
    data = matrix(c(data), numrep, dimat, byrow = T)
    corrmat <- "CorrelationMat_st"
  }
  if (type == 3) {
    fname <- "LogNormDenRestr"
    const <- const - log(2 * pi)
  }
  if (type == 4) 
    fname <- "LogNormDenStand"
  if (type == 5) {
    corrmat <- "CorrelationMat_tap"
    if (spacetime) 
      corrmat <- "CorrelationMat_st_tap"
    fname <- "LogNormDenTap"
    corr <- double(numpairs)
    tapcorr <- double(numpairs)
    tapmod <- setup$tapmodel
    tp = .C(corrmat, tapcorr = tapcorr, as.integer(tapmod), 
            as.double(c(0, 0, 1)), as.double(1), DUP = TRUE, 
            NAOK = TRUE)$tapcorr
    setup$taps <- tp
  }
  if (optimizer == "L-BFGS-B") {
    param = param[!(names(param) %in% c("mean","nugget"))]
    param["scale"] <- max(min(param[["scale"]],100),0.01)
    #param["sill"] <- max(min(param[["sill"]],0.9),0.1)
    lower = rep(-Inf,length(param)); lower[names(param)=="sill"]=0
    upper = rep(Inf,length(param)); upper[names(param)=="sill"]=1
    
    #===== For Compare
    param_exact = param
    param_exact["mean"] = 0
    param_exact["scale"] = 0.1
    param_exact["sill"] = 0.7
    param_exact["nugget"] = 0.3
    param_exact[paste0("beta",1:12)] = 0
    if(T){
    loglik.exact <- loglik(const = const, corr = corr, corrmat = corrmat, corrmodel = corrmodel,
                            data=data, covariate = covariate,
                            dimat=dimat, fixed=fixed, fname=fname, grid=grid, ident = ident, 
                            model=model, namescorr=namescorr,  namescov=namescovariate, namesnuis = namesnuis,
                            param = param_exact, setup=setup)
    }
    #======End Compare
    print("Start Optimizationing.....")
    Likelihood <- optim(param, loglik, const = const, corr = corr, 
                        corrmat = corrmat, corrmodel = corrmodel, control = list(fnscale = -1, 
                                                                                 factr = 1, pgtol = 1e-14, maxit = 1e+08), data = data, 
                        dimat = dimat, fixed = fixed, fname = fname, grid = grid, covariate = covariate,
                        ident = ident, lower = lower, method = optimizer, 
                        model = model, namescorr = namescorr, namesnuis = namesnuis, namescov = namescovariate,
                        upper = upper, setup = setup)
    
    print("End Optimizationing.....")
  }else Likelihood <- optim(param, loglik, const = const, corr = corr, 
                           corrmat = corrmat, corrmodel = corrmodel, control = list(fnscale = -1, 
                                                                                    factr = 1, pgtol = 1e-14, maxit = 1e+08), data = data, 
                           dimat = dimat, fixed = fixed, fname = fname, grid = grid, covariate = covariate,
                           ident = ident, method = optimizer, model = model, namescorr = namescorr, namescov = namescovariate,
                           namesnuis = namesnuis, setup = setup)
  param <- Likelihood$par
  Likelihood$par["nugget"] = 1- param[["sill"]]
  numparam <- length(param)
  Likelihood$clic <- NULL
  if (type == 3 || type == 4) 
    Likelihood$clic <- -2 * (Likelihood$value - numparam)
  if (Likelihood$convergence == 0) 
    Likelihood$convergence <- "Successful"
  else if (Likelihood$convergence == 1) 
    Likelihood$convergence <- "Iteration limit reached"
  else Likelihood$convergence <- "Optimization may have failed"
  if (varest) {
    gnames <- namesparam
    if (flagnuis[1]) {
      numparam <- numparam - 1
      gnames <- gnames[gnames != "mean"]
    }
    param <- c(param, fixed)
    paramcorr <- param[namescorr]
    numparamcorr <- length(paramcorr)
    nuisance <- param[namesnuis]
    eps <- (.Machine$double.eps)^(1/3)
    numfish <- numparam * (numparam - 1)/2 + numparam
    cc = .C(corrmat, corr = corr, as.integer(corrmodel), 
            as.double(nuisance), as.double(paramcorr), DUP = TRUE, 
            NAOK = TRUE)$corr
    corr <- cc
    varian <- corr * nuisance["sill"]
    if (!spacetime) 
      dname <- "DCorrelationMat"
    else dname <- "DCorrelationMat_st"
    numparamcorr <- length(paramcorr[flagcor == 1])
    namescorr <- namescorr[flagcor == 1]
    if (type == 3 || type == 4) {
      varcov <- (nuisance["nugget"] + nuisance["sill"]) * 
        ident
      varcov[lower.tri(varcov)] <- varian
      varcov <- t(varcov)
      varcov[lower.tri(varcov)] <- varian
      cholcov <- chol(varcov)
      invar <- chol2inv(cholcov)
      fish <- double(numfish)
      if (type == 3) 
        P <- invar - array(rowSums(invar), c(dimat, 1)) %*% 
        colSums(invar)/sum(invar)
      gradient <- array(0, dim = c(dimat, numparam, dimat))
      colnames(gradient) <- gnames
      dcorr <- double(numpairstot * numparamcorr)
      dc = .C(dname, as.integer(corrmodel), dcorr = dcorr, 
              as.double(eps), as.integer(flagcor), as.integer(numparamcorr), 
              as.double(paramcorr), corr, DUP = TRUE, NAOK = TRUE)$dcorr
      dcorr <- dc
      dim(dcorr) <- c(numpairstot, numparamcorr)
      if (flagnuis[2]) 
        gradient[, namesnuis[2], ] <- ident
      if (flagnuis[3]) {
        R <- ident
        R[lower.tri(R)] <- corr
        R <- t(R)
        R[lower.tri(R)] <- corr
        gradient[, namesnuis[3], ] <- R
      }
      for (i in 1:numparamcorr) {
        grad <- matrix(0, dimat, dimat)
        grad[lower.tri(grad)] <- dcorr[, i]
        grad <- t(grad)
        grad[lower.tri(grad)] <- dcorr[, i]
        if (flagcor[namescorr][i]) 
          gradient[, namescorr[i], ] <- grad
      }
      i <- 1
      k <- 1
      for (i in 1:numparam) for (j in i:numparam) {
        if (type == 3) 
          fish[k] <- sum(diag(P %*% gradient[, i, ] %*% 
                                P %*% gradient[, j, ]))/2
        if (type == 4) 
          fish[k] <- sum(diag(invar %*% gradient[, i, 
                                                 ] %*% invar %*% gradient[, j, ]))/2
        k <- k + 1
      }
      fisher <- diag(0, numparam)
      fisher[lower.tri(fisher, diag = TRUE)] <- fish
      fisher <- t(fisher)
      fisher[lower.tri(fisher, diag = TRUE)] <- fish
      if (flagnuis[1]) {
        if (type == 4) 
          fishmean <- sum(invar)
        zeros <- rep(0, numparam)
        fisher <- rbind(c(fishmean, zeros), cbind(zeros, 
                                                  fisher))
      }
      invfisher <- try(solve(fisher), silent = TRUE)
      if (!is.matrix(invfisher)) 
        invfisher <- NULL
      Likelihood$sensmat <- NULL
      Likelihood$varimat <- NULL
    }
    if (type == 5) {
      varcov <- varian
      varcov[varcov == (nuisance["sill"])] <- nuisance["nugget"] + 
        nuisance["sill"]
      spamvar <- new("spam", entries = varcov, colindices = setup$ja, 
                     rowpointers = setup$ia, dimension = as.integer(rep(dimat, 
                                                                        2)))
      varcov <- as.matrix(spamvar)
      covtap <- new("spam", entries = varcov[setup$idx] * 
                      setup$taps, colindices = setup$ja, rowpointers = setup$ia, 
                    dimension = as.integer(rep(dimat, 2)))
      cholcovtap <- chol(covtap)
      invtap <- as.matrix(spam::solve.spam(cholcovtap))
      covtap <- as.matrix(covtap)
      HH <- double(numfish)
      JJ <- double(numfish)
      gradient <- array(0, c(numpairs, numparam))
      colnames(gradient) <- gnames
      if (!spacetime) 
        dname <- "DCorrelationMat_tap"
      else dname <- "DCorrelationMat_st_tap"
      dcorr <- double(numpairs * numparamcorr)
      dc = .C(dname, as.integer(corrmodel), dcorr = dcorr, 
              as.double(eps), as.integer(flagcor), as.integer(numparamcorr), 
              as.double(paramcorr), corr, DUP = TRUE, NAOK = TRUE)$dcorr
      dcorr <- dc
      dim(dcorr) <- c(numpairs, numparamcorr)
      if (flagnuis[2]) 
        gradient[, namesnuis[2]] <- ident[setup$idx]
      if (flagnuis[3]) 
        gradient[, namesnuis[3]] <- corr
      gradient[, namescorr] <- dcorr
      k <- 1
      H <- diag(0, numparam)
      J <- diag(0, numparam)
      for (i in 1:numparam) for (j in i:numparam) {
        gradtapI <- gradtapJ <- bigI <- bigJ <- spamvar
        slot(gradtapI, "entries") <- gradient[, i] * 
          setup$taps
        slot(gradtapJ, "entries") <- gradient[, j] * 
          setup$taps
        HH[k] <- -0.5 * sum(diag(gradtapI %*% invtap %*% 
                                   gradtapJ %*% invtap))
        slot(bigI, "entries") <- (invtap %*% gradtapI %*% 
                                    invtap)[setup$idx] * setup$taps
        slot(bigJ, "entries") <- (invtap %*% gradtapJ %*% 
                                    invtap)[setup$idx] * setup$taps
        JJ[k] <- 0.5 * sum(diag(bigI %*% varcov %*% bigJ %*% 
                                  varcov))
        k <- k + 1
      }
      H[lower.tri(H, diag = TRUE)] <- HH
      H <- t(H)
      H[lower.tri(H, diag = TRUE)] <- HH
      J[lower.tri(J, diag = TRUE)] <- JJ
      J <- t(J)
      J[lower.tri(J, diag = TRUE)] <- JJ
      if (flagnuis[1]) {
        slot(spamvar, "entries") <- invtap[setup$idx] * 
          setup$taps
        zeros <- rep(0, numparam)
        fishmH <- sum(spamvar)
        H <- rbind(c(fishmH, zeros), cbind(zeros, H))
        fishmJ <- sum(spamvar %*% varcov %*% spamvar)
        J <- rbind(c(fishmJ, zeros), cbind(zeros, J))
      }
      invH <- try(solve(H), silent = TRUE)
      Likelihood$sensmat <- H
      Likelihood$varimat <- J
      if (!is.matrix(invH) || !is.matrix(H)) {
        invfisher <- NULL
        Likelihood$clic <- NULL
      }
      else {
        invfisher <- invH %*% J %*% invH
        Likelihood$clic <- -2 * (Likelihood$value - sum(diag(J %*% 
                                                               invH)))
      }
    }
    Likelihood$varcov <- invfisher
    if (is.null(Likelihood$varcov)) {
      warning("Asymptotic information matrix is singular")
      Likelihood$varcov <- "none"
      Likelihood$stderr <- "none"
    }
    else {
      dimnames(Likelihood$varcov) <- list(namesparam, namesparam)
      Likelihood$stderr <- diag(Likelihood$varcov)
      if (any(Likelihood$stderr < 0)) 
        Likelihood$stderr <- "none"
      else {
        Likelihood$stderr <- sqrt(Likelihood$stderr)
        names(Likelihood$stderr) <- namesparam
      }
    }
  }
  Likelihood$loglik.exact = loglik.exact
  return(Likelihood)
}
