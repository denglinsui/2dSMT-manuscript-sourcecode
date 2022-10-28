WlsInit2 <- function (coordx, coordy, coordt,covariate, corrmodel, data, distance, 
          fcall, fixed, grid, likelihood, margins, maxdist, maxtime, 
          model, numblock, param, parscale, paramrange, replicates, 
          start, taper, tapsep, threshold, type, varest, vartype, weighted, 
          winconst, winstp){
  initparam <- InitParam2(coordx, coordy, coordt, covariate,corrmodel, 
                         data, distance, fcall, fixed, grid, likelihood, margins, 
                         maxdist, maxtime, model, numblock, param, parscale, paramrange, 
                         replicates, start, taper, tapsep, threshold, "WLeastSquare", 
                         type, varest, vartype, weighted,
                         winconst, winstp
                         )
  if (!is.null(initparam$error)) 
    stop(initparam$error)
  initparam$type <- CheckType(type)
  if (substr(model, 1, 6) == "Binary") 
    return(initparam)
  if (is.null(start)) 
    start <- NA
  else start <- unlist(start)
  if (is.null(fixed)) 
    fixed <- NA
  else fixed <- unlist(fixed)
  if (initparam$numstart == initparam$numparam) {
    if (model == "Gaussian" & (type %in% c("Standard", "Pairwise", 
                                           "Tapering"))) {
      if (is.na(fixed["mean"])) {
        if (is.na(start["mean"])) 
          initparam$param <- c(initparam$fixed["mean"], 
                               initparam$param)
        else initparam$param <- c(start["mean"], initparam$param)
        initparam$namesparam <- sort(names(initparam$param))
        initparam$param <- initparam$param[initparam$namesparam]
        initparam$numparam <- initparam$numparam + 1
        initparam$flagnuis["mean"] <- 1
        initparam$numfixed <- initparam$numfixed - 1
        if (initparam$numfixed > 0) 
          initparam$fixed <- fixed
        else initparam$fixed <- NULL
      }
      else initparam$fixed["mean"] <- fixed["mean"]
    }
    return(initparam)
  }
  if (initparam$numstart > 0) {
    initparam$param <- initparam$param[seq(1, initparam$numparam)[-pmatch(initparam$namesstart, 
                                                                          initparam$namesparam)]]
    initparam$fixed <- c(initparam$fixed, initparam$start)
  }
  Lsquare <- function(bins, bint, corrmodel, fixed, fun, lenbins, 
                      moments, namescorr, namesnuis, numbins, numbint, param) {
    param <- c(param, fixed)
    paramcorr <- param[namescorr]
    nuisance <- param[namesnuis]
    result <- .C(fun, as.double(bins), as.double(bint), as.integer(corrmodel), 
                 as.double(lenbins), as.double(moments), as.integer(numbins), 
                 as.integer(numbint), as.double(nuisance), as.double(paramcorr), 
                 res = double(1), DUP = TRUE, NAOK = TRUE)$res
    return(result)
  }
  fname <- NULL
  numbins <- as.integer(13)
  bins <- double(numbins)
  numvario <- numbins - 1
  moments <- double(numvario)
  lenbins <- integer(numvario)
  if (model == "ExtGauss" || model == "BrowResn" || model == 
      "ExtT") {
    fname <- "Binned_Madogram"
    data <- Dist2Dist(data, to = "Uniform")
  }
  if (initparam$spacetime) {
    numbint <- initparam$numtime - 1
    bint <- double(numbint)
    momentt <- double(numbint)
    lenbint <- integer(numbint)
    numbinst <- numvario * numbint
    binst <- double(numbinst)
    momentst <- double(numbinst)
    lenbinst <- integer(numbinst)
    if (type == "Tapering") 
      fname <- "Binned_Variogram_st2"
    else fname <- "Binned_Variogram_st"
    EV = .C(fname, bins = bins, bint = bint, as.double(initparam$data), 
            lenbins = lenbins, lenbinst = lenbinst, lenbint = lenbint, 
            moments = moments, momentst = momentst, momentt = momentt, 
            as.integer(numbins), as.integer(numbint), DUP = TRUE, 
            NAOK = TRUE)
    bins = EV$bins
    bint = EV$bint
    lenbins <- EV$lenbins
    lenbint <- EV$lenbint
    lenbinst <- EV$lenbinst
    moments <- EV$moments
    momentst <- EV$momentst
    momentt <- EV$momentt
    indbin <- lenbins > 0
    bins <- bins[indbin]
    moments <- moments[indbin]
    lenbins <- lenbins[indbin]
    numbins <- sum(indbin)
    indbint <- lenbint > 0
    bint <- bint[indbint]
    momentt <- momentt[indbint]
    lenbint <- lenbint[indbint]
    numbint <- sum(indbint)
    indbinst <- lenbinst > 0
    momentst <- momentst[indbinst]
    lenbinst <- lenbinst[indbinst]
    numbinst <- sum(indbinst)
    moment <- matrix(momentst, nrow = numbins, ncol = numbint, 
                     byrow = TRUE)
    lenbin <- matrix(lenbinst, nrow = numbins, ncol = numbint, 
                     byrow = TRUE)
    moment <- rbind(momentt, moment)
    moment <- cbind(c(0, moments), moment)
    lenbin <- rbind(lenbint, lenbin)
    lenbin <- cbind(c(1, lenbins), lenbin)
    moments <- moment
    lenbins <- lenbin
    bins <- c(-bins[1], bins)
    bint <- c(0, bint)
    numbins <- numbins + 1
    numbint <- numbint + 1
  }
  else {
    if (type == "Tapering") 
      fname <- "Binned_Variogram2"
    else fname <- "Binned_Variogram"
    numbint <- 1
    bint <- double(numbint)
    momentt <- double(1)
    momentst <- double(1)
    lenbint <- integer(1)
    lenbinst <- integer(1)
    EV = .C(fname, bins = bins, as.double(data), lenbins = lenbins, 
            moments = moments, as.integer(numbins), DUP = TRUE, 
            NAOK = TRUE)
    bins <- EV$bins
    lenbins <- EV$lenbins
    moments <- EV$moments
    centers <- bins[1:numvario] + diff(bins)/2
    indbin <- lenbins > 0
    bins <- bins[indbin]
    centers <- centers[indbin]
    moments <- moments[indbin]
    lenbins <- lenbins[indbin]
    numbins <- sum(indbin)
    variogram <- moments/lenbins
    if (!is.na(initparam$param["scale"])) 
      initparam$param["scale"] <- centers[max(variogram) == 
                                            variogram]
  }
  if (model == "Gaussian") 
    fname <- "LeastSquare_G"
  if (model == "ExtGauss") 
    fname <- "LeastSquare_MEG"
  if (model == "BrowResn") 
    fname <- "LeastSquare_MBR"
  if (model == "ExtT") 
    fname <- "LeastSquare_MET"
  fitted <- optim(initparam$param, Lsquare, bins = bins, bint = bint, 
                  corrmodel = initparam$corrmodel, fixed = initparam$fixed, 
                  fun = fname, lenbins = lenbins, moments = moments, namescorr = initparam$namescorr, 
                  namesnuis = initparam$namesnuis, numbins = numbins, numbint = numbint, 
                  control = list(fnscale = -1, reltol = 1e-14, maxit = 1e+08), 
                  hessian = FALSE)
  initparam$param <- fitted$par
  
  #initparam$param[paste0("beta",1:dim(covariate)[2])] <- rep(1,dim(covariate)[2])
  if (initparam$numstart > 0) {
    initparam$param <- c(initparam$param, initparam$start)
    initparam$fixed <- initparam$fixed[seq(1, initparam$numstart + 
                                             initparam$numfixed)[-pmatch(initparam$namesstart, 
                                                                         names(initparam$fixed))]]
    initparam$namesparam <- sort(names(initparam$param))
    initparam$param <- initparam$param[initparam$namesparam]
  }
  if (model == "Gaussian" & (type %in% c("Standard", "Pairwise", 
                                         "Tapering"))) {
    if (is.na(fixed["mean"])) {
      if (is.na(start["mean"])) 
        initparam$param <- c(initparam$fixed["mean"], 
                             initparam$param)
      else initparam$param <- c(start["mean"], initparam$param)
      initparam$namesparam <- sort(names(initparam$param))
      initparam$param <- initparam$param[initparam$namesparam]
      initparam$numparam <- initparam$numparam + 1
      initparam$flagnuis["mean"] <- 1
      initparam$numfixed <- initparam$numfixed - 1
      if (initparam$numfixed > 0) 
        initparam$fixed <- fixed
      else initparam$fixed <- NULL
    }
    else initparam$fixed["mean"] <- fixed["mean"]
  }
  return(initparam)
}
