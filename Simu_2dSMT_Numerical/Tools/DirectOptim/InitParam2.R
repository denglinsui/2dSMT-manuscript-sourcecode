InitParam2 <- function (coordx, coordy, coordt,covariate, corrmodel, data, distance, 
          fcall, fixed, grid, likelihood, margins, maxdist, maxtime, 
          model, numblock, param, parscale, paramrange, replicates, 
          start, taper, tapsep, threshold, type, typereal, varest, 
          vartype, weighted, winconst, winstp
          ){
  CheckSpaceTime <- function(corrmodel) {
    CheckSpaceTime <- NULL
    if (corrmodel >= 1 & corrmodel <= 20) 
      CheckSpaceTime <- FALSE
    else CheckSpaceTime <- TRUE
    return(CheckSpaceTime)
  }
  CheckDistance <- function(distance) {
    CheckDistance <- NULL
    CheckDistance <- switch(distance, eucl = 0, Eucl = 0, 
                            chor = 1, Chor = 1, geod = 2, Geod = 2, proj = 3, 
                            Proj = 3)
    return(CheckDistance)
  }
  SetRangeParam <- function(namesparam, numparam) {
    low <- 1e-12
    lower <- NULL
    upper <- NULL
    for (i in 1:numparam) {
      if (namesparam[i] == "mean") {
        lower <- c(lower, -Inf)
        upper <- c(upper, Inf)
      }
      if (namesparam[i] == "nugget") {
        lower <- c(lower, 0)
        upper <- c(upper, Inf)
      }
      if (namesparam[i] == "power") {
        lower <- c(lower, low)
        upper <- c(upper, 2)
      }
      if (namesparam[i] == "power_s") {
        lower <- c(lower, low)
        upper <- c(upper, 2)
      }
      if (namesparam[i] == "power_t") {
        lower <- c(lower, low)
        upper <- c(upper, 2)
      }
      if (namesparam[i] == "power1") {
        lower <- c(lower, low)
        upper <- c(upper, 2)
      }
      if (namesparam[i] == "power2") {
        lower <- c(lower, low)
        upper <- c(upper, Inf)
      }
      if (namesparam[i] == "scale") {
        lower <- c(lower, low)
        upper <- c(upper, Inf)
      }
      if (namesparam[i] == "scale_s") {
        lower <- c(lower, low)
        upper <- c(upper, Inf)
      }
      if (namesparam[i] == "scale_t") {
        lower <- c(lower, low)
        upper <- c(upper, Inf)
      }
      if (namesparam[i] == "sep") {
        lower <- c(lower, low)
        upper <- c(upper, 1)
      }
      if (namesparam[i] == "sill") {
        lower <- c(lower, low)
        upper <- c(upper, Inf)
      }
      if (namesparam[i] == "smooth") {
        lower <- c(lower, low)
        upper <- c(upper, Inf)
      }
    }
    return(list(lower = lower, upper = upper))
  }
  error <- NULL
  namesfixed <- namesstart <- namessim <- NULL
  numfixed <- numstart <- 0
  namesnuis <- NuisanceParam(model)
  model <- CheckModel(model)
  flagnuis <- NULL
  namescorr <- CorrelationParam(corrmodel)
  numparamcorr <- length(namescorr)
  paramcorr <- rep(1, numparamcorr)
  names(paramcorr) <- namescorr
  flagcorr <- NULL
  corrmodel <- CheckCorrModel(corrmodel)
  spacetime <- CheckSpaceTime(corrmodel)
  if (is.null(coordy)) {
    coordy <- coordx[, 2]
    coordx <- coordx[, 1]
  }
  if (grid) {
    numcoordx <- length(coordx)
    numcoordy <- length(coordy)
    numcoord <- numcoordx * numcoordy
  }else numcoord <- numcoordx <- numcoordy <- length(coordx)
  tapering <- ia <- idx <- ja <- integer(1)
  nozero <- NULL
  tapmodel = NULL
  cutoff <- FALSE
  distance <- CheckDistance(distance)
  if (fcall == "Fitting") {
    nuisance <- NULL
    likelihood <- CheckLikelihood(likelihood)
    vartype <- CheckVarType(vartype)
    type <- CheckType(type)
    if (model == 1) {
      mu <- mean(data)
      if (any(type == c(1, 3, 6))) 
        if (is.list(fixed)) 
          fixed$mean <- mu
        else fixed <- list(mean = mu)
        nuisance <- c(mu, 0, var(c(data)))
        if (likelihood == 2 && CheckType(typereal) == 5) 
          tapering <- 1
    }
    if (model == 2) {
      if (is.null(threshold) || !is.numeric(threshold)) 
        threshold <- 0
      p <- mean(data)
      mu <- threshold + qnorm(p)
      nuisance <- c(mu, 0, 1)
      if (!is.null(start$nugget)) 
        if (length(start) > 1) 
          start <- start[!names(start) %in% "nugget"]
      else start <- NULL
      if (is.list(fixed)) 
        fixed$nugget <- 0
      else fixed <- list(nugget = 0)
      if (!is.null(fixed$sill)) 
        fixed$nugget <- 1 - fixed$sill
    }
    if (model > 2) {
      if (model == 3) {
        if (is.list(fixed)) 
          fixed$sill <- 1
        else fixed <- list(sill = 1)
        if (corrmodel == 4) 
          fixed$power1 <- 2
      }
      if (model == 5) 
        nuisance <- c(nuisance, 1)
      nuisance <- c(nuisance, 0.5)
      if (margins == "Gev") 
        data <- Dist2Dist(data)
    }
    names(nuisance) <- namesnuis
    namescovariate <- paste0("beta",1:dim(covariate)[2])
    paramcovariate <- rep(1,dim(covariate)[2])
    names(paramcovariate) <- namescovariate
    namesparam <- sort(c(namescorr, namesnuis, namescovariate))
    param <- c(nuisance, paramcorr, paramcovariate)
    param <- param[namesparam]
    numparam <- length(param)
    flag <- rep(1, numparam)
    namesflag <- namesparam
    names(flag) <- namesflag
    if (!is.null(fixed)) {
      fixed <- unlist(fixed)
      namesfixed <- names(fixed)
      numfixed <- length(namesfixed)
      if (numfixed == numparam) {
        error <- "the are not parameters left to estimate\n"
        return(list(error = error))
      }
      flag[pmatch(namesfixed, namesflag)] <- 0
      param <- param[-pmatch(namesfixed, namesparam)]
      numparamcorr <- numparamcorr - sum(namesfixed %in% 
                                           namescorr)
      namesparam <- names(param)
      numparam <- length(param)
    }
    flagcorr <- flag[namescorr]
    flagnuis <- flag[namesnuis]
    if (!is.null(start)) {
      start <- unlist(start)
      namesstart <- names(start)
      if (any(type == c(1, 3, 6))) 
        if (any(namesstart == "mean")) 
          start <- start[!namesstart == "mean"]
      namesstart <- names(start)
      numstart <- length(start)
      param[pmatch(namesstart, namesparam)] <- start
    }
    if (paramrange) 
      paramrange <- SetRangeParam(namesparam, numparam)
    else paramrange <- list(lower = NULL, upper = NULL)
    if (missing(winconst) || !is.numeric(winconst)) 
      winconst <- 0
    if (missing(winstp) || !is.numeric(winstp)) 
      winstp <- 0
    if (spacetime) {
      numtime <- length(coordt)
      if (grid) 
        if (replicates > 1) {
          dim(data) <- c(numcoord, numtime, replicates)
          data <- aperm(data, c(2, 1, 3))
        }
      else data <- matrix(data, ncol = numcoord, nrow = numtime, 
                          byrow = TRUE)
      if (typereal == "Tapering") {
        tapering <- 1
        idx <- integer((numcoord * numtime)^2)
        ja <- integer((numcoord * numtime)^2)
        ia <- integer(numcoord * numtime + 1)
        tapmodel <- CheckCorrModel(taper)
      }
    }
    else {
      numtime <- 1
      coordt <- 0
      if (grid) 
        data <- matrix(data, ncol = numcoord, nrow = replicates, 
                       byrow = TRUE)
      else data <- matrix(data, ncol = numcoord, nrow = replicates)
      if (typereal == "Tapering") {
        tapering <- 1
        idx <- integer((numcoord * numtime)^2)
        ja <- integer((numcoord * numtime)^2)
        ia <- integer(numcoord * numtime + 1)
        tapmodel <- CheckCorrModel(taper)
      }
    }
  }
  if (fcall == "Simulation") {
    namesnuis <- sort(unique(c(namesnuis, NuisanceParam("Gaussian"))))
    param <- unlist(param)
    numparam <- length(param)
    namesparam <- names(param)
    if (model %in% c(3, 5)) 
      if (is.null(numblock)) 
        numblock <- as.integer(500)
    else numblock <- as.integer(numblock)
    if (model == 3) {
      param["df"] <- 0
      if (corrmodel == 2) 
        param["scale"] <- 2 * log(numblock) * param["scale"]
      if (corrmodel %in% c(3, 4)) 
        param["scale"] <- 2 * sqrt(log(numblock)) * param["scale"]
      if (corrmodel == 6) 
        param["scale"] <- 2 * (log(numblock))^(1/param["power"]) * 
          param["scale"]
    }
    namessim <- c("mean", "sill", "nugget", "scale", namescorr[!namescorr == 
                                                                 "scale"])
    if (spacetime) 
      numtime <- length(coordt)
    else {
      numtime <- 1
      coordt <- 0
    }
    if (typereal == "Tapering" && type == "Tapering") {
      tapering <- 1
      idx <- integer((numcoord * numtime)^2)
      ja <- integer((numcoord * numtime)^2)
      ia <- integer(numcoord * numtime + 1)
      tapmodel <- CheckCorrModel(taper)
    }
    if (model == 2) {
      if (is.null(threshold) || !is.numeric(threshold)) 
        threshold <- 0
    }
  }
  numpairs <- integer(1)
  srange <- double(1)
  trange <- double(1)
  if (is.null(maxdist)) 
    srange <- c(srange, double(1))
  else {
    srange <- c(srange, as.double(maxdist))
  }
  if (is.null(maxtime)) 
    trange <- c(trange, double(1))
  else {
    trange <- c(trange, as.double(maxtime))
  }
  isinit <- as.integer(1)
  SG = .C("SetGlobalVar", as.double(coordx), as.double(coordy), 
          as.double(coordt), as.integer(grid), ia = ia, idx = idx, 
          isinit = isinit, ja = ja, as.integer(numcoord), as.integer(numcoordx), 
          as.integer(numcoordy), numpairs = numpairs, as.integer(replicates), 
          srange = srange, as.double(tapsep), as.integer(numtime), 
          trange = trange, as.integer(tapering), as.integer(tapmodel), 
          as.integer(distance), as.integer(weighted), DUP = TRUE, 
          NAOK = TRUE)
  srange <- SG$srange
  trange <- SG$trange
  ja <- SG$ja
  ia <- SG$ia
  idx <- SG$idx
  isinit <- SG$isinit
  numpairs <- SG$numpairs
  if (tapering) {
    nozero <- numpairs/(numcoord * numtime)^2
    idx <- idx[1:numpairs]
    ja <- ja[1:numpairs]
  }
  return(list(coordx = coordx, coordy = coordy, coordt = coordt, 
              corrmodel = corrmodel, data = data, covariate = covariate, distance = distance, 
              error = error, flagcorr = flagcorr, flagnuis = flagnuis, 
              fixed = fixed, likelihood = likelihood, lower = paramrange$lower, 
              model = model, namescorr = namescorr, namesfixed = namesfixed, 
              namesnuis = namesnuis, namesparam = namesparam, namessim = namessim, namescovariate = namescovariate,
              namesstart = namesstart, numblock = numblock, numcoord = numcoord, 
              numcoordx = numcoordx, numcoordy = numcoordy, numfixed = numfixed, 
              numpairs = numpairs, numparam = numparam, numparamcorr = numparamcorr, 
              numrep = replicates, numstart = numstart, numtime = numtime, 
              param = param, setup = list(ia = ia, idx = idx, ja = ja, 
                                          nozero = nozero, tapmodel = tapmodel, tapsep = tapsep), 
              spacetime = spacetime, srange = srange, start = start, 
              upper = paramrange$upper, type = type, threshold = threshold, 
              trange = trange, vartype = vartype, winconst = winconst, 
              winstp = winstp))
}
