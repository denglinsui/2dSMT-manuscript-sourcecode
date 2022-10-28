#source('~/project/2DSpatial/Tools/DirectOptim/Likelihood2.R')
source('~/project/2DSpatial/Tools/DirectOptim/LikelihoodTest.R')
source('~/project/2DSpatial/Tools/DirectOptim/CheckInput2.R')
source('~/project/2DSpatial/Tools/DirectOptim/WlsInit2.R')
source('~/project/2DSpatial/Tools/DirectOptim/InitParam2.R')
source('~/project/2DSpatial/Tools/DirectOptim/loglik2.R')

FitComposite2 <- function (data, coordx, covariate, coordy = NULL, coordt = NULL, corrmodel, 
          distance = "Eucl", fixed = NULL, grid = FALSE, likelihood = "Marginal", 
          margins = "Gev", maxdist = NULL, maxtime = NULL, model = "Gaussian", 
          optimizer = "Nelder-Mead", replicates = 1, start = NULL, 
          taper = NULL, tapsep = NULL, threshold = NULL, type = "Pairwise", 
          varest = FALSE, vartype = "SubSamp", weighted = FALSE, winconst, 
          winstp) 
{
  call <- match.call()
  checkinput <- CheckInput2(coordx, coordy, coordt, covariate,corrmodel, data, distance, 
                           "Fitting", fixed, grid, likelihood, margins, maxdist, maxtime, 
                           model, NULL, optimizer, NULL, replicates, start, taper, 
                           tapsep, threshold, type, varest, vartype, 
                           weighted)
  
  if (!is.null(checkinput$error)) 
    stop(checkinput$error)
  FitComposite <- NULL
  score <- sensmat <- varcov <- varimat <- parscale <- NULL
  initparam <- WlsInit2(coordx, coordy, coordt,covariate, corrmodel, data, distance, 
                       "Fitting", fixed, grid, likelihood, margins, 
                       maxdist, maxtime, model, NULL, NULL, parscale, optimizer == 
                         "L-BFGS-B", replicates, start, taper, tapsep, threshold, 
                       type, varest, vartype, weighted, winconst, winstp
                       )
  if (!is.null(initparam$error)) 
    stop(initparam$error)
  if (likelihood == "Full") {
    #initparam$param <- initparam$param[-c("nugget")]
    fitted <- Likelihood2(initparam$corrmodel, initparam$data, initparam$covariate,
                         initparam$fixed, 
                         initparam$flagcorr, initparam$flagnuis, 
                         grid, initparam$lower, initparam$model, initparam$namescorr, 
                         initparam$namesnuis, initparam$namesparam, initparam$namescovariate, initparam$numcoord, 
                         initparam$numpairs, initparam$numparamcorr, initparam$numrep, 
                         initparam$numtime, optimizer, initparam$param, initparam$setup, 
                         initparam$spacetime, varest, taper, initparam$type, 
                         initparam$upper)
  }
  if (likelihood == "Marginal" || likelihood == "Conditional") {
    fitted <- CompLikelihood(initparam$coordx, initparam$coordy, 
                             initparam$corrmodel, initparam$data, initparam$distance, 
                             initparam$flagcorr, initparam$flagnuis, initparam$fixed, 
                             grid, initparam$likelihood, initparam$lower, initparam$model, 
                             initparam$namescorr, initparam$namesnuis, initparam$namesparam, 
                             initparam$numparam, initparam$numparamcorr, optimizer, 
                             initparam$param, initparam$spacetime, initparam$threshold, 
                             initparam$type, initparam$upper, varest, initparam$vartype, 
                             initparam$winconst, initparam$winstp)
  }
  .C("DeleteGlobalVar", DUP = TRUE, NAOK = TRUE)
  FitComposite <- list(clic = fitted$clic, coordx = initparam$coordx, 
                       coordy = initparam$coordy, coordt = initparam$coordt, 
                       convergence = fitted$convergence, corrmodel = corrmodel, 
                       data = initparam$data, distance = distance, fixed = initparam$fixed, 
                       grid = grid, iterations = fitted$counts, likelihood = likelihood, 
                       logCompLik = fitted$value, message = fitted$message, 
                       model = model, numcoord = initparam$numcoord, numrep = initparam$numrep, 
                       numtime = initparam$numtime, param = fitted$par, nozero = initparam$setup$nozero, 
                       score = fitted$score, srange = initparam$srange, stderr = fitted$stderr, 
                       sensmat = fitted$sensmat, varcov = fitted$varcov, varimat = fitted$varimat, 
                       vartype = vartype, trange = initparam$trange, threshold = initparam$threshold, 
                       type = type, winconst = fitted$winconst, winstp = fitted$winstp,
                       loglik.exact = fitted$loglik.exact,
                       init.param = initparam$param)
  structure(c(FitComposite, call = call), class = c("FitComposite"))
}
