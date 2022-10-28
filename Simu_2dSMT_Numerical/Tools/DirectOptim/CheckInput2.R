CheckInput2 <- function (coordx, coordy, coordt, covariate, corrmodel, data, distance,  
                         fcall, fixed, grid, likelihood, margins, maxdist, maxtime, 
                         model, numblock, optimizer, param, replicates, start, taper, 
                         tapsep, threshold, type, varest, vartype, weighted) {
  error <- NULL
  CheckParamRange <- function(param) {
    if (!is.na(param["df"])) 
      if (param["df"] <= 0) 
        return(FALSE)
    if (!is.na(param["nugget"])) 
      if (param["nugget"] < 0 & param["nugget"] > 1 ) 
        return(FALSE)
    if (!is.na(param["power"])) 
      if (param["power"] <= 0 || param["power"] > 2) 
        return(FALSE)
    if (!is.na(param["power_s"])) 
      if (param["power_s"] <= 0 || param["power_s"] > 2) 
        return(FALSE)
    if (!is.na(param["power_t"])) 
      if (param["power_t"] <= 0 || param["power_t"] > 2) 
        return(FALSE)
    if (!is.na(param["power1"])) 
      if (param["power1"] <= 0 || param["power1"] > 2) 
        return(FALSE)
    if (!is.na(param["power2"])) 
      if (param["power2"] <= 0) 
        return(FALSE)
    if (!is.na(param["sep"])) 
      if (param["sep"] < 0 || param["sep"] > 1) 
        return(FALSE)
    if (!is.na(param["scale"])) 
      if (param["scale"] <= 0) 
        return(FALSE)
    if (!is.na(param["scale_s"])) 
      if (param["scale_s"] <= 0) 
        return(FALSE)
    if (!is.na(param["scale_t"])) 
      if (param["scale_t"] <= 0) 
        return(FALSE)
    if (!is.na(param["sill"])) 
      if (param["sill"] <= 0 & param["sill"] > 1) 
        return(FALSE)
    if (!is.na(param["smooth"])) 
      if (param["smooth"] <= 0) 
        return(FALSE)
    if (!is.na(param["smooth_s"])) 
      if (param["smooth_s"] <= 0) 
        return(FALSE)
    if (!is.na(param["smooth_t"])) 
      if (param["smooth_t"] <= 0) 
        return(FALSE)
    return(TRUE)
  }
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
  if (fcall != "Kriging") {
    if (missing(coordx) || !is.numeric(coordx)) {
      error <- "insert a suitable set of numeric coordinates\n"
      return(list(error = error))
    }
    if (!is.null(coordy) & !is.numeric(coordy)) {
      error <- "insert a suitable set of numeric coordinates\n"
      return(list(error = error))
    }
    if (missing(corrmodel) || !is.character(corrmodel)) {
      error <- "insert the correlation model\n"
      return(list(error = error))
    }
    if (!is.null(grid) & !is.logical(grid)) {
      error <- "the parameter grid need to be a logic value\n"
      return(list(error = error))
    }
    if (!is.null(model) & !is.character(model)) {
      error <- "insert the name of the random field\n"
      return(list(error = error))
    }
    if (is.null(CheckModel(model))) {
      error <- "the model name of the random field is not correct\n"
      return(list(error = error))
    }
    if (is.null(replicates) || (abs(replicates - round(replicates)) > 
                                0) || replicates < 1) {
      error <- "the parameter replicates need to be a positive integer\n"
      return(list(error = error))
    }
    if (is.null(CheckDistance(distance))) {
      error <- "the name of the distance is not correct\n"
      return(list(error = error))
    }
    if (is.null(CheckCorrModel(corrmodel))) {
      error <- "the name of the correlation model is not correct\n"
      return(list(error = error))
    }
    if (fcall == "Fitting") {
      if (missing(data) || !is.numeric(data)) {
        error <- "insert a numeric vector or matrix of data\n"
        return(list(error = error))
      }
      if (!is.null(fixed) & !is.list(fixed)) {
        error <- "insert fixed values as a list of parameters\n"
        return(list(error = error))
      }
      if (!is.null(fixed)) {
        namfixed <- names(fixed)
        if (!all(namfixed %in% c(NuisanceParam(model), 
                                 CorrelationParam(corrmodel),
                                 paste0("beta",1:dim(covariate)[2])))) {
          error <- "some names of the fixed parameters is/are not correct\n"
          return(list(error = error))
        }
        if (!CheckParamRange(unlist(fixed))) {
          error <- "some fixed values are out of the range\n"
          return(list(error = error))
        }
      }
      if (!is.null(likelihood) & !is.character(likelihood)) {
        error <- "insert the type of likelihood objects\n"
        return(list(error = error))
      }
      if (!is.null(maxdist)) {
        error <- "insert a positive numeric value for the maximum spatial distance\n"
        if (!is.numeric(maxdist)) 
          return(list(error = error))
        else if (maxdist < 0) 
          return(list(error = error))
      }
      if (!is.null(maxtime)) {
        error <- "insert a positive numeric value for the maximum time interval\n"
        if (!is.numeric(maxtime)) 
          return(list(error = error))
        else if (maxtime < 0) 
          return(list(error = error))
      }
      if (model == "BinaryGauss") {
        if (length(unique(c(data))) != 2) {
          error <- "the data are not binary"
          return(list(error = error))
        }
        if (!is.null(start$sill)) 
          if (start$sill > 1) {
            error <- "some starting values are out of the range\n"
            return(list(error = error))
          }
        if (!is.null(fixed$sill)) 
          if (fixed$sill > 1) {
            error <- "some starting values are out of the range\n"
            return(list(error = error))
          }
      }
      if (model %in% c("BrownResn", "ExtGauss", "ExtT")) 
        if (!margins %in% c("Frechet", "Gev")) {
          error <- "insert the correct type of marginal distributions\n"
          return(list(error = error))
        }
      if (model == "BrowResn") 
        if (CheckSpaceTime(CheckCorrModel(corrmodel))) {
          if (!corrmodel %in% c("exp_exp", "exp_gauss", 
                                "stable_stable")) {
            error <- "the correlation model is not adequate for the Brown-Resnick process\n"
            return(list(error = error))
          }
        }
      else if (!corrmodel %in% c("exponential", "gauss", 
                                 "gencauchy", "stable")) {
        error <- "the correlation model is not adequate for the Brown-Resnick process\n"
        return(list(error = error))
      }
      if (!is.null(optimizer) & !is.character(optimizer)) {
        error <- "insert the type of maximising algorithm\n"
        return(list(error = error))
      }
      if (!is.null(varest) & !is.logical(varest)) {
        error <- "the parameter std.err need to be a logical value\n"
        return(list(error = error))
      }
      if (type == "Tapering") {
        if (is.null(taper) || is.null(maxdist)) {
          error <- "tapering need a taper correlation model and/or a compact support\n"
          return(list(error = error))
        }
        if (!taper %in% c("Wendland1", "Wendland2", "Wendland3", 
                          "Wendland1_Wendland1", "Wendland1_Wendland2", 
                          "Wendland1_Wendland3", "Wendland2_Wendland1", 
                          "Wendland2_Wendland2", "Wendland2_Wendland3", 
                          "Wendland3_Wendland1", "Wendland3_Wendland2", 
                          "Wendland3_Wendland3", "qt_time", "qt_space")) {
          error <- "insert a correct name for the taper correlation model\n"
          return(list(error = error))
        }
      }
      if (!is.null(type) & !is.character(type)) {
        error <- "insert the configuration of the likelihood objects\n"
        return(list(error = error))
      }
      if (is.null(CheckType(type))) {
        error <- "the type name of the likelihood objects is not correct\n"
        return(list(error = error))
      }
      if (is.null(CheckLikelihood(likelihood))) {
        error <- "the setting name of the likelihood objects is not correct\n"
        return(list(error = error))
      }
      if (likelihood == "Full") {
        if (!any(type == c("Restricted", "Standard", 
                           "Tapering"))) {
          error <- "insert a type name of the likelihood objects compatible with the full likelihood\n"
          return(list(error = error))
        }
      }
      if (likelihood == "Marginal") {
        if (!any(type == c("Difference", "Pairwise"))) {
          error <- "insert a type name of the likelihood objects compatible with the composite-likelihood\n"
          return(list(error = error))
        }
      }
      if (varest & (likelihood == "Conditional" || likelihood == 
                    "Marginal") & (!is.null(vartype) & !is.character(vartype))) {
        error <- "insert the type of estimation method for the variances\n"
        return(list(error = error))
      }
      if (varest & is.null(CheckVarType(vartype)) & (likelihood == 
                                                     "Conditional" || likelihood == "Marginal")) {
        error <- "the name of the estimation type for the variances is not correct\n"
        return(list(error = error))
      }
      if (!is.null(tapsep)) {
        error <- "insert a numeric value between 0 and 1 for the separability parameter of the space time taper\n"
        if (!is.numeric(tapsep)) 
          return(list(error = error))
        else if (tapsep < 0 || tapsep > 1) 
          return(list(error = error))
      }
      if (!is.null(start)) {
        if (!is.list(start)) {
          error <- "insert starting values as a list of parameters\n"
          return(list(error = error))
        }
        namstart <- names(start)
        if (!all(namstart %in% c(NuisanceParam(model), 
                                 CorrelationParam(corrmodel),
                                 paste0("beta",1:dim(covariate)[2])))) {
          error <- "some names of the starting parameters is/are not correct\n"
          return(list(error = error))
        }
        if (any(namstart == "mean") & (type == "Difference" || 
                                       type == "Restricted")) {
          error <- "the mean parameter is not allow with the difference composite likelihood\n"
          return(list(error = error))
        }
        if (corrmodel == "gencauchy" && param["power1"] != 
            2) {
          error <- "the parameter power1 need to be equal to 2 for the Brown-Resnick process\n"
          return(list(error = error))
        }
        if (any(namstart == "sill") & (model == "BrowResn")) {
          error <- "the sill parameter is not allow with Brown-Renick model\n"
          return(list(error = error))
        }
        if (!CheckParamRange(unlist(start))) {
          error <- "some starting values are out of the range\n"
          return(list(error = error))
        }
        if (!is.null(fixed)) 
          if (any(namstart %in% namfixed)) {
            error <- "some fixed parameter name/s is/are matching with starting parameter name/s\n"
            return(list(error = error))
          }
      }
      if (!is.null(weighted) & !is.logical(weighted)) {
        error <- "insert if the composite likelihood need to be weighted"
        return(list(error = error))
      }
      dimdata <- dim(data)
      if (is.null(coordt)) {
        if (CheckSpaceTime(CheckCorrModel(corrmodel))) {
          error <- "temporal coordinates are missing\n"
          return(list(error = error))
        }
        if (replicates > 1) {
          if (grid) {
            if (is.null(dimdata)) {
              error <- c("insert an array d x d of n iid spatial observations\n")
              return(list(error = error))
            }
            if (length(dimdata) != 3) {
              error <- c("the dimension of the data matrix is not correct\n")
              return(list(error = error))
            }
            if (length(coordx) != dimdata[1] || length(coordy) != 
                dimdata[2]) {
              error <- c("the number of coordinates does not match with the number of spatial observations\n")
              return(list(error = error))
            }
            if (dimdata[3] != replicates) {
              error <- c("the number of replications does not match with the data observations\n")
              return(list(error = error))
            }
          }
          else {
            if (is.null(dimdata)) {
              error <- c("insert a matrix n x d of spatial observations\n")
              return(list(error = error))
            }
            if (length(dimdata) != 2) {
              error <- c("the dimension of the data matrix is not correct\n")
              return(list(error = error))
            }
            if (is.null(coordy)) {
              if (is.null(dim(coordx))) {
                error <- c("insert a matrix d x 2 of spatial coordinates\n")
                return(list(error = error))
              }
              if (dimdata[2] != nrow(coordx) || ncol(coordx) != 
                  2) {
                error <- c("the number of coordinates does not match with the number of spatial observations\n")
                return(list(error = error))
              }
            }
            else if (length(coordx) != length(coordy)) {
              error <- c("the number of the two coordinates does not match\n")
              return(list(error = error))
            }
            if (dimdata[1] != replicates) {
              error <- c("the number of replications does not match with the data observations\n")
              return(list(error = error))
            }
          }
        }
        else {
          if (grid) {
            if (is.null(dimdata)) {
              error <- c("insert a matrix d x d of spatial observations\n")
              return(list(error = error))
            }
            if (length(dimdata) != 2) {
              error <- c("the dimension of the data matrix is not correct\n")
              return(list(error = error))
            }
            if (length(coordx) != dimdata[1] || length(coordy) != 
                dimdata[2]) {
              error <- c("the number of coordinates does not match with the number of spatial observations\n")
              return(list(error = error))
            }
          }
          else {
            numsite <- length(data)
            if (is.null(numsite)) {
              error <- c("insert a vector of spatial observations\n")
              return(list(error = error))
            }
            if (is.null(coordy)) {
              dimcoord <- dim(coordx)
              if (is.null(dimcoord)) {
                error <- c("insert a suitable set of coordinates\n")
                return(list(error = error))
              }
              else {
                if (dimcoord[1] != numsite || dimcoord[2] != 
                    2) {
                  error <- c("the number of coordinates does not match with the number of spatial observations\n")
                  return(list(error = error))
                }
              }
            }
            else {
              if (length(coordx) != length(coordy)) {
                error <- c("the number of the two coordinates does not match\n")
                return(list(error = error))
              }
              if (length(coordx) != numsite) {
                error <- c("the number of coordinates does not match with the number of spatial observations\n")
                return(list(error = error))
              }
            }
          }
        }
      }
      else {
        if (!is.numeric(coordt)) {
          error <- "insert a numerical vector of temporal coordinates\n"
          return(list(error = error))
        }
        if (length(coordt) <= 1) {
          error <- "insert a numerical vector of temporal coordinates\n"
          return(list(error = error))
        }
        if (replicates > 1) {
          if (grid) {
            if (is.null(dimdata)) {
              error <- c("insert an array d x d x t of n iid spatial-temporal observations\n")
              return(list(error = error))
            }
            if (length(dimdata) != 4) {
              error <- c("the dimension of the data matrix is not correct\n")
              return(list(error = error))
            }
            if (length(coordx) != dimdata[1] || length(coordy) != 
                dimdata[2]) {
              error <- c("the number of coordinates does not match with the number of spatial observations\n")
              return(list(error = error))
            }
            if (length(coordt) != dimdata[3]) {
              error <- c("the number of the temporal coordinate does not match with the third dimensiona of the data matrix\n")
              return(list(error = error))
            }
            if (dimdata[4] != replicates) {
              error <- c("the number of replications doen not match with the data observations\n")
              return(list(error = error))
            }
          }
          else {
            if (is.null(dimdata)) {
              error <- c("insert an array t x d of n iid spatial-temporal observations\n")
              return(list(error = error))
            }
            if (length(dimdata) != 3) {
              error <- c("the dimension of the data matrix is not correct\n")
              return(list(error = error))
            }
            if (is.null(coordy)) {
              if (dimdata[2] != nrow(coordx) || ncol(coordx) != 
                  2) {
                error <- c("the number of coordinates does not match with the number of spatial observations\n")
                return(list(error = error))
              }
            }
            else {
              if (length(coordx) != length(coordy)) {
                error <- c("the number of the two spatial coordinates does not match\n")
                return(list(error = error))
              }
              if (length(coordx) != dimdata[2]) {
                error <- c("the number of coordinates does not match with the number of spatial observations\n")
                return(list(error = error))
              }
            }
            if (dimdata[1] != length(coordt)) {
              error <- c("the time coordinate does not match with the number of rows of the data array\n")
              return(list(error = error))
            }
            if (dimdata[3] != replicates) {
              error <- c("the number of replications does not match with the data observations\n")
              return(list(error = error))
            }
          }
        }
        else {
          if (grid) {
            if (is.null(dimdata)) {
              error <- c("insert an array of d x d x t spatial-temporal observations\n")
              return(list(error = error))
            }
            if (length(dimdata) != 3) {
              error <- c("the dimension of the data matrix is not correct\n")
              return(list(error = error))
            }
            if (length(coordx) != dimdata[1] || length(coordy) != 
                dimdata[2]) {
              error <- c("the number of coordinates does not match with the number of spatial observations\n")
              return(list(error = error))
            }
            if (dimdata[3] != length(coordt)) {
              error <- c("the time coordinate does not match with the third dimension of the data array\n")
              return(list(error = error))
            }
          }
          else {
            if (is.null(dimdata)) {
              error <- c("insert a matrix of t x d spatial-temporal observations\n")
              return(list(error = error))
            }
            if (length(dimdata) != 2) {
              error <- c("the dimension of the data matrix is not correct\n")
              return(list(error = error))
            }
            if (is.null(coordy)) {
              if (dimdata[2] != nrow(coordx) || ncol(coordx) != 
                  2) {
                error <- c("the number of coordinates does not match with the number of spatial observations\n")
                return(list(error = error))
              }
            }
            else {
              if (length(coordx) != length(coordy)) {
                error <- c("the number of the two spatial coordinates does not match\n")
                return(list(error = error))
              }
              if (length(coordx) != dimdata[2]) {
                error <- c("the number of coordinates does not match with the number of spatial observations\n")
                return(list(error = error))
              }
            }
            if (dimdata[1] != length(coordt)) {
              error <- c("the time coordinate does not match with the number of the matrix rows\n")
              return(list(error = error))
            }
          }
        }
      }
    }
    if (fcall == "Simulation") {
      if (type == "Tapering") {
        if (!taper %in% c("Wendland1", "Wendland2", "Wendland3", 
                          "Wendland1_Wendland1", "Wendland1_Wendland2", 
                          "Wendland1_Wendland3", "Wendland2_Wendland1", 
                          "Wendland2_Wendland2", "Wendland2_Wendland3", 
                          "Wendland3_Wendland1", "Wendland3_Wendland2", 
                          "Wendland3_Wendland3", "qt_time", "qt_space")) {
          error <- "insert a correct name for the taper correlation model\n"
          return(list(error = error))
        }
        if (!is.null(maxdist)) {
          error <- "insert a positive numeric value for the maximum spatial distance\n"
          if (!is.numeric(maxdist)) 
            return(list(error = error))
          else if (maxdist < 0) 
            return(list(error = error))
        }
        if (!is.null(maxtime)) {
          error <- "insert a positive numeric value for the maximum temporal distance\n"
          if (!is.numeric(maxtime)) 
            return(list(error = error))
          else if (maxtime < 0) 
            return(list(error = error))
        }
        if (!is.null(tapsep)) {
          error <- "separability parameter of spacetime taper must be between 0  and 1\n"
          if (!is.numeric(tapsep)) 
            return(list(error = error))
          else if (tapsep < 0 || tapsep > 1) 
            return(list(error = error))
        }
      }
      if (is.null(param) || !is.list(param)) {
        error <- "insert the parameters as a list\n"
        return(list(error = error))
      }
      if (length(param) != length(c(unique(c(NuisanceParam("Gaussian"), 
                                             NuisanceParam(model))), CorrelationParam(corrmodel)))) {
        error <- "some parameters are missing or does not match with the declared model\n"
        return(list(error = error))
      }
      if (!all(names(param) %in% c(unique(c(NuisanceParam("Gaussian"), 
                                            NuisanceParam(model))), CorrelationParam(corrmodel),
                                   paste0("beta",1:dim(covariate)[2])))) {
        error <- "some names of the parameters are not correct\n"
        return(list(error = error))
      }
      if (!CheckParamRange(unlist(param))) {
        error <- "some parameters are out of the range\n"
        return(list(error = error))
      }
      if (model %in% c("BrowResn", "ExtT") && !is.null(numblock)) {
        error <- "insert the number of observation in each block\n"
        if (!is.numeric(numblock)) 
          return(list(error = error))
        else if (numblock < 0) 
          return(list(error = error))
      }
      if (model == "BrowResn") {
        if (!corrmodel %in% c("exponential", "gauss", 
                              "gencauchy", "stable", "stable_stable")) {
          error <- "the correlation model is not adequate for the Brown-Resnick process\n"
          return(list(error = error))
        }
        if (corrmodel == "gencauchy" && param["power1"] != 
            2) {
          error <- "the parameter power1 need to be equal to 2 for the Brown-Resnick process\n"
          return(list(error = error))
        }
      }
    }
  }
  else {
    if (class(corrmodel) != "CovMat") {
      error <- "covmatrix is not a object of class covmatrix\n"
      return(list(error = error))
    }
    if (missing(coordx)) {
      error <- "spatial locations must be a matrix of dimension 2\n"
      return(list(error = error))
    }
    else {
      if (is.vector(coordx) && !length(coordx) == 2) {
        error <- "spatial locations must be a vector of dimension 2\n"
        return(list(error = error))
      }
      if (is.matrix(coordx) && !ncol(coordx) == 2) {
        error <- "spatial locations must be  a matrix of dimension 2\n"
        return(list(error = error))
      }
    }
    if ((!is.null(coordt) && !is.numeric(coordt))) {
      error <- "time  must be a vector\n"
      return(list(error = error))
    }
    if (!type %in% c("simple", "ordinary", "Simple", "Ordinary")) {
      error <- "kriging type can be  simple or Ordinary\n"
      return(list(error = error))
    }
    if (missing(data) || !is.numeric(data)) {
      error <- "insert a numeric vector of data\n"
      return(list(error = error))
    }
    if (length(data) != nrow(corrmodel$covmatrix)) {
      error <- "data dimension does not correspond to the covmatrix dimension\n"
      return(list(error = error))
    }
  }
}

