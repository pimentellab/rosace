NULL

########## rosace.R ##########

#' Create a Rosace object
#'
#' Create a Rosace object from raw data
#'
#' @param object Either a matrix-like object with
#' unnormalized raw count or an Assay object
#' @param var.data Variants's meta data. If empty,
#' create a data frame with variants' name only.
#' @param project.name Name of the Rosace object
#' @param ... Arguments passed to CreateAssayObject
#'
#' @return A Rosace object
#'
#' @rdname CreateRosaceObject
#' @export
#'
CreateRosaceObject <- function(object,
                               var.data,
                               project.name = 'RosaceProject',
                               ...) {
  UseMethod(generic = 'CreateRosaceObject', object = object)
}

########## preprocessing.R ##########

#' Normalize count (AssayGrowth)
#'
#' Normalize (filtered and imputed) count data (in Assay)
#' and get normalized count
#'
#' @param object An object: matrix, AssayGrowth, or Rosace
#' @param normalization.method The normalization method to use: "wt" or "total"
#' @param ... Additional arguments to be passed to the normalization method.
#'
#' @return An object with normalized count
#'
#' @rdname NormalizeData
#' @export
#'
NormalizeData <- function(object, normalization.method, ...) {
  UseMethod(generic = 'NormalizeData', object = object)
}

#' Filter out variants in raw count data
#'
#' @param object An object: AssayGrowth (na.rm), or Rosace
#' @param ... Additional arguments to be passed to the filtering method.
#'
#' @return An object with filtered count
#'
#' @rdname FilterData
#' @export
#'
FilterData <- function(object, ...) {
  UseMethod(generic = 'FilterData', object = object)
}

#' Impute raw count data and get imputed count
#'
#'
#' @param object An object: matrix (AssayGrowth), AssayGrowth, or Rosace
#' @param impute.method The imputation method to use: "knn" or "zero"
#' @param ... Additional arguments to be passed to the imputation method.
#'
#' @return An object with imputed count
#'
#' @rdname ImputeData
#' @export
#'
ImputeData <- function(object, impute.method, ...) {
  UseMethod(generic = 'ImputeData', object = object)
}

#' Integrate Assay into Assays
#'
#' @param object An object: Assay, AssaySet, or Rosace
#' @param ... Additional arguments to be passed to the integration method.
#'
#' @return An object with integrated count (AssaySet)
#'
#' @rdname IntegrateData
#' @export
#'
IntegrateData <- function(object, ...) {
  UseMethod(generic = 'IntegrateData', object = object)
}

########## runBASE.R ##########
#' Do simple linear regression on normalized count  (Growth)
#'
#' @param object An object: AssayGrowth, AssaySetGrowth, or Rosace
#' @param ... Additional arguments to be passed to the regression method.
#'
#' @return A Score object with regression result
#'
#' @rdname runSLR
#' @export
#'
runSLR <- function(object, ...) {
  UseMethod(generic = 'runSLR', object = object)
}


##########  runstan.R
#' Generate input for Rosace model from Assay or AssaySet
#'
#' @param object Assay or AssaySet
#' @param save.input if given save input to the save.input path
#' @param ... Additional arguments to be passed to the input generation method.
#'
#' @return A list of input for Rosace1 model
#'
#' @rdname GenRosaceInput
#' @export
#'
GenRosaceInput <- function(object, save.input, ...) {
  UseMethod(generic = 'GenRosaceInput', object = object)
}

#' Create score object from MCMC result
#'
#' @param object Assay or AssaySet
#' @param main.score functional score
#' @param var.map variant map
#' @param param.post (optional) posterior of parameters other than "functional score"
#' @param diags (optional) diagnostics of MCMC
#'
#' @return A Score object
#'
#' @rdname MCMCCreateScore
#' @export
#'
MCMCCreateScore <- function(object, main.score, var.map, param.post, diags) {
  UseMethod(generic = 'MCMCCreateScore', object = object)
}

# TODO: Build DEBUG Option
#' Run Rosace on an Assay/AssaySet object
#'
#' @param object Rosace, Assay/AssaySet (default)
#' @param savedir directory to save the output
#' @param mc.cores integer, number of cores to use for parallel computing
#' @param debug logical, if TRUE, return a list of cmdstanfit and Score object.
#' @param install logical, if TRUE, install or update cmdstanr
#' There's a chance for Github download error if running multiple instances on a server.
#' @param ... Additional arguments to be passed to the run method.
#'
#' @return A object with Rosace result (Score object) if debug = FALSE
#'
#' @rdname RunRosace
#' @export
#'
RunRosace <- function(object, savedir, mc.cores, debug, install, ...) {
  UseMethod(generic = 'RunRosace', object = object)
}

##########  rosette.R ##########
#' Create a Rosette object
#'
#' @param object An object: numeric, Score, or Rosace
#' @param project.name Name of the Rosette object
#' @param ... Additional arguments to be passed to the creation method.
#'
#' @return A Rosette object
#'
#' @rdname CreateRosetteObject
#' @export
#'
CreateRosetteObject <- function(object, project.name, ...) {
  UseMethod(generic = 'CreateRosetteObject', object = object)
}


########## simsummary.R ##########
#' Raw count dispersion plot
#'
#' @param object An object: matrix, Assay, or Rosace
#' @param t The time point to compute dispersion (>= 2). If missing, use all time points.
#' @param ... Additional arguments to be passed to the plot method.
#'
#' @rdname PlotDisp
#' @export
#'
PlotDisp <- function(object, t, ...) {
  UseMethod(generic = 'PlotDisp', object = object)
}

#' Estimate dispersion from raw count (use all time points)
#'
#' @param object An object: AssayGrowth, or Rosace
#' @param ... Additional arguments to be passed to the estimation method.
#'
#' @rdname EstimateDisp
#' @export
#'
EstimateDisp <- function(object, ...) {
  UseMethod(generic = 'EstimateDisp', object = object)
}

#' Estimate dispersion from count0
#'
#' @param object An object: AssayGrowth, or Rosace
#' @param ... Additional arguments to be passed to the estimation method.
#'
#' @rdname EstimateDispStart
#' @export
#'
EstimateDispStart <- function(object, ...) {
  UseMethod(generic = 'EstimateDispStart', object = object)
}

