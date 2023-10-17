#' @importFrom methods setClass new
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Assay Class
#' @slot counts A matrix of raw counts
#' (variants x time points in growth screen)
#' @slot var.names A vector of variant names
#' @slot key A character string for the experiment
#' @slot rep The identity of the replicate
#'
#' @rdname Assay-class
#' @exportClass Assay
#'
Assay <- methods::setClass(
  Class = 'Assay',
  slots = c(
    counts = 'matrix',
    var.names = 'vector', # vector of character
    key = 'character',
    rep = 'numeric'
  )
)

#' The AssayGrowth Class
#' @slot norm.counts A matrix of normalized counts (variants x time points)
#' @slot norm.var.names A vector of variant names for normalized counts
#' @slot rounds Number of rounds = Number of time points - 1
#'
#' @rdname Assay-class
#' @exportClass AssayGrowth
#'
AssayGrowth <- methods::setClass(
  Class = 'AssayGrowth',
  slots = c(
    norm.counts = 'matrix',
    norm.var.names = "vector", # vector of character
    rounds = "numeric"
  ),
  contains = "Assay"
)

#' The AssaySet Class
#' @slot combined.counts A matrix of counts (variants x cols in reps)
#' @slot var.names A vector of variant names corresponding
#' to the row of combined.counts.
#' @slot reps A vector of replicate identity
#' @slot key A character string for the experiment
#'
#' @rdname AssaySet-class
#' @exportClass AssaySet
#'
AssaySet <- methods::setClass(
  Class = 'AssaySet',
  slots = c(
    combined.counts = 'matrix',
    var.names = 'vector',  # vector of character
    reps = "vector", # vector of numeric
    key = 'character'
  )
)

#' The AssaySetGrowth Class
#' @slot raw.counts A matrix of raw counts (variants x cols in reps)
#' @slot rounds A vector of number of rounds
#'
#' @rdname AssaySet-class
#' @exportClass AssaySetGrowth
#'
AssaySetGrowth <- methods::setClass(
  Class = 'AssaySetGrowth',
  slots = c(
    raw.counts = 'matrix',
    rounds = "vector"  # vector of numeric
  ),
  contains = "AssaySet"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create an Assay object
#'
#' Create an Assay object from a raw count matrix.
#'
#' Non-unique variant names are not allowed. Please make unique before
#' calling this function. By default, variants with more than 0.5 missing (NA)
#' values will be filtered out. You can change the percentage by specifying
#' `na.rmax` or call `FilterData` after creating the rosace object.
#'
#' @param counts Unnormalized data such as raw counts
#' @param var.names A vector of variant names
#' @param key A character string for the experiment
#' @param rep Identity of the replicate
#' @param type Type of experiment: growth
#' @param ... Additional arguments for create assay object.
#' @return An Assay object
#'
#' @rdname CreateAssayObject
#' @export
#'
CreateAssayObject <- function(counts, var.names, key, rep, type, ...) {

  # Check var.names
  if (missing(var.names)) {
    warning("No variants' names provided, use the row names",
            call. = FALSE, immediate. = TRUE)
    if (is.null(rownames(counts))) {
      stop("No variant names (rows) present in the input matrix")
    }
    var.names <- rownames(counts)
  } else {
    if (length(var.names) != nrow(counts)) {
      stop("Length of variant names is not equal to number of rows in input matrix")
    }
    if (!is.null(rownames(counts))) {
      if (any(var.names != rownames(counts))) {
        warning("'var.names' does not match with the row names of the input matrix",
                call. = FALSE, immediate. = TRUE)
      }
    }
  }

  # Check duplicated var.names
  if (anyDuplicated(var.names)) {
    warning(
      "Non-unique variants (rownames) present in the input matrix, making unique",
      call. = FALSE,
      immediate. = TRUE
    )
    var.names <- make.unique(names = var.names)
  }

  # create new object
  if (type == "growth") {
    assay <- CreateAssayGrowthObject(counts = counts, var.names = var.names,
                                     key = key, rep = rep, ...)
  } else {
    stop("Currently not support assay type: ", type)
  }

  return(assay)
}


#' Create an AssayGrowth object
#'
#' The expected format of the input matrix is variants x time points.
#'
#' @param na.rmax Maximum ratio of NA values allowed for a variant
#'  to be included in the assay (AssayGrowth)
#' @return An AssayGrowth object
#'
#' @rdname CreateAssayObject
#'
CreateAssayGrowthObject <- function(counts, var.names, key, rep, na.rmax = 0.5) {
  dimnames(counts) <- NULL # no dimension allowed for slot count

  # create assay
  assay <- methods::new(
    Class = 'AssayGrowth',
    counts = counts,
    var.names = var.names,
    key = key,
    rep = rep,
    rounds = ncol(counts) - 1
  )

  # filter data
  assay <- FilterData(object = assay, na.rmax = na.rmax)

  return(assay)
}

#' Create an AssaySet object
#'
#' @param combined.counts A matrix of counts (variants x time points in all reps)
#' @param var.names A vector of variant names
#' @param reps A vector of replicate identity
#' @param key A character string for the experiment
#' @param type Type of experiment: growth
#' @param ... Additional arguments for create AssaySet object
#' @return An AssaySet object
#'
#' @rdname CreateAssaySetObject
#'
CreateAssaySetObject <- function(combined.counts, var.names, reps, key, type, ...) {

  if (length(var.names) != nrow(combined.counts)) {
    stop("length of 'var.names' does not match with rows of count matrix")
  }
  if (anyDuplicated(var.names)) {
    stop("duplicated variants' name")
  }

  if (type == "growth") {
    assayset <- CreateAssaySetGrowthObject(combined.counts = combined.counts,
                                           var.names = var.names, reps = reps,
                                           key = key, ...)
  } else {
    stop("Currently not support assay type: ", type)
  }

  return(assayset)
}

#' Create an AssaySetGrowth object
#'
#' @param raw.counts raw.counts A matrix of raw counts (variants x cols in reps)
#' @param rounds A vector of number of rounds
#' @return An AssaySetGrowth object
#'
#' @rdname CreateAssaySetObject
#'
CreateAssaySetGrowthObject <- function(combined.counts, var.names, reps, key, raw.counts, rounds) {
  rownames(combined.counts) <- NULL # do not store rownames in counts

  if (sum(rounds + 1) != ncol(combined.counts)) {
    stop("length of 'rounds' does not match with columns of count matrix")
  }

  if (length(rounds) != length(reps)) {
    stop("length of 'rep.label' does not match with length of 'time.label'")
  }

  rownames(combined.counts) <- NULL # do not store rownames in counts
  rownames(raw.counts) <- NULL # do not store rownames in counts

  # create AssaySetGrowth object
  assayset <- methods::new(
    Class = 'AssaySetGrowth',
    combined.counts = combined.counts,
    raw.counts = raw.counts,
    var.names = var.names,
    rounds = rounds,
    reps = reps,
    key = key
  )

  return(assayset)
}

#' Transform an Assay object to AssaySet object
#'
#' @param assay An Assay object
#' @return An AssaySet object
#'
Assay2AssaySet <- function(assay) {

  # check assay type
  if (isa(assay, "AssayGrowth")) {

    if (length(assay@norm.var.names) == 0) {
      stop("Need to normalize the first assay.")
    }

    combined.counts <- assay@norm.counts 
    colnames(combined.counts) <- paste("rep", assay@rep, "X", 1:(assay@rounds+1), sep = "")
    rownames(combined.counts) <- NULL

    raw.counts <- assay@counts[apply(array(assay@norm.var.names), 1, FUN = function(x) {which(x == assay@var.names)}), ]
    colnames(raw.counts) <- paste("rep", assay@rep, "X", 1:(assay@rounds+1), sep = "")
    rownames(raw.counts) <- NULL

    assayset <- CreateAssaySetGrowthObject(combined.counts = combined.counts,
                                           raw.counts = raw.counts,
                                           var.names = assay@norm.var.names,
                                           reps = c(assay@rep),
                                           key = assay@key,
                                           rounds = c(assay@rounds))
  } else {
    stop("This is very wrong LOL.")
  }

  return(assayset)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Return the name of the Assay object
#' @param x An Assay object
#' @return A character string of the name of the Assay object
#'
#' @rdname names
#' @export
#'
names.Assay <- function(x) {
  return(paste(x@key, x@rep, sep = "_"))
}

#' Return the name of the AssaySet object
#' @param x An AssaySet object
#' @return A character string of the name of the AssaySetobject
#'
#' @rdname names
#' @export
#'
names.AssaySet <- function(x) {
  return(x@key)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

