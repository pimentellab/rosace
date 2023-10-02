#' @importFrom methods setClass new
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Score Class
#' 
#' @slot method A character string of the method used to calculate the score
#' @slot type A character string of the type: Assay, AssaySet
#' @slot assay.name A character string of the Assay/AssaySet name
#' @slot score A data frame of the score: 
#' The first three columns are variants, score, test-statistics.
#' @slot optional.score A data frame of the optional score 
#' (row corresponding to score)
#' @slot misc A list of other information
#' 
#' @exportClass Score
#'
Score <- methods::setClass(
  Class = 'Score',
  slots = c(
    method = 'character', # SLR, ROSACE, ENRICH2
    type = 'character', # AssayGrowth, AssayGrowthSet
    assay.name = "character",
    score = 'data.frame', # variants, score, test-statistics
    optional.score = 'data.frame', 
    misc = 'list' # optional slot
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Extract variants' names from a Score object
#' 
#' @param score A Score object
#' @return A character vector of variants' names
#'
#' @export
#'
ExtractVarNames <- function(score) {
  return(score@score[[1]])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a Score object
#' 
#' @param method A character string of the method used to calculate the score: SLR, ROSACE, ENRICH2
#' @param type A character string of the type: Assay, Assays
#' @param assay.name A character string of the Assay/Assays name
#' @param score A data frame of the score: The first three columns are variants, score, test-statistics.
#' @param optional.score A data frame of the optional score: The first column is variants.
#' @param misc A list of other information
#'
#' @export
#'
CreateScoreObject <- function(
    method,
    type,
    assay.name,
    score,
    optional.score,
    misc
) {

  if (missing(misc)) {
    misc <- list()
  }

  if (missing(optional.score)) {
    optional.score <- data.frame()
  } else {
    if (nrow(score) != nrow(optional.score)){
      stop('Number of variants in input score and optional score are different.')
    }
  }

  score <- methods::new(
    Class = 'Score',
    method = method,
    type = type,
    assay.name = assay.name,
    score = score,
    optional.score = optional.score,
    misc = misc
  )

  return(score)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Return the name of the Score object
#' 
#' @param x A Score object
#' @return A character string of the name of the Score object
#' 
#' @rdname names
#' @export
#'
names.Score <- function(x) {
  return(paste(x@assay.name, x@method, sep = "_"))
}
