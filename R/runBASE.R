#' @importFrom stats lm
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Internal function for simple linear regression
#'
#' @param y numeric vector of response
#' @param x numeric vector of predictor
#'
#' @return numeric vector of regression coefficients
#'
SLR <- function(y, x) {
  lm_fit <- stats::lm(y ~ x)
  return(summary(lm_fit)$coefficients[2,])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname runSLR
#' @method runSLR AssayGrowth
#' @export
runSLR.AssayGrowth <- function(object, ...) {
  CheckDots(...)

  if (length(object@norm.counts) == 0) {
      stop("Run normalization before regression.")
  }

  rounds <- object@rounds
  xt <- seq(0,rounds)/rounds

  result <- t(apply(object@norm.counts, 1, SLR, xt))
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')

  scores <- data.frame(variant = object@norm.var.names,
                       result[, c("Estimate", "Std. Error", "Pr(>|t|)")],
                       row.names = NULL)
  colnames(scores) <- c('variants', 'estimate', 'std.error', 'p-value')

  optional_scores <- data.frame(result[, c("t value")],
                                row.names = NULL)
  colnames(optional_scores) <- c('t-stats')

  score <- CreateScoreObject(method = "SLR",
                             type = "AssayGrowth",
                             assay.name = names(object),
                             score = scores,
                             optional.score = optional_scores,
                             misc = list(rounds = rounds))

  return(score)
}

#' @rdname runSLR
#' @method runSLR AssaySetGrowth
#' @export
runSLR.AssaySetGrowth <- function(object, ...) {
  CheckDots(...)

  denom <- max(object@rounds)
  xt <- unlist(lapply(object@rounds, function(x) seq(0, x)/denom))

  result <- data.frame(t(apply(object@combined.counts, 1, SLR, xt)))
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')

  scores <- data.frame(variant = object@var.names,
                       result[, c("Estimate", "Std. Error", "Pr(>|t|)")],
                       row.names = NULL)
  colnames(scores) <- c('variants', 'estimate', 'std.error', 'p-value')

  optional_scores <- data.frame(result[, c("t value")],
                                row.names = NULL)
  colnames(optional_scores) <- c('t-stats')

  # label which replicate(s) the regression is performed on
  cols <- cumsum(unlist(lapply(object@rounds, function(x) x+1)))
  label <- apply(object@combined.counts, 1,
                 function(row) object@reps[-which(is.na(row[cols]))])

  converted_label <- lapply(label, function(x) {
    if (length(x) > 0) {
      paste('rep', x, collapse = ",")
    } else {
      "all"
    }})

  # scores$which_rep <- unlist(converted_label)
  optional_scores$which_rep <- unlist(converted_label)

  score <- CreateScoreObject(method = "SLR",
                             type = "AssaySetGrowth",
                             assay.name = names(object),
                             score = scores,
                             optional.score = optional_scores,
                             misc = list(max.rounds = denom))

  return(score)
}

#' @param name character string of Assay/AssaySet name
#' @param type Assay or AssaySet
#'
#' @rdname runSLR
#' @method runSLR Rosace
#' @export
runSLR.Rosace <- function(object, name, type, ...) {
  CheckDots(...)

  if (type == "Assay") {
    score <- runSLR(ExtractAssay(object, name))
  } else if (type == "AssaySet") {
    score <- runSLR(ExtractAssaySet(object, name))
  } else {
    stop("Unsupported class type. Provide Assay or AssaySet.")
  }

  object <- AddScoreData(object = object, score = score)

  return(object)

}

