#' @import dplyr
#' @importFrom methods setClass new
#' @importFrom readr write_tsv read_tsv
#' @importFrom stats pnorm
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
#' The first four columns are variants, score, sd, test-statistics.
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
    score = 'data.frame', # variants, score, sd, test-statistics
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
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname OutputScore
#' @method OutputScore Score
#' @export
OutputScore.Score <- function(object, pos.info = FALSE, sig.test = 0.05, ...) {
  CheckDots(...)

  # for Rosace method only
  # process lfsr since the beta is normally distributed
  # warning: if beta is not normally distributed, please don't do this trasnformation
  if (object@method == "ROSACE") {
    df <- object@score %>% dplyr::rowwise() %>%
      mutate(lfsr.neg = stats::pnorm(0, mean = .data$mean, sd = .data$sd, lower.tail = FALSE),
             lfsr.pos = stats::pnorm(0, mean = .data$mean, sd = .data$sd, lower.tail = TRUE),
             lfsr = min(.data$lfsr.neg, .data$lfsr.pos),
             test.neg = (.data$lfsr.neg <= sig.test/2),
             test.pos = (.data$lfsr.pos <= sig.test/2),
             label = ifelse(.data$test.neg, "Neg", "Neutral"),
             label = ifelse(.data$test.pos, "Pos", .data$label)) %>%
      dplyr::ungroup()

    if (pos.info == TRUE) {
      df_pos <- cbind(df,
                      object@optional.score %>%
                        dplyr::select(.data$ctrl, .data$pos,
                                      .data$phi_mean, .data$phi_sd,
                                      .data$sigma2_mean, .data$sigma2_sd))
      df_pos <- df_pos %>%
        dplyr::filter(.data$ctrl == FALSE) %>%
        dplyr::select(.data$pos, .data$phi_mean, .data$phi_sd, .data$sigma2_mean, .data$sigma2_sd) %>%
        unique() %>%
        arrange(.data$pos)

      df_pos <- df_pos %>% dplyr::rowwise() %>%
        mutate(lfsr.neg = stats::pnorm(0, mean = .data$phi_mean, sd = .data$sigma2_mean, lower.tail = FALSE),
               lfsr.pos = stats::pnorm(0, mean = .data$phi_mean, sd = .data$sigma2_mean, lower.tail = TRUE),
               lfsr = min(.data$lfsr.neg, .data$lfsr.pos),
               test.neg = (.data$lfsr.neg <= sig.test/2),
               test.pos = (.data$lfsr.pos <= sig.test/2),
               label = ifelse(.data$test.neg, "Neg", "Neutral"),
               label = ifelse(.data$test.pos, "Pos", .data$label)) %>%
        dplyr::ungroup()

      df <- cbind(df,
                  object@optional.score %>%
                    dplyr::select(.data$phi_mean, .data$phi_sd,
                                  .data$sigma2_mean, .data$sigma2_sd))

      return(list(df_variant = df,
                  df_position = df_pos))
    }

  } else {
    df <- object@score
    warning("OutputScore only supports hypothesis testing labeling for Rosace method.
            Return the score data frame.")
  }

  return(df)
}

#' @param name The targeted name of Score object
#' @rdname OutputScore
#' @method OutputScore Rosace
#' @export
OutputScore.Rosace <- function(object, pos.info = FALSE, sig.test = 0.05, name, ...) {
  CheckDots(...)
  score <- ExtractScore(object, name)
  var <- ExtractVarScore(object, name)
  if (pos.info == FALSE) {
    score <- OutputScore(score, pos.info = pos.info, sig.test = sig.test)
    df <- cbind(var, score[, -1])
    return(df)
  } else {
    score_list <- OutputScore(score, pos.info = pos.info, sig.test = sig.test)
    score_list$df_variant <- cbind(var, score_list$df_variant[, -1])
    return(score_list)
  }
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
