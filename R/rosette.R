#' @import dplyr
#' @importFrom readr write_tsv
#' @importFrom tidyr pivot_wider
#' @importFrom methods setClass new
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Rosette Class
#'
#' A Rosette object is a container for inferring summary statistics from a
#' AssayGrowth of a Rosace experiment.
#'
#' @slot score.df A data.frame with columns score, pos, mut, ctrl, (var_group),
#'  (mut_group)
#' @slot mut.label A data.frame with columns mut, mut_group
#' @slot var.dist A data.frame with columns var_group, mean, sd, count
#' @slot mut.var.alpha A data.frame with columns mut_group, var_group, alpha,
#' count
#' @slot disp Dispersion of raw count (sequencing step)
#' @slot disp.start Dispersion of starting variant library
#' @slot Gvar Number of variant groups
#' @slot Gmut Number of mutant groups
#' @slot rounds Number of rounds
#' @slot project.name A character string indicating the name of the project
#'
#' @exportClass Rosette
#'
Rosette <- methods::setClass(
  Class = 'Rosette',
  slots = c(
    score.df = 'data.frame', # score, pos, mut, ctrl, (var_group, mut_group)
    rounds = "numeric",
    project.name = 'character',
    Gmut = 'numeric',
    mut.label = 'data.frame', # mut, mut_group
    Gvar = 'numeric',
    var.dist = 'data.frame', # var_group, mean, sd, count
    mut.var.alpha = 'data.frame',
    disp = 'numeric',
    disp.start = "numeric"
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Add dispersion to a Rosette object
#' @param object A Rosette object
#' @param disp A numeric value
#' @return A Rosette object with dispersion added
#' @export
#'
AddDisp <- function(object, disp) {
  object@disp <- disp
  return(object)
}

#' Return wide dataframe of score in a Rosette object
#'
#' Row is position and column is mutant. The entry is scores.
#'
#' @param object A Rosette object
#' @return A "wide" data.frame
#'
#' @export
#'
GetScoreWide <- function(object) {
  widedf <- object@score.df %>%
    dplyr::select(.data$score, .data$pos, .data$mut) %>%
    tidyr::pivot_wider(names_from = .data$mut, values_from = .data$score) %>%
    dplyr::arrange(.data$pos)

  return(widedf)
}

#' Generate summary statistics from a Rosette object
#'
#' Including var.dist, mut.var.alpha, and mut.perc.
#'
#' @param object A Rosette object
#' @param save.files A logical value indicating whether to save the output files
#' @param save.dir A directory to save the output files
#'
#' @export
#'
GenerateOutput <- function(object, save.files = TRUE, save.dir) {

  if (nrow(object@mut.var.alpha) == 0) {
    stop("Object has empty slot 'mut.var.alpha'. Run 'PMVCountDist' first.")
  }
  if (nrow(object@var.dist) == 0) {
    stop("Object has empty slot 'var.dist'. Run 'GenVarLabel' first.")
  }
  if (nrow(object@mut.label) == 0) {
    stop("Object has empty slot 'mut.label'. Run 'GenMutLabel' first.")
  }

  # var.dist
  var.count <- object@score.df %>%
    dplyr::group_by(.data$var_group) %>%
    dplyr::summarise(count = n())
  var.dist <- object@var.dist %>%
    dplyr::left_join(var.count)

  # alpha
  mut.var.count <- object@score.df %>%
    dplyr::group_by(.data$mut_group, .data$var_group) %>%
    dplyr::summarise(count = n())
  mut.var.alpha <- object@mut.var.alpha %>%
    dplyr::left_join(mut.var.count)

  # mut.perc
  mut.perc <- object@mut.label %>%
    dplyr::group_by(.data$mut_group) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(perc = .data$count/sum(.data$count))

  # save
  if (save.files) {
    save.dir <- file.path(save.dir, object@project.name)
    if (!dir.exists(save.dir)) {
      dir.create(save.dir, recursive = TRUE)
    }
    readr::write_tsv(var.dist, file.path(save.dir, "var_dist.tsv"))
    readr::write_tsv(mut.var.alpha, file.path(save.dir, "mut_var_alpha.tsv"))
    readr::write_tsv(mut.perc, file.path(save.dir, "mut_perc.tsv"))
  }

  return(list(var.dist = var.dist,
              mut.var.alpha = mut.var.alpha,
              mut.perc = mut.perc))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param score.name A character string indicating the name of the score
#' @param pos.col A character string indicating the name of the position column in var.data
#' @param mut.col A character string indicating the name of the mutant column in var.data
#' @param ctrl.col A character string indicating the name of the control column in var.data
#' @param ctrl.name A character string indicating the label of the control in the control column
#'
#' @rdname CreateRosetteObject
#' @method CreateRosetteObject Rosace
#'
#' @export
#'
CreateRosetteObject.Rosace <- function(object, project.name, score.name,
                                       pos.col, mut.col,
                                       ctrl.col, ctrl.name,
                                       ...) {
  CheckDots(...)

  score <- ExtractScore(object, name = score.name)
  if (score@type != "AssayGrowth") {
    stop("Rosette currently only support assay of type 'AssayGrowth'.")
  }

  # estimate dispersion
  assay <- ExtractAssay(object, name = score@assay.name)
  var.data <- ExtractVarAssay(object, name = score@assay.name, norm = FALSE)
  ctrl.label <- var.data[[ctrl.col]] == ctrl.name
  disp <- EstimateDisp(object = assay, ctrl.label = ctrl.label)
  disp.start <- EstimateDispStart(object = assay)

  # extract score
  var.data <- ExtractVarScore(object, name = score.name)
  pos.label <- var.data[[pos.col]]
  mut.label <- var.data[[mut.col]]
  ctrl.label <- var.data[[ctrl.col]] == ctrl.name

  rosette <- CreateRosetteObject(object = score,
                                 project.name = project.name,
                                 pos.label = pos.label,
                                 mut.label = mut.label,
                                 ctrl.label = ctrl.label,
                                 disp = disp,
                                 disp.start = disp.start)

  return(rosette)
}

#' @param pos.label A vector of position labels
#' @param mut.label A vector of mutant labels
#' @param ctrl.label A vector of control labels (boolean: ctrl or not)
#' @param disp Dispersion of raw count
#' @param disp.start Dispersion of variant library
#'
#' @rdname CreateRosetteObject
#' @method CreateRosetteObject Score
#' @export
#'
CreateRosetteObject.Score <- function(object, project.name, pos.label, mut.label, ctrl.label, disp, disp.start, ...) {
  CheckDots(...)

  if (object@type != "AssayGrowth") {
    stop("Rosette currently only support assay of type 'AssayGrowth'.")
  }

  df <- data.frame(score = object@score[[2]], pos = pos.label, mut = mut.label, ctrl = ctrl.label)
  rounds <- object@misc$rounds

  rosette <- methods::new(
    Class = 'Rosette',
    score.df = df,
    rounds = rounds,
    project.name = project.name,
    disp = disp,
    disp.start = disp.start
  )

  return(rosette)
}
