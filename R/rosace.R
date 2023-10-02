#' @import dplyr
#' @importFrom methods setClass new 
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Rosace Class
#' @slot assays A list of Assay object
#' @slot assay.sets A list of AssaySet object
#' @slot var.data Contains meta-information about each variant.
#' First column has to contain variants' names
#' @slot scores A list of Score object
#' @slot project.name Name of the project
#' @slot misc A list of other information
#'
#' @rdname Rosace-class
#' @exportClass Rosace
#'
Rosace <- methods::setClass(
  Class = 'Rosace',
  slots = c(
    assays = 'list',
    assay.sets = "list",
    var.data = 'data.frame',
    scores = "list",
    project.name = 'character',
    misc = 'list' # optional list
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get the Names of Specific Object Type
#'
#' Get the names of Assay, AssaySet, Score objects
#'
#' @param object A Rosace object
#' @return The names of all component objects in this Rosace object
#'
#' @rdname RosaceObjectName
#' @export
#'
GetAssayName <- function(object) {
  return(names(object@assays))
}

#' @rdname RosaceObjectName
#' @export
#'
GetAssaySetName <- function(object) {
  return(names(object@assay.sets))
}

#' @rdname RosaceObjectName
#' @export
#'
GetScoreName <- function(object) {
  return(names(object@scores))
}

#' Extract Specific Object Type
#'
#' Get the names of Assay, AssaySet, Score objects
#'
#' @param object A Rosace object
#' @param name The targeted name of Assay/AssaySet/Score object
#' @return The Assay/AssaySet/Score object
#'
#' @rdname RosaceObjectExtract
#' @export
#'
ExtractAssay <- function(object, name) {
  if (!name %in% names(object@assays)) {
    print(names(object@assays))
    stop("Invalid Assay name.")
  }

  return(object@assays[[name]])
}

#' @rdname RosaceObjectExtract
#' @export
#'
ExtractAssaySet <- function(object, name) {
  if (!name %in% names(object@assay.sets)) {
    print(names(object@assay.sets))
    stop("Invalid AssaySet name.")
  }

  return(object@assay.sets[[name]])
}

#' @rdname RosaceObjectExtract
#' @export
#'
ExtractScore <- function(object, name) {
  if (!name %in% names(object@scores)) {
    print(names(object@scores))
    stop("Invalid 'scores' name.")
  }

  return(object@scores[[name]])
}

#' Apply function to all Assay objects with the same key
#'
#' @param object A Rosace object
#' @param key The targeted key of Assay object
#' @param fxn function to be applied to all Assay objects with the same key
#' @param ... Additional arguments to be passed to the function
#' @return Rosace object with updated Assay objects
#'
#' @export
#'
ApplyAssayKey <- function(object, key, fxn, ...) {

  valid.key <- FALSE
  for (name in names(object@assays)) {
    assay <- object@assays[[name]]
    if (key == assay@key) {
      valid.key <- TRUE
      object@assays[[name]] <- fxn(assay, ...)
    }
  }

  if (!valid.key) {
    stop("The key provided does not match any Assay object.")
  }

  return(object)
}

#' Extract the Meta Data for Variants in the Specific Object Type
#'
#' Extract the meta data for variants in Assay, AssaySet, Score object
#'
#' @param object A Rosace object
#' @param name The targeted name of Assay/AssaySet/Score object
#' @param norm (for AssayGrowth only) use the 'var.names' or 'norm.var.names'
#' @return Meta data for variants in Assay, AssaySet, Score object
#'
#' @rdname VardataExtract
#' @export
#'
ExtractVarAssay <- function(object, name, norm = TRUE) {
  assay <- ExtractAssay(object, name)

  if (norm && isa(assay, "AssayGrowth")) {
      df <- data.frame(variants = assay@norm.var.names)
  } else {
    df <- data.frame(variants = assay@var.names)
  }
  
  res <- df %>%
    dplyr::left_join(object@var.data,
              by = c("variants" = colnames(object@var.data)[1]))
  return(res)
}

#' @rdname VardataExtract
#' @export
#'
ExtractVarAssaySet <- function(object, name) {
  assayset <- ExtractAssaySet(object, name)

  df <- data.frame(variants = assayset@var.names)

  res <- df %>%
    dplyr::left_join(object@var.data,
              by = c("variants" = colnames(object@var.data)[1]))

  return(res)
}

#' @rdname VardataExtract
#' @export
#'
ExtractVarScore <- function(object, name) {
  score <- ExtractScore(object, name)

  df <- data.frame(variants = ExtractVarNames(score))
  res <- df %>%
    dplyr::left_join(object@var.data,
              by = c("variants" = colnames(object@var.data)[1]))

  return(res)
}

#' Add Assay Data to the Rosace Object
#'
#' @param object A Rosace object
#' @param assay An Assay object
#' @param var.data Variant's meta data for the Assay object
#' @return The Rosace object with added Assay
#' 
#' @export
#'
AddAssayData <- function(object, assay, var.data) {

  if (!isa(object, "Rosace")) {
    stop("object argument is not a 'Rosace' object")
  }

  if (!isa(assay, "Assay")) {
    stop("assay argument is not an 'Assay' object")
  }

  if (missing(var.data)) {
    var.data <- data.frame(variants = assay@var.names)
  }
  merged.var.data <- rbind(object@var.data, var.data)
  merged.var.data <- dplyr::distinct(merged.var.data)
  if (any(duplicated(merged.var.data[[1]]))) {
    stop("Merging failed for variants' data. Information conflict between original object and new data.")
  }
  object@var.data <- merged.var.data

  # check if the variants in the assay are in the var.data
  if (!all(assay@var.names %in% object@var.data[[1]])) {
    stop("Some variants in the count matrix are no in the first column of 'var.data'")

  }

  assay.list <- list(assay)
  names(assay.list) <- paste(assay@key, assay@rep, sep = "_")
  if (names(assay.list) %in% names(object@assays)) {
    warning("Adding assay that has duplicated name with one of assays in the object. Update.")
    object@assays[[names(assay.list)]] <- assay
  } else {
    object@assays <- append(object@assays, assay.list)
  }

  return(object)
}

#' Add Assayset Object to the Rosace Object
#' 
#' @param object A Rosace object
#' @param assayset A AssaySet object
#' 
#' @return The Rosace object with added AssaySet
#' 
AddAssaySetData <- function(object, assayset) {
  if (!isa(object, "Rosace")) {
    stop("object argument is not a 'Rosace' object")
  }
  if (!isa(assayset, "AssaySet")) {
    stop("assayset argument is not a 'Assayset' object")
  }

  assayset.list <- list(assayset)
  names(assayset.list) <- names(assayset)

  if (names(assayset.list) %in% names(object@assay.sets)) {
    warning(
      "AssaySet object with the same name already exists. Update.",
      call. = FALSE,
      immediate. = TRUE
    )
    object@assay.sets[[names(assayset)]] <- assayset
  } else {
    object@assay.sets <- append(object@assay.sets, assayset.list)
  }

  return(object)
}


#' Add Score Object to the Rosace Object
#' 
#' @param object A Rosace object
#' @param score A Score object
#' 
#' @return The Rosace object with added Score
#' @export
#'
AddScoreData <- function(object, score) {

  if (!isa(object, "Rosace")) {
    stop("object argument is not a 'Rosace' object")
  }
  if (!isa(score, "Score")) {
    stop("score argument is not a 'Score' object")
  }

  score.list <- list(score)
  names(score.list) <- names(score)

  if (names(score.list) %in% names(object@scores)) {
    warning(
      "Score object already exists. Update.",
      call. = FALSE,
      immediate. = TRUE
    )
    object@scores[[names(score.list)]] <- score
  } else {
    object@scores <- append(object@scores, score.list)
  }

  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param type Type of assay to create: "growth"
#'
#' @rdname CreateRosaceObject
#' @method CreateRosaceObject matrix
#' @export
#'
CreateRosaceObject.matrix <- function(
    object,
    var.data,
    project.name = 'RosaceProject',
    type,
    ...
) {

  assay <- CreateAssayObject(counts = object, type = type, ...)

  if (missing(var.data)) {
    var.data <- data.frame(variants = assay@var.names)
  } 

  object <- CreateRosaceObject(
    object = assay,
    var.data = var.data,
    project.name = project.name
  )

  return(object)
}

#' @rdname CreateRosaceObject
#' @method CreateRosaceObject Assay
#' @export
#'
CreateRosaceObject.Assay <- function(
    object,
    var.data,
    project.name = 'RosaceProject',
    ...
) {

  CheckDots(...)

  if (missing(var.data)) {
    var.data <- data.frame(variants = object@var.names)
  } else {
    if (!all(object@var.names %in% var.data[[1]])) {
      stop("Some variants in the count matrix are no in the first column of 'var.data'")
    }
  }

  assay.list <- list(object)
  names(assay.list) <- paste(object@key, object@rep, sep = "_")

  object <- methods::new(
    Class = 'Rosace',
    assays = assay.list,
    var.data = var.data,
    project.name = project.name
  )

  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


