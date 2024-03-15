#' @import dplyr
#' @importFrom impute impute.knn
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mat2df <- function(mat, var.names, rep, set.col) {
  df <- data.frame(var.names = var.names, mat)

  if (set.col) {
    colnames(df)[2:ncol(df)] <- paste("rep", rep, colnames(df)[2:ncol(df)], sep = "")
  }

  return(df)
}

matCombine <- function(mat1, mat2,
                       var.names1, var.names2,
                       rep1, rep2,
                       set.col1, set.col2) {

  df1 <- mat2df(mat1, var.names1, rep1, set.col = set.col1)
  df2 <- mat2df(mat2, var.names2, rep2, set.col = set.col2)
  df.join <- dplyr::full_join(df1, df2)

  combined.counts <- as.matrix(df.join[, -1])
  rownames(combined.counts) <- df.join$var.names
  var.names <- df.join$var.names

  return(list(combined.counts = combined.counts,
              var.names = var.names))
}

zeroImpute <- function(mat) {
  mat[is.na(mat)] <- 0
  return(mat)
}

knnImpute <- function(mat, na.rmax = 0.5) {

  if (max(rowSums(is.na(mat))) > (na.rmax * ncol(mat))) {
    stop("Some variants have more than 'na.rmax' missing data. Filter variants first.")
  }

  mat_impute <- impute::impute.knn(
    mat,
    k = 10,
    rowmax = na.rmax,
    maxp = nrow(mat),
    rng.seed = 362436069
  )
  mat_impute <- mat_impute$data

  return(mat_impute)
}

wtNormalize <- function(mat, var.names, wt.var.names, wt.rm) {

  wt_idx <- which(var.names %in% wt.var.names)
  if (length(wt_idx) == 0) {
    stop("No wild-type variants available. Stop wild-type normalization.")
  } else if (length(wt_idx) == 1){
    wt <- mat[wt_idx, ]
  } else {
    wt <- colSums(mat[wt_idx, ])
  }

  if (wt.rm) {
    norm.var.names <- var.names[-wt_idx]
    mat <- mat[-wt_idx, ]
  } else {
    norm.var.names <- var.names
  }

  mat <- log(mat + 1/2)
  wt <- log(wt + 1/2)

  norm.mat <- t(apply(mat, 1, `-`, wt))
  norm.mat <- apply(norm.mat, 2, '-', norm.mat[, 1])
  rownames(norm.mat) <- norm.var.names

  return(norm.mat)
}

totalNormalize <- function(mat, var.names) {
  total <- colSums(mat, na.rm = TRUE)

  total <- log(total + 1/2)
  mat <- log(mat + 1/2)

  norm.mat <- t(apply(mat, 1, `-`, total))
  norm.mat <- apply(norm.mat, 2, '-', norm.mat[, 1])
  rownames(norm.mat) <- var.names

  return(norm.mat)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param assay An object from class 'Assay'
#' @rdname IntegrateData
#' @method IntegrateData Assay
#' @export
IntegrateData.Assay <- function(object, assay, ...) {
  CheckDots(...)

  # AssayGrowth
  if (class(assay) != class(object)) {
    stop("Two assays are not from the same class.")
  }

  # same key and different rep
  if (assay@key != object@key) {
    stop("Two assays do not have the same 'key'.")
  }
  if (assay@rep == object@rep) {
    stop("Two assays have the same 'rep'.")
  }

  # integrate object starts now
  if (isa(object, "AssayGrowth") && isa(assay, "AssayGrowth")) {
    if (length(object@norm.var.names) == 0) {
      stop("Need to normalize the first assay.")
    }
    if (length(assay@norm.var.names) == 0) {
      stop("Need to normalize the second assay.")
    }

    assayset_list <- matCombine(object@norm.counts, assay@norm.counts,
                                object@norm.var.names, assay@norm.var.names,
                                object@rep, assay@rep,
                                TRUE, TRUE)
    assayset_list_raw <- matCombine(object@counts, assay@counts,
                                    object@var.names, assay@var.names,
                                    object@rep, assay@rep,
                                    TRUE, TRUE)
    raw.counts <- assayset_list_raw$combined.counts
    raw.counts <- raw.counts[apply(array(assayset_list$var.names), 1, FUN = function(x) {which(x == assayset_list_raw$var.names)}), ]

    assayset <- CreateAssaySetObject(
      combined.counts = assayset_list$combined.counts,
      raw.counts = raw.counts, 
      var.names = assayset_list$var.names,
      reps = c(object@rep, assay@rep),
      key = object@key,
      type = "growth",
      rounds = c(object@rounds, assay@rounds)
    )

  } else {
    stop("Unmatched or unknown assay class.")
  }

  return(assayset)
}

#' @param assay An object from class 'Assay'
#' @rdname IntegrateData
#' @method IntegrateData AssaySet
#' @export
IntegrateData.AssaySet <- function(object, assay, ...) {
  CheckDots(...)

  # both the same key
  if (object@key != assay@key) {
    stop("Existing 'Assays' object and 'new.assay' do not have the same 'key'.")
  }

  if (assay@rep %in% object@reps) {
    stop("Existing AssaySet already contains the 'rep' of the new assay.")
  }

  # AssayGrowth
  if (isa(object, "AssaySetGrowth") && isa(assay, "AssayGrowth")) {

    if (length(assay@norm.var.names) == 0) {
      stop("Need to normalize the new assay.")
    }

    assayset_list <- matCombine(object@combined.counts, assay@norm.counts,
                                object@var.names, assay@norm.var.names,
                                object@reps, assay@rep,
                                FALSE, TRUE)
    assayset_list_raw <- matCombine(object@raw.counts, assay@counts,
                                    object@var.names, assay@var.names,
                                    object@reps, assay@rep,
                                    FALSE, TRUE)
    raw.counts <- assayset_list_raw$combined.counts
    raw.counts <- raw.counts[apply(array(assayset_list$var.names), 1, FUN = function(x) {which(x == assayset_list_raw$var.names)}), ]

    assayset <- CreateAssaySetObject(
      combined.counts = assayset_list$combined.counts,
      raw.counts = raw.counts,
      var.names = assayset_list$var.names,
      reps = c(object@reps, assay@rep),
      key = object@key,
      type = "growth",
      rounds = c(object@rounds, assay@rounds)
    )
  } else {
    stop("Unmatched or unknown assay class.")
  }

  return(assayset)

}

#' @param key The key for the experiment in the AssaySet object
#' @rdname IntegrateData
#' @method IntegrateData Rosace
#' @export
IntegrateData.Rosace <- function(object, key, ...) {
  CheckDots(...)

  # find assay in the assays with key
  assayset <- AssaySet()
  assay <- Assay()
  for (i in 1:length(object@assays)) {

    if (object@assays[[i]]@key == key) {
      assay <- object@assays[[i]]

      if (length(assayset@key) == 0) {
        assayset <- Assay2AssaySet(assay)
      } else {
        assayset <- IntegrateData(assayset, assay)
      }
    }
  }

  object <- AddAssaySetData(object = object, assayset = assayset)
  return(object)
}

#' @rdname ImputeData
#' @method ImputeData AssayGrowth
#' @export
ImputeData.AssayGrowth <- function(object, impute.method, ...) {

  if (impute.method == "knn") {
    imputed.mat <- knnImpute(mat = object@counts, ...)
  } else if (impute.method == "zero") {
    CheckDots(...)
    imputed.mat <- zeroImpute(mat = object@counts)
  }

  object@counts <- imputed.mat

  if (nrow(object@norm.counts) > 0) {
    warning(
      "clearing out 'norm.counts' and 'norm.var.names' in the assay",
      call. = FALSE,
      immediate. = TRUE)
    object@norm.counts <- matrix(nrow = 0, ncol = 0)
    object@norm.var.names <- logical(0)
  }

  return(object)
}

#' @param name character string of AssayGrowth name
#' @param key key of AssayGrowth
#' @rdname ImputeData
#' @method ImputeData Rosace
#' @export
ImputeData.Rosace <- function(object, impute.method, name, key, ...) {

  if (!missing(name)) {

    if (!missing(key)) {
      stop("Either name or key must be provided. Not both.")
    }
    assay <- ExtractAssay(object = object, name = name)
    object@assays[[name]] <- ImputeData(assay, impute.method = impute.method, ...)

  } else if (!missing(key)) {
    object <- ApplyAssayKey(object = object, key = key,
                            fxn = ImputeData, impute.method = impute.method, ...)
  } else {
    stop("Either name or key must be provided.")
  }

  return(object)
}

#' @param na.rmax The maximum ratio of NA values allowed in a variant
#' @param min.count The minimum number of counts needed for a variant
#'
#' @rdname FilterData
#' @method FilterData AssayGrowth
#' @export
FilterData.AssayGrowth <- function(object, na.rmax = 0.5, min.count = 20, ...) {

  CheckDots(...)

  # na.rmax remove
  idx <- which(rowSums(is.na(object@counts)) > (ncol(object@counts) * na.rmax))
  if (length(idx) > 0) {
    warning(
      paste("filtering ", length(idx), " variants that have more than ",
            na.rmax, " missing data\n", sep = ""),
      call. = FALSE, immediate. = TRUE)

    object@counts <- object@counts[-idx, ]
    object@var.names <- object@var.names[-idx]
  } 

  # min.count remove
  idx <- which(rowSums(object@counts, na.rm = TRUE) < min.count)
  if (length(idx) > 0) {
    warning(
      paste("filtering ", length(idx), " variants that have less than ",
            min.count, " count\n", sep = ""),
      call. = FALSE, immediate. = TRUE)

    object@counts <- object@counts[-idx, ]
    object@var.names <- object@var.names[-idx]
  } 

  if (nrow(object@norm.counts) > 0) {
    warning(
      "clearing out 'norm.counts' and 'norm.var.names' in the assay",
      call. = FALSE,
      immediate. = TRUE)
    object@norm.counts <- matrix(nrow = 0, ncol = 0)
    object@norm.var.names <- logical(0)
  }

  return(object)
}

#' @param name character string of AssayGrowth name
#' @param key key of AssayGrowth
#' @rdname FilterData
#' @method FilterData Rosace
#' @export
FilterData.Rosace <- function(object, name, key, ...) {

  if (!missing(name)) {

    if (!missing(key)) {
      stop("Either name or key must be provided. Not both.")
    }
    assay <- ExtractAssay(object = object, name = name)
    object@assays[[name]] <- FilterData(assay, ...)

  } else if (!missing(key)) {
    object <- ApplyAssayKey(object = object, key = key, fxn = FilterData, ...)
  } else {
    stop("Either name or key must be provided.")
  }

  return(object)

}

#' @param var.names a vector of variant names to be normalized
#' @rdname NormalizeData
#' @method NormalizeData matrix
NormalizeData.matrix <- function(object, normalization.method, var.names, ...) {

  if (missing(var.names)) {
    warning(
      "No variants' names provided, use the row names",
      call. = FALSE,
      immediate. = TRUE
    )
    if (is.null(rownames(object))) {
      stop("No variant names (rows) present in the input matrix")
    }
    var.names = rownames(object)
  }

  if (sum(is.na(object)) > 0) {
    stop("Counts containing NA, impute the data first.")
  }

  if (normalization.method == "wt") {
    normalized.data <- wtNormalize(
      mat = object,
      var.names = var.names,
      ...
    )
  } else if (normalization.method == "total") {
    CheckDots(...)
    normalized.data <- totalNormalize(
      mat = object,
      var.names = var.names
    )
  } else {
    stop("Unknown normalization method: ", normalization.method)
  }

  return(normalized.data)
}

#' @rdname NormalizeData
#' @method NormalizeData AssayGrowth
#' @export
NormalizeData.AssayGrowth <- function(object, normalization.method, ...) {
  normalized.data <- NormalizeData(object = object@counts,
                                   normalization.method = normalization.method,
                                   var.names = object@var.names,
                                   ...)

  object@norm.var.names <- rownames(normalized.data)
  dimnames(normalized.data) <- NULL
  object@norm.counts <- normalized.data

  return(object)
}

#' @param name character string of AssayGrowth name
#' @param key key of AssayGrowth
#'
#' @rdname NormalizeData
#' @method NormalizeData Rosace
#' @export
NormalizeData.Rosace <- function(object, normalization.method, name, key, ...) {

  if (!missing(name)) {

    if (!missing(key)) {
      stop("Either name or key must be provided. Not both.")
    }
    assay <- ExtractAssay(object = object, name = name)
    object@assays[[name]] <- NormalizeData(assay,
                                           normalization.method = normalization.method,
                                           ...)

  } else if (!missing(key)) {
    object <- ApplyAssayKey(object = object, key = key, fxn = NormalizeData,
                            normalization.method = normalization.method, ...)
  } else {
    stop("Either name or key must be provided.")
  }

  return(object)
}

