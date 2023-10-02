#' @importFrom utils argsAnywhere isS3stdGeneric methods
#'
NULL


#' Check the use of the dots
#'
#' Function to check the use of unused arguments passed to \code{...}; this
#' function is designed to be called from another function to see if an
#' argument passed to \code{...} remains unused and alert the user if so.
#' Adapted from Seurat::CheckDots
#'
#' @param ... Arguments passed to other methods
#' @param args A list of arguments to check for unused arguments
#' @param fxns A list of functions to check for unused arguments
#'
#' @return NULL
#'
CheckDots <- function(..., args = NULL, fxns = NULL) {
  args.names <- names(list(...))

  # no arguments passed
  if (length(list(...)) == 0) {
    return(invisible(x = NULL))
  }
  # no named arguments passed
  if (is.null(args.names)) {
    stop("No named arguments passed")
  }

  # check if args are list of characters
  for (a in args) {
    if (!(is.character(a))) {
      stop("CheckDots only works on args that are characters, not ", class(x = f))
    }
  }
  # only one function passed, convert to list
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  # check if fxns are list of characters or functions
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }

  # start checking arguments
  if (is.null(fxns)) {
    if (is.null(args)) {
      # no arguments needed, check for unused arguments
      unused <- args.names
    } else {
      # check if any of the arguments passed is in args
      unused <- Filter(
        f = function(x) { return(!x %in% args)},
        x = args.names
      )
    }
  } else {
    # get arguments of functions in fxns list
    fxn.args <- suppressWarnings(expr = sapply(
      X = fxns,
      FUN = function(x) {
        x <- tryCatch(
          expr = if (isS3stdGeneric(f = x)) {
            as.character(x = methods(generic.function = x)) }
            else { x },
          error = function(...) { return(x)}
        )
        x <- if (is.character(x = x)) {
          sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
        } else if (length(x = x) <= 1) {
          list(x)
        }
        return(sapply(
          X = x,
          FUN = function(f) { return(names(x = formals(fun = f)))},
          simplify = FALSE, USE.NAMES = TRUE
        ))
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    ))
    fxn.args <- unlist(x = fxn.args, recursive = FALSE)

    # check if function passed are valid
    fxn.null <- vapply(
      X = fxn.args,
      FUN = is.null,
      FUN.VALUE = logical(length = 1L)
    )
    if (all(fxn.null) && !is.null(x = fxns)) {
      stop("None of the functions passed could be found", call. = FALSE)
    } else if (any(fxn.null)) {
      warning(
        "The following functions passed could not be found: ",
        paste(names(x = which(x = fxn.null)), collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
      fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
    }

    # check if any of the functions passed accept the dots
    dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
    names(x = dfxns) <- names(x = fxn.args)
    for (i in 1:length(x = fxn.args)) {
      dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
    }

    # if any of the functions passed accept the dots, get their names
    if (any(dfxns)) {
      dfxns <- names(x = which(x = dfxns)) # get names of functions that accept the dots
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(dfxns, collapse = ', ')
      )
      return(invisible(x = NULL))
    } else {
      # if none of the functions passed accept the dots, check for unused arguments
      unused <- Filter(
        f = function(x) {
          return(!x %in% c(unlist(fxn.args), args))
        },
        x = args.names
      )
    }
  }

  # If there is unused parameters, warn the user
  if (length(x = unused) > 0) {
    msg <- paste0(
      "The following arguments are not used: ",
      paste(unused, collapse = ', ')
    )
    switch(
      EXPR = getOption(x = "Rosace.checkdots", default = 'warn'),
      "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
      "stop" = stop(msg),
      "silent" = NULL,
      stop("Invalid Rosace.checkdots option. Please choose one of warn, stop, silent")
    )
  }
  return(invisible(x = NULL))
}
