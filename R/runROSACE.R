#' @import cmdstanr
#' @import dplyr
#' @importFrom readr write_tsv
#' @importFrom stringr str_detect
#' @importFrom stats sd quantile
#' @importFrom posterior default_convergence_measures
#' @importFrom tidyr drop_na
#' @importFrom impute impute.knn
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check Stan Setup
#'
#' @param mc.cores Number of cores to use for parallel builds
#' @param install.update whether to update CmdStan
#'
#' @return None
#'
#' @export
#'
CheckStanSetup <- function(mc.cores, install.update = TRUE) {
  check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
  if (install.update) {
    install_cmdstan(cores = mc.cores, quiet = TRUE)
  }
  set_cmdstan_path()
  cmdstan_path()
  cmdstan_version()
}

#' Compile Stan model
#'
#' @param file Stan model file path
#' @param print whether to print model out model
#'
#' @return CmdStanModel object
#'
#' @export
#'
CompileModel <- function(file, print = TRUE) {
  mod <- cmdstanr::cmdstan_model(file)
  if (print) {
    mod$print()
  }
  mod$exe_file()

  return(mod)
}

#' Run Stan MCMC with input list and compiled model
#'
#' @param input list of input data
#' @param mod compiled model
#' @param seed random seed
#' @param refresh number of iterations between progress updates
#'
#' @return CmdStanMCMC object
#'
#' @export
#'
MCMCRunStan <- function(input, mod, seed = 100, refresh = 100) {
  fit <- mod$sample(
    data = input,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    refresh = refresh
  )

  diagnostics <- fit$diagnostic_summary()
  print(diagnostics)

  return(fit)
}

#' Extract MCMC diagnostics from CmdStanMCMC object
#'
#' @param fit CmdStanMCMC object
#' @param sampler whether to get diagnostics summary or
#' extract detailed sampler diagnostics
#'
#' @return list of diagnostics summary or data.frame of sampler diagnostics
#'
#' @export
#'
MCMCDiagnostics <- function(fit, sampler = FALSE) {

  if (!sampler) {
    diags_summary <- fit$diagnostic_summary()
    return(diags_summary)
  } else {
    diags_sampler <- fit$sampler_diagnostics(format = "df")
    return(diags_sampler)
  }

}

#' Extract "functional score" posterior distribution from CmdStanMCMC object
#'
#' @param fit CmdStanMCMC object
#' @param param.key parameter key to extract
#' @param param.post posterior distribution of parameters (optional)
#' If not given, extract from fit object.
#' @param savefile If given, save tsv output to file
#' @param output.lfsr whether to compute lfsr
#'
#' @return data.frame of "functional score" posterior distribution
#'
#' @export
#'
MCMCScoreDf <- function(fit, param.key, param.post, savefile, output.lfsr = TRUE){

  if (missing(param.post)) {
    # warning("No input of 'param.post'. By default calling summary measures
    #         through 'cmdstanr' summary. Might be slow.",
    #         call. = FALSE,
    #         immediate. = TRUE)

    param.name <- fit$metadata()$variables
    param.name <- param.name[stringr::str_detect(param.name,
                                                 paste("^", param.key, "\\[", sep = ""))]
    if (length(param.name) == 0) {
      param.name <- fit$metadata()$variables
      param.name <- param.name[stringr::str_detect(param.name, paste("^", param.key, sep = ""))]
    }
    if (length(param.name) == 0) {
      stop("No parameter found. Check the spelling of param.key.")
    }

    df <-
      fit$summary(param.name, mean, stats::sd,
                  ~stats::quantile(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    colnames(df)[3] <- "sd"

  } else {
    df <- param.post %>%
      dplyr::filter(stringr::str_detect(.data$variable,
                                        paste("^", param.key, "\\[", sep = "")))
    if (nrow(df) == 0) {
      stop("'param.post' does not contain variables with 'param.key'.")
    }
  }

  if (output.lfsr) {
    lfsr <- MCMCLfsr(fit, param.key = param.key)
    df$lfsr <- lfsr
  }

  if (!missing(savefile)) {
    readr::write_tsv(df, file = savefile)
  }

  return(df)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param pos.label vector of true position of variants
#' @param ctrl.label vector of whether variant is control
#' @param thred integer, threshold for number of variants per position label (index)
#' Passed to function "varPosIndexMap"
#‘
#' @rdname GenRosaceInput
#' @method GenRosaceInput AssayGrowth
#' @export
GenRosaceInput.AssayGrowth <- function(object, save.input, pos.label, ctrl.label, thred = 10, ...) {
  CheckDots(...)

  # generate variants - mean count group mapping
  # TODO: change 25 heuristics here
  row_idx <- apply(array(object@norm.var.names), 1, FUN = function(x) {which(x == object@var.names)})
  raw.counts <- object@counts[row_idx, ]
  vMAPm <- ceiling(rank(rowSums(raw.counts, na.rm = TRUE))/25)


  if (!is.na(pos.label[1])) {
    # generate variants - position mapping
    df_map <-
      varPosIndexMap(var.names = object@norm.var.names,
                     pos.label = pos.label,
                     ctrl.label = ctrl.label,
                     thred = thred)
    if (max(df_map$index) != length(unique(df_map$index))) {
      stop("Error when generating position index.")
    }

    # generate input list
    input <- list(m = object@norm.counts,
                  T = object@rounds + 1,
                  t = seq(0, object@rounds)/object@rounds,
                  V = length(object@norm.var.names),
                  vMAPp = df_map$index,
                  P = length(unique(df_map$index)),
                  vMAPm = vMAPm,
                  M = max(vMAPm))

  } else {
    # generate variants "map"
    df_map <- data.frame(variants = object@norm.var.names)

    # generate input list
    input <- list(m = object@norm.counts,
                  T = object@rounds + 1,
                  t = seq(0, object@rounds)/object@rounds,
                  V = length(object@norm.var.names),
                  vMAPm = vMAPm,
                  M = max(vMAPm))
  }

  # save and return
  if (!missing(save.input)) {
    save(input, df_map, file = save.input)
  }
  return(list(input = input, df_map = df_map))

}

#' @param pos.label vector of true position of variants
#' @param ctrl.label vector of whether variant is control
#' @param thred integer, threshold for number of variants per position label (index)
#' Passed to function "varPosIndexMap"
#' @rdname GenRosaceInput
#' @method GenRosaceInput AssaySetGrowth
#' @export
GenRosaceInput.AssaySetGrowth <- function(object, save.input, pos.label, ctrl.label, thred = 10, ...) {
  CheckDots(...)

  # generate variants - mean count group mapping
  # TODO: change 25 heuristics here
  vMAPm <- ceiling(rank(rowSums(object@raw.counts, na.rm = TRUE))/25)

  # impute the combined.counts matrix if not complete
  impute.output <- imputeAssaysCountKNN(object@combined.counts, object@rounds)
  counts <- impute.output$counts
  rounds <- impute.output$rounds

  if (!is.na(pos.label[1])) {
    # generate variants - position mapping
    df_map <-
      varPosIndexMap(var.names = object@var.names,
                    pos.label = pos.label,
                    ctrl.label = ctrl.label,
                    thred = thred)

    if (max(df_map$index) != length(unique(df_map$index))) {
      stop("Error when generating position index.")
    }

    # generate input list
    input <- list(m = counts,
                  T = sum(rounds + 1),
                  t = unlist(lapply(rounds, function(x) seq(0, x)/max(rounds))),
                  V = length(object@var.names),
                  vMAPp = df_map$index,
                  P = length(unique(df_map$index)),
                  vMAPm = vMAPm,
                  M = max(vMAPm))
  } else {
    df_map <- data.frame(variants = object@var.names)

    # generate input list
    input <- list(m = counts,
                  T = sum(rounds + 1),
                  t = unlist(lapply(rounds, function(x) seq(0, x)/max(rounds))),
                  V = length(object@var.names),
                  vMAPm = vMAPm,
                  M = max(vMAPm))
  }

  if (!missing(save.input)) {
    save(input, df_map, file = save.input)
  }
  return(list(input = input, df_map = df_map))

}

#' @rdname MCMCCreateScore
#' @method MCMCCreateScore Assay
#' @export
MCMCCreateScore.Assay <- function(object, main.score,
                                  param.post, diags) { # optional, can be missing

  # score_all <- cbind(var.map, main.score)
  score_all <- main.score

  score <- score_all %>%
    dplyr::select(.data$variants, .data$mean, .data$sd, .data$lfsr)
  optional.score <- score_all %>%
    dplyr::select(-.data$variants, -.data$mean, -.data$sd, -.data$lfsr)

  misc <- list()
  if (!missing(param.post)) {
    misc <- append(misc, list(param.post = param.post))
  }
  if (!missing(diags)) {
    misc <- append(misc, list(diags = diags))
  }
  if (isa(object, "AssayGrowth")) { # for Rosette Object
    misc <- append(misc, list(rounds = object@rounds))
  }

  score <- CreateScoreObject(method = "ROSACE",
                             type = class(object)[1],
                             assay.name = names(object),
                             score = score,
                             optional.score = optional.score,
                             misc = misc)
}

#' @rdname MCMCCreateScore
#' @method MCMCCreateScore AssaySet
#' @export
MCMCCreateScore.AssaySet <- function(object, main.score,
                                     param.post, diags) { # optional, can be missing

  # score_all <- cbind(var.map, main.score)
  score_all <- main.score

  score <- score_all %>% dplyr::select(.data$variants, .data$mean, .data$sd, .data$lfsr)
  optional.score <- score_all %>% dplyr::select(-.data$variants, -.data$mean, -.data$sd, -.data$lfsr)

  misc <- list()
  if (!missing(param.post)) {
    misc <- append(misc, list(param.post = param.post))
  }
  if (!missing(diags)) {
    misc <- append(misc, list(diags = diags))
  }

  score <- CreateScoreObject(method = "ROSACE",
                             type = class(object)[1],
                             assay.name = names(object),
                             score = score,
                             optional.score = optional.score,
                             misc = misc)
}

#' @param pos.label vector of true position of variants
#' @param ctrl.label vector of whether variant is in the control group (NA if none provided)
#'
#' @rdname RunRosace
#' @method RunRosace AssayGrowth
#' @export
#'
RunRosace.AssayGrowth <- function(object, savedir, mc.cores = 4, debug = FALSE, install = TRUE,
                                  pos.label, ctrl.label, ...) {
  CheckDots(..., args = "thred")
  return(helperRunRosaceGrowth(object = object,
                               savedir = savedir,
                               mc.cores = mc.cores,
                               pos.label = pos.label,
                               ctrl.label = ctrl.label,
                               debug = debug,
                               install = install,
                               ...))
}

#' @param pos.label vector of true position of variants
#' @param ctrl.label vector of whether variant is in the control group (NA if none provided)
#'
#' @rdname RunRosace
#' @method RunRosace AssaySetGrowth
#' @export
#'
RunRosace.AssaySetGrowth <- function(object, savedir, mc.cores = 4, debug = FALSE, install = TRUE,
                                     pos.label, ctrl.label, ...) {
  CheckDots(..., args = "thred")
  return(helperRunRosaceGrowth(object = object,
                               savedir = savedir,
                               mc.cores = mc.cores,
                               pos.label = pos.label,
                               ctrl.label = ctrl.label,
                               debug = debug,
                               install = install,
                               ...))
}

#' @param name Name of the object to be analyzed
#' @param type "Assay" or "AssaySet"
#' @param pos.col For Growth screen, the column name for position in the var.data
#' (optional in no_pos mode)
#' @param ctrl.col For Growth screen, optional for control to have one position index
#' @param ctrl.name For Growth screen, optional, the name of the control type
#'
#' @rdname RunRosace
#' @method RunRosace Rosace
#' @export
#'
RunRosace.Rosace <- function(object, savedir, mc.cores = 4, debug = FALSE, install = TRUE,
                             name, type,
                             pos.col, ctrl.col, ctrl.name, ...) {

  # Extract Assay
  if (type == "Assay") {
    sub_object <- ExtractAssay(object, name)
  } else if (type == "AssaySet") {
    sub_object <- ExtractAssaySet(object, name)
  } else {
    stop("Unsupported class type. Provide Assay or AssaySet.")
  }

  # Generate Score Object
  if (isa(sub_object, "AssayGrowth") || isa(sub_object, "AssaySetGrowth")) {
    # Growth screen
    # pos.col optional
    if (missing(pos.col)) {
      warning("position column not provided. run the no-position model.")
      pos.label <- NA
    } else {
      if (type == "Assay") {
        pos.label <- ExtractVarAssay(object, name, norm = TRUE)[[pos.col]]
      } else {
        pos.label <- ExtractVarAssaySet(object, name)[[pos.col]]
      }
    }

    # ctrl.col optional
    if (missing(ctrl.col) || missing(ctrl.name)) {
      warnings("control column or name not provided.")
      ctrl.label <- NA
    } else {
      if (type == "Assay") {
        ctrl.label <- ExtractVarAssay(object, name, norm = TRUE)[[ctrl.col]] == ctrl.name
      } else {
        ctrl.label <- ExtractVarAssaySet(object, name)[[ctrl.col]] == ctrl.name
      }
    }

    # run Rosace
    # AssayGrowth/AssaySetGrowth: pos.label (could be NA), ctrl.label (could be NA), thred (optional)
    score <- RunRosace(object = sub_object,
              savedir = savedir,
              mc.cores = mc.cores,
              debug = debug,
              install = install,
              pos.label = pos.label,
              ctrl.label = ctrl.label,
              ...)
  } else {
    stop("THIS IS VERY WRONG. Check ExtractAssay and ExtractAssaySet.")
  }

  if (debug) {
    fit <- score$fit
    score <- score$score
    return(list(fit = fit, score = score))

  } else {
    # Add Score Object to Rosace if not debugging
    object <- AddScoreData(object = object, score = score)
    return(object)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ...: thred (optional)
helperRunRosaceGrowth <- function(object, savedir, mc.cores, pos.label, ctrl.label,
                                  debug = FALSE, install = TRUE, ...) {
  # create directory if not exists
  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }

  # stan check
  CheckStanSetup(mc.cores = mc.cores, install.update = install)

  # model
  if (is.na(pos.label[1])) {
    mod.file <- WriteStanModel(type = "growth_nopos")
  } else {
    mod.file <- WriteStanModel(type = "growth_pos")
  }
  mod <- CompileModel(file = mod.file, print = FALSE)

  # growth: pos.label, thred in ...
  input.list <- GenRosaceInput(object = object,
                               save.input =
                                 paste(savedir, "/input_", names(object),
                                       ".RData", sep = ""),
                               pos.label = pos.label,
                               ctrl.label = ctrl.label,
                               ...)
  input <- input.list$input
  df_map <- input.list$df_map

  # MCMC sampling
  # WARNING: FIT is a temporary environment!!!
  fit <- MCMCRunStan(input, mod, seed = 100, refresh = 10)

  # MCMC diagnostics
  diags <- MCMCDiagnostics(fit, sampler = FALSE)

  # MCMC score
  main.score <- MCMCScoreDf(fit, param.key = "beta",
                            savefile = paste(savedir, "/beta.tsv", sep = ""),
                            output.lfsr = TRUE)
  main.score <- cbind(df_map, main.score)

  epsilon <-  MCMCScoreDf(fit, param.key = "epsilon2",
                          savefile = paste(savedir, "/epsilon2.tsv", sep = ""),
                          output.lfsr = FALSE)

  if (!is.na(pos.label[1])) {
    sigma <-  MCMCScoreDf(fit, param.key = "sigma2",
                          savefile = paste(savedir, "/sigma2.tsv", sep = ""),
                          output.lfsr = FALSE)
    phi <-  MCMCScoreDf(fit, param.key = "phi",
                        savefile = paste(savedir, "/phi.tsv", sep = ""),
                        output.lfsr = FALSE)

    ##### mapping the phi and sigma to main.score
    phi <- phi %>%
      dplyr::mutate(index = 1:nrow(phi)) %>%
      dplyr::select(.data$index, phi_mean = .data$mean, phi_sd = .data$sd)
    sigma <- sigma %>% dplyr::select(sigma2_mean = .data$mean, sigma2_sd = .data$sd)
    df_pos_index <- cbind(phi, sigma)
    main.score <- main.score %>%
      dplyr::left_join(df_pos_index, by = c("index" = "index"))
  }

  # Create Score Object
  score <- MCMCCreateScore(object = object, main.score = main.score, diags = diags)

  if (debug) {
    return(list(score = score, fit = fit))
  } else {
    return(score)
  }
}

#' Compute lfsr for variables from CmdStanMCMC object
#'
#' @param fit CmdStanMCMC object
#' @param param.key character string to match parameter names (functinal score)
#'
#' @return numeric vector of lfsr for each parameter
#'
MCMCLfsr <- function(fit, param.key = "beta") {

  param.name <- fit$metadata()$variables
  param.name <- param.name[stringr::str_detect(param.name,
                                               paste("^", param.key, sep = ""))]
  param.draw <- fit$draws(variables = param.name, format = "draws_matrix")
  lfsr <- apply(param.draw, MARGIN = 2,
                function(x) {min(mean(x > .Machine$double.eps),
                                 mean(x < .Machine$double.eps))})
  return(lfsr)
}

#' Map each variant to its position label (index)
#'
#' Threshold is tricky to choose.
#'
#' @param var.names character vector of variant names
#' @param pos.label character or numeric vector of true positions
#' @param ctrl.label vector of whether variant is control
#' @param thred integer, threshold for number of variants per position label (index)
#'
#' @return data.frame with columns 'variant', 'pos', 'index'
#'
varPosIndexMap <- function(var.names, pos.label, ctrl.label, thred = 10) {

  df_map <- data.frame(variants = var.names, pos = pos.label)
  n_pos <- df_map %>%
    dplyr::group_by(.data$pos) %>%
    dplyr::summarise(n_pos =  n()) 
  n_pos$index <- 0

  # group position into index
  n_pos <- n_pos[order(n_pos$pos), ] # order by position
  curr_index <- 1
  counter <- 0
  for (i in 1:nrow(n_pos)) {
    n_pos$index[i] <- curr_index
    counter <- counter + n_pos$n_pos[i]
    if (counter >= thred) {
      # Reset counter and increment curr_index
      counter <- 0
      curr_index <- curr_index + 1
    }
  }
  if (counter < thred) {
    n_pos$index[n_pos$index == curr_index] <- curr_index - 1
    counter <- 0
  }
  df_map <- df_map %>% dplyr::left_join(n_pos)

  # map synonymous mutation index
  n_syn_group <- max(n_pos$n_pos) - 1
  if (!is.na(ctrl.label[1])) {
    df_map$ctrl <- ctrl.label

    counter <- 0
    for (i in which(df_map$ctrl)[order(df_map$pos[df_map$ctrl])]) {
      df_map$index[i] <- curr_index
      counter <- counter + 1
      if (counter >= n_syn_group) {
        counter <- 0
        curr_index <- curr_index + 1
      }
    }
    if ((counter < thred) && (sum(df_map$ctrl) >= thred)) {
      df_map$index[df_map$index == curr_index] <- curr_index - 1
    }
    # df_map <- df_map %>% dplyr::mutate(index = ifelse(.data$ctrl, curr_idx, .data$index))
  } else {
    df_map$ctrl <- FALSE
  }

  return(df_map)
}

#' Impute Assays count
#'
#' If replicates are with different number of rounds,
#' truncate counts beyond the miminum number of rounds.
#' Otherwise, impute with mean counts of replicates.
#'
#' @param counts A matrix of counts
#' @param rounds A vector of number of rounds
#'
#' @return A list of imputed counts and rounds
#'
imputeAssaysCount <- function(counts, rounds) {

  if (sum(is.na(counts)) == 0) {
    print("No count imputation needed for combined.assays.")
    return(list(counts = counts,
                rounds = rounds))
  }

  # if replicates are with different number of rounds
  if (length(unique(rounds)) >= 2) {
    warnings("Impute assays counts with different number of rounds.
             Truncate counts beyond miminum number of rounds.
             Ex: [2, 3, 3] to [2, 2, 2].",
             immediate. = TRUE)

    end <- cumsum(rounds + 1)
    n_del <- rounds - min(rounds)
    idx_del <- c()
    for (i in 1:length(rounds)) {
      if (n_del[i] != 0) {
        idx_del <- c(idx_del, (end[i] - n_del[i] + 1):end[i])
      }
    }
    rm(end, n_del)

    counts <- counts[, -idx_del]
    rounds <- rep(min(rounds), length(rounds))

    if (ncol(counts) != (sum(rounds + 1))) {
      stop("imputeAssaysCount: error in count column removal.")
    }
  }

  # impute with new counts and rounds
  miss_idx <- which(rowSums(is.na(counts)) > 0)
  for (i in miss_idx) {
    row <- counts[i, ]
    # row to matrix
    data <- matrix(row, ncol = max(rounds + 1), byrow = TRUE)
    # column mean input
    for(j in 1:ncol(data)){
      data[is.na(data[, j]), j] <- mean(data[, j], na.rm = TRUE)
    }
    # matrix to row
    counts[i, ] <- c(t(data))
  }

  if (sum(is.na(counts))) {
    stop("imputeAssaysCount: error in imputation. NA still exists.")
  }

  return(list(counts = counts,
              rounds = rounds))
}

#' Impute Assays count (KNN Package)
#'
#' @param counts A matrix of counts
#' @param rounds A vector of number of rounds
#'
#' @return A list of imputed counts and rounds
#'
imputeAssaysCountKNN <- function(counts, rounds) {

  # TODO: change row max heristics 0.99 here
  mat_impute <- impute::impute.knn(
    counts,
    k = 10,
    rowmax = 0.99,
    maxp = nrow(counts),
    rng.seed = 362436069
  )
  mat_impute <- mat_impute$data

  return(list(counts = mat_impute,
              rounds = rounds))
}

