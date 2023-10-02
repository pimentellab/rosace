#' @import dplyr
#' @import ggplot2
#' @importFrom stats as.dist cutree hclust optim rect.hclust
#' @importFrom tidyr pivot_wider unite separate
#' @importFrom philentropy JSD
#' @importFrom grDevices png dev.off
#' @importFrom MASS fitdistr
#' @importFrom mclust Mclust mclustBIC
#' @importFrom compositions fitDirichlet
#'
NULL

# library(dplyr)
# library(ggplot2)
# library(philentropy) # JSD distance
# library(MASS) # fit distribution
# library(mclust) # GMM clustering
# library(compositions) # fit dirichlet

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Hierarchical clustering the mutant in Rosette object
#'
#' Determine the number of mutant group from the plot.
#'
#' @param object An Rosette object
#' @param save.plot If TRUE, save the plot to the directory save.dir
#' @param save.dir The directory to save the plot
#'
#' @export
#'
HclustMutant <- function(object, save.plot = FALSE, save.dir = NULL) {

  dist <- JSDDist(object = object)
  hclust_ <- stats::hclust(dist, method = "complete")

  if (!save.plot) {
    plot(hclust_)
  } else {
    grDevices::png(file.path(save.dir, "hclustCheck.png"),
      width = 240, height = 160, res = 1200, units = 'mm')
    plot(hclust_)
    grDevices::dev.off()
  }

  return(hclust_)
}

#' Generate mutant group label in Rosette object
#'
#' @param object An Rosette object
#' @param Gmut Number of mutant groups
#' @param hclust Input from function HclustMutant. If not given, will call HclustMutant.
#' @param save.plot If TRUE, save the plot to the directory save.dir
#' @param save.dir The directory to save the plot
#'
#' @return An Rosette object with mut.label and Gmut
#'
#' @export
#'
GenMutLabel <- function(object, Gmut, hclust, save.plot = FALSE, save.dir = NULL) {

  if (Gmut == 1) {
    # one group, no need for clustering
    mut <- object@score.df$mut
    mut_group <- tibble(mut = unique(mut), mut_group = "mut1")
  } else {
    if (missing(hclust)) {
      hclust <- HclustMutant(object, save.plot = FALSE)
    }

    # save plot
    if (save.plot) {
      grDevices::png(paste(save.dir, "hclust_k", Gmut, ".png", sep = ""),
        width = 240, height = 160, res = 1200, units = 'mm')
      plot(hclust)
      stats::rect.hclust(hclust, k = Gmut, border = 2:6)
      grDevices::dev.off()
    } else {
      plot(hclust)
      stats::rect.hclust(hclust, k = Gmut, border = 2:6)
    }

    # cut_jsd into variant group table
    cut <- stats::cutree(hclust, k = Gmut)
    mut_group <- dplyr::tibble(mut = names(cut),
                               mut_group = paste("mut", cut, sep = ""))
  }

  # update Gmut, mut.label and score.df in rosette
  object@mut.label <- mut_group
  if (!is.null(object@score.df$mut_group)) {
    object@score.df <- object@score.df %>% dplyr::select(-.data$mut_group)
  }
  object@score.df <- object@score.df %>%
    dplyr::left_join(mut_group, by = c("mut" = "mut"))
  object@Gmut <- Gmut

  # remove score.df if mut_group is not defined
  object@score.df <- object@score.df %>% dplyr::filter(!is.na(mut_group))

  return(object)
}

#' Plot histogram of scores in Rosette object
#'
#' Plot histograms with variant group label and/or mutant group label
#'
#' @param object An Rosette object
#' @param var.group If TRUE, plot histograms with variant group label
#' @param mut.group If TRUE, plot histograms with mutant group label
#'
#' @return A ggplot object
#'
#' @export
#'
PlotScoreHist <- function(object, var.group = FALSE, mut.group = FALSE) {

  df <- object@score.df

  if (!var.group & !mut.group) {
    p <- ggplot(df, aes(.data$score)) +
      geom_histogram(color = "grey", bins = 30) +
      theme_bw()
  } else if (var.group & !mut.group) {
    p <- ggplot(df, aes(.data$score, fill = .data$var_group)) +
      geom_histogram(color = "grey", bins = 30,
                     alpha=0.5, position="identity") +
      theme_bw()
  } else if (!var.group & mut.group) {
    p <- ggplot(df, aes(.data$score, fill = .data$mut_group)) +
      geom_histogram(color = "grey", bins = 30) +
      facet_wrap(vars(.data$mut_group))
  } else {
    p <- ggplot(df, aes(.data$score, fill = .data$var_group)) +
      geom_histogram(color = "grey", bins = 30,
                     alpha=0.5, position="identity") +
      facet_wrap(vars(.data$mut_group))
  }

  return(p)
}

#' Generate variant group label in Rosette object
#'
#' @param object An Rosette object
#' @param Gvar Number of variant groups
#'
#' @return An Rosette object with var.dist and Gvar
#'
#' @export
#'
GenVarLabel <- function(object, Gvar) {
  score <- object@score.df$score
  ctrl <- object@score.df$ctrl

  res <- gmmCluster(score, ctrl, Gvar)
  object@var.dist <- res$var.dist
  object@score.df$var_group <- ifelse(ctrl, "ctrl", res$var.group)

  object@Gvar <- Gvar
  return(object)
}

#' Infer distribution for number of variants for each variant and
#' mutant group label at each position
#'
#' Fit a Dirichlet distribution for number of variants in
#' (mutant group, variant group) at each position
#'
#' @param object A Rosette object
#' @param pos.missing Percentage of variants missing at each position allowed
#'
#' @return A Rosette object with mut.var.alpha
#'
#' @export
#'
PMVCountDist <- function(object, pos.missing = 0.2) {

  if (pos.missing > 1 || pos.missing < 0) {
    stop("Percentage of variants missing at each position is [0, 1].")
  }

  score.df <- object@score.df
  if (is.null(score.df$mut_group) || is.null(score.df$var_group)) {
    stop("Must generate mut.label and var.label before calling the function.")
  }

  M <- length(unique(score.df$mut))
  M_thred <- M * (1 - pos.missing)

  mat <- mutvarCount(pos = score.df$pos,
                     mut_group = score.df$mut_group,
                     var_group = score.df$var_group,
                     M_thred = M_thred)

  options(robust = FALSE)
  fit <- compositions::fitDirichlet(mat)
  alpha <- data.frame(mut_var_group = colnames(mat), alpha = fit$alpha)
  alpha <- alpha %>% tidyr::separate(.data$mut_var_group, into = c("mut_group", "var_group"))

  object@mut.var.alpha <- alpha
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Rosace-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param ctrl.label Boolean vector indicating control variants
#'
#' @rdname PlotDisp
#' @method PlotDisp AssayGrowth
#' @export
#'
PlotDisp.AssayGrowth <- function(object, t, ctrl.label, ...) {
  CheckDots(...)

  mat <- object@counts[ctrl.label, ]

  n_arr <- seq(10, nrow(mat), length.out = 50)
  if (missing(t)) {
    phi_arr <- sapply(n_arr, function(n) inferPerVarPhi(mat, n))
  } else {
    phi_arr <- sapply(n_arr, function(n) inferPerVarPhi(mat, n, t))
  }

  data <- data.frame(n_arr, phi_arr)
  p <- ggplot(data, aes(x = n_arr, y = phi_arr)) +
    geom_line() +
    xlab("number of negative control variants randomly chosen") +
    ylab("per variant phi")

  return(p)
}

#' @param name Name of AssayGrowth
#' @param ctrl.col Name of column in score.df that contains type of variants
#' @param ctrl.name Label of control variants in type.col
#'
#' @rdname PlotDisp
#' @method PlotDisp Rosace
#' @export
#'
# main function to use
PlotDisp.Rosace <- function(object, t, name, ctrl.col, ctrl.name, ...) {
  CheckDots(...)

  var.data <- ExtractVarAssay(object, name = name, norm = FALSE)
  ctrl.label <- (var.data[[ctrl.col]] == ctrl.name)
  assay <- ExtractAssay(object, name = name)

  if (missing(t)) {
    p <- PlotDisp(object = assay, ctrl.label = ctrl.label)
  } else {
    p <- PlotDisp(object = assay, t = t, ctrl.label = ctrl.label)
  }

  return(p)
}

#' @param ctrl.label Boolean vector indicating control variants
#'
#' @rdname EstimateDisp
#' @method EstimateDisp AssayGrowth
#' @export
EstimateDisp.AssayGrowth <- function(object, ctrl.label, ...) {
  CheckDots(...)

  mat <- object@counts[ctrl.label, ]
  return(inferPerVarPhi(mat))
}

#' @param name Name of AssayGrowth
#' @param ctrl.col Name of column in score.df that contains type of variants
#' @param ctrl.name Label of control variants in ctrl.col
#'
#' @rdname EstimateDisp
#' @method EstimateDisp Rosace
#' @export
#'
EstimateDisp.Rosace <- function(object, name, ctrl.col, ctrl.name, ...) {
  CheckDots(...)

  var.data <- ExtractVarAssay(object, name = name, norm = FALSE)
  ctrl.label <- var.data[[ctrl.col]] == ctrl.name
  assay <- ExtractAssay(object, name = name)

  if (!isa(assay, "AssayGrowth")) {
    stop("EstimateDisp only supports AssayGrowth object.")
  }

  return(EstimateDisp(object = assay, ctrl.label = ctrl.label))
}

#' @rdname EstimateDispStart
#' @method EstimateDispStart AssayGrowth
#' @export
EstimateDispStart.AssayGrowth <- function(object, ...) {
  CheckDots(...)

  count0 <- object@counts[1, ]
  return(inferDispLibrary(count0))
}

#' @param name Name of AssayGrowth
#' @rdname EstimateDispStart
#' @method EstimateDispStart Rosace
#' @export
#'
EstimateDispStart.Rosace <- function(object, name, ...) {
  CheckDots(...)

  assay <- ExtractAssay(object, name = name)

  if (!isa(assay, "AssayGrowth")) {
    stop("EstimateDisp only supports AssayGrowth object.")
  }

  return(EstimateDispStart(object = assay))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Compute Empirical JSD distance
JSDDist <- function(object) {
  df_wide <- GetScoreWide(object)

  # remove columns with more than 20% missing
  idx_rm <- (colSums(is.na(df_wide)) > 0.2 * nrow(df_wide))
  if (sum(idx_rm) > 0) {
    warning("removing column", colnames(df_wide)[idx_rm],
      call. = FALSE,  immediate. = TRUE)
    df_wide <- df_wide[, !idx_rm]
  }

  # matrix approximation and counter
  M <- ncol(df_wide) - 1
  mat_approx <- round(df_wide[, 2:(M + 1)], digits = 1)
  counter <- dplyr::tibble(Value = 0)
  for (i in 1:M) {
    tmp <- as.data.frame(table(mat_approx[, i]))
    colnames(tmp) <- c("Value", colnames(mat_approx)[i])
    tmp$Value <- round(as.numeric(as.character(tmp$Value)), digits = 1)
    counter <- dplyr::full_join(counter, tmp)
  }
  counter[is.na(counter)] <- 0
  rm(i, tmp)

  # compute distance using package philentropy
  mat_jsd <-
    philentropy::JSD(t(as.matrix(counter[, 2:(M + 1)])), unit = "log2", est.prob = "empirical")
  colnames(mat_jsd) <- colnames(counter)[-1]
  rownames(mat_jsd) <- colnames(counter)[-1]
  mat_jsd <- stats::as.dist(mat_jsd)

  return(mat_jsd)
}

# Log likelihood function for dirichlet multinomial
dirichlet_multinomial_log_pdf <- function(x, Phi, p) {
  alpha <- Phi * p
  n <- sum(x)
  alpha_0 <- sum(alpha)
  res <- log(n)
  res <- res + lbeta(alpha_0, n)

  for (i in 1:length(x)) {
    x_k <- x[i]
    alpha_k <- alpha[i]
    if (x_k > 0) {
      res <- res - log(x_k)
      res <- res - lbeta(alpha_k, x_k)
    }
  }
  return(res)
}

# Optim wrapper for negative dirichlet multinomial log likelihood
dirichlet_multinomial_negative_ll_optim_wrapper <- function(x, p) {
  func <- function(Phi) {

    if (is.null(ncol(x)) || (ncol(x) == 1)) {
      return(-dirichlet_multinomial_log_pdf(x, Phi, p))
    }
    else {
      sum <- 0
      for (i in 1:ncol(x)) {
        sum = sum - dirichlet_multinomial_log_pdf(x[, i], Phi, p)
      }
      return(sum)
    }
  }
  return(func)
}

# Infer per variant phi (dispersion) for counts
# Dirichlet alpha = phi * p
inferPerVarPhi <- function(counts, sample, t) {

  # add pseudo-count
  counts <- counts + 0.5

  if (missing(sample)) {
    n <- nrow(counts)
    p <- counts[, 1] / sum(counts[, 1])
    if (missing(t)) {
      x <- counts[, -1]
    } else {
      if (t < 2) {
        stop("t must be greater than 1.")
      }
      x <- counts[, t]
    }
  } else {
    n <- as.integer(sample)
    indices <- sample(seq(1, nrow(counts)), size = n)
    p <- counts[indices, 1] / sum(counts[indices, 1])
    if (missing(t)) {
      x <- counts[indices, -1]
    } else {
      if (t < 2) {
        stop("t must be greater than 1.")
      }
      x <- counts[indices, t]
    }
  }
  Phi_init <- 2.5 * n
  res <- stats::optim(Phi_init, dirichlet_multinomial_negative_ll_optim_wrapper(x, p), method = "BFGS")
  return(res$par[1]/n)
}

# count0: vector input (starting count)
inferDispLibrary <- function(count0) {
  n <- length(count0)
  Phi_init <- 2.5 * n
  p <- rep(1/n, n)
  res <- stats::optim(Phi_init, dirichlet_multinomial_negative_ll_optim_wrapper(count0, p), method = "BFGS")
  return(res$par[1]/n)
}

# Gaussian mixture model clustering for variant groups
gmmCluster <- function(score, ctrl, Gvar = 2) {

  score_ctrl <- score[ctrl]
  fit <- MASS::fitdistr(score_ctrl, "normal")
  var.dist <- data.frame(var_group = "ctrl", mean = fit$estimate[1], sd = fit$estimate[2])

  fit <- mclust::Mclust(score, modelNames = "V", G = Gvar)
  var.group <- paste("var", fit$classification, sep = "")
  var.dist <- rbind(var.dist,
                    data.frame(var_group = paste("var", 1:fit$G, sep = ""),
                              mean = fit$parameters$mean,
                              sd = fit$parameters$variance$scale))

  if (Gvar != 2) {
    warning("Gvar is not 2, might generate ill behaviors from GMM. Needs to be fixed.")
  } else {
    df <- data.frame(order = 1:length(score), score = score, var.group = var.group) %>%
      dplyr::arrange(.data$score)
    change_point <- which(df$var.group != c(df$var.group[-1], df$var.group[nrow(df)]))
    if (length(change_point) == 2) {
      # detect the mass and assign the new group
      if (change_point[1] > nrow(df) - change_point[2]) {
        change_point <- change_point[1]
        df$var.group[(change_point + 1):nrow(df)] <- "var2"
      } else {
        change_point <- change_point[2]
        df$var.group[1:change_point] <- "var1"
      }

      # extract the new group
      df <- df %>% dplyr::arrange(.data$order)
      var.group <- df$var.group

      # update the var.dist
      fit1 <- MASS::fitdistr(df$score[df$var.group == "var1"], "normal")
      fit2 <- MASS::fitdistr(df$score[df$var.group == "var2"], "normal")
      var.dist[var.dist$var_group == "var1", 2] <- fit1$estimate[1]
      var.dist[var.dist$var_group == "var1", 3] <- fit1$estimate[2]
      var.dist[var.dist$var_group == "var2", 2] <- fit2$estimate[1]
      var.dist[var.dist$var_group == "var2", 3] <- fit2$estimate[2]

    } else if (length(change_point) > 2) {
      stop("Weird behavior of GMM. Check manually.")
    }
  }

  return(list(var.dist = var.dist, var.group = var.group))
}

# Compute per positions, number of variants in each (mutant, variant) group
mutvarCount <- function(pos, mut_group, var_group, M_thred) {

  df <- data.frame(pos = pos,
                   mut_group = mut_group,
                   var_group = var_group)
  df <- df %>%
    dplyr::filter(var_group != "ctrl") %>%
    dplyr::group_by(.data$pos, .data$mut_group, .data$var_group) %>%
    dplyr::summarise(countPMV = n(), .groups = 'drop') %>%
    dplyr::mutate(countPMV = as.double(.data$countPMV))

  df <- df %>%
    tidyr::unite("mut_var_group", .data$mut_group, .data$var_group) %>%
    dplyr::arrange(.data$mut_var_group) %>%
    tidyr::pivot_wider(names_from = .data$mut_var_group, values_from = .data$countPMV)

  df <- df[rowSums(df[ , -1], na.rm = T) >= M_thred, -1]
  df[is.na(df)] <- 0.01 # pseudocount 0.01

  return(as.matrix(df))
}


