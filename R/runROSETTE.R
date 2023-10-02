#' @import dplyr
#' @importFrom tidyr separate
#' @importFrom readr write_tsv
#' @importFrom stats rmultinom rnorm rpois rbinom
#' @importFrom compositions rDirichlet.rcomp
#' @importFrom rjson toJSON
#' @importFrom utils write.table
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Simulate data from a config file
#'
#' @param config A list of config
#' @param save.tsv A boolean indicating whether to save the data in tsv format
#' @param save.rosace A boolean indicating whether to save the Rosace object
#' @param save.enrich2 A boolean indicating whether to save the data for enrich2
#'
#' @return None
#'
#' @export
#'
runRosette <- function(config,
                       save.tsv = TRUE,
                       save.rosace = TRUE,
                       save.enrich2 = FALSE){

  # create save directory
  save.dir <- config[['sim']][['save.sim']]
  save.dir <- file.path(save.dir,
                        paste(config[["sim"]][["type.sim"]],
                              "_rep",  config[["exp"]][["n.rep"]],
                              "_rd", config[["exp"]][["n.round"]],
                              "_", config[["sim"]][["mode.sim"]],
                              sep = ""))
  if (!dir.exists(save.dir)) {
    dir.create(save.dir, recursive = TRUE)
  }

  # start simulation
  sim_rounds <- config[["sim"]][["n.sim"]]
  for(i in 1:sim_rounds){
    set.seed(i)
    print(paste("Starting simulation round", i))

    effects_list <- GenerateEffect(cfg = config) # effects, expected_effects, expected_group
    counts_list <- GenerateCount(cfg = config, effects = effects_list$effects) # counts, sequencing

    sub.save.dir <- file.path(save.dir, paste("sim", i, sep = ""))
    if (!dir.exists(sub.save.dir)) {
      dir.create(sub.save.dir, recursive = TRUE)
    }

    output_tsv(effects_list = effects_list, counts_list = counts_list, save.dir = sub.save.dir)
    output_rosace(effects_list = effects_list, counts_list = counts_list,
                  Nrep = config[["exp"]][["n.rep"]], Nround = config[["exp"]][["n.round"]],
                  save.dir = sub.save.dir)
    output_enrich2(counts_list = counts_list,
                   Nrep = config[["exp"]][["n.rep"]], Nround = config[["exp"]][["n.round"]],
                   mode.sim = config[['sim']][["mode.sim"]], save.dir = sub.save.dir)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GenerateEffect <- function(cfg) {
  Nround <- cfg[["exp"]][["n.round"]]
  Nrep <- cfg[["exp"]][["n.rep"]]
  Npos <- cfg[["exp"]][["n.pos"]]
  Nmut <- cfg[["exp"]][["n.mut"]]
  type_exp <- cfg[["sim"]][["type.sim"]]
  Gvar <- length(cfg[["effect"]][["var.dist"]][["var_group"]])
  Gmut <- length(cfg[["effect"]][["mut.perc"]][['mut_group']])

  # expected effect distribution
  expected_effect_list <-
    expected_effect_distribution_list(var_para = cfg[["effect"]][["var.dist"]],
                                      var_shrink_normal = cfg[["effect"]][["var.shrink"]],
                                      effect_scale = Nround/cfg[["effect"]][["rounds"]],
                                      null_group = cfg[["effect"]][["null.var.group"]])

  # generate expected effect
  expected_result <- generate_expected_effect(expected_effect_list = expected_effect_list,
                                              mut_para = cfg[["effect"]][["mut.var.alpha"]],
                                              mut_count = cfg[["effect"]][["mut.perc"]]$count,
                                              Gvar = Gvar, Gmut = Gmut,
                                              Npos = Npos, Nmut = Nmut,
                                              var_names = cfg[["effect"]][["var.dist"]]$var_group[-1],
                                              pos_flag = cfg[["effect"]][["pos.flag"]])

  # calculate true effect depending on experiments
  effects <- generate_true_effect(expected_effect = expected_result$expected_effect,
                                  type_exp = type_exp,
                                  wt_effect = cfg[["effect"]][["wt.effect"]],
                                  Nround = Nround, Nrep = Nrep, Npos = Npos, Nmut = Nmut)

  # add noise to true effect
  if(!is.null(cfg[["noise"]][["noise.freq"]])){
    effects <- add_effect_noise(effects = effects,
                                listn = expected_effect_list[['ctrl']],
                                noise_freq = cfg[["noise"]][["noise.freq"]],
                                type_exp = type_exp,
                                wt = cfg[["effect"]][["wt.effect"]],
                                Nround = Nround)
  }

  return(list(effects = effects,
              expected_effect = expected_result$expected_effect,
              expected_group = expected_result$expected_Gvar))
}

GenerateCount <- function(cfg, effects) {

  Nround <- cfg[["exp"]][["n.round"]]
  Nrep <- cfg[["exp"]][["n.rep"]]
  Npos <- cfg[["exp"]][["n.pos"]]
  Nmut <- cfg[["exp"]][["n.mut"]]
  Ncell_pop <- cfg[["pop"]][["pop.size"]]
  type_exp <- cfg[["sim"]][["type.sim"]]

  # starting count distribution
  dist_start <- starting_count_distribution(Ncell_pop = Ncell_pop, Npos = Npos, Nmut = Nmut,
                                            disp_start = cfg[["pop"]][["lib.disp"]] * cfg[["pop"]][["lib.shrink"]])

  # generate starting counts
  counts <- generate_starting_counts(dist_start = dist_start,
                                     Nround = Nround, Nrep = Nrep, Npos = Npos, Nmut = Nmut,
                                     replicate_mode = cfg[["exp"]][["mode.rep"]])

  # run experiment
  if(type_exp == 'binding') {
    # TODO: test binding screen
    stop("Function not tested yet.")
    counts <- run_binding(counts = counts,
                          effects = effects,
                          Ncell_pop = Ncell_pop, Nround = Nround, Nrep = Nrep,
                          messages = TRUE)
  } else if (type_exp == 'growth') {
    counts <- run_growth(counts = counts,
                         effects = effects,
                         Ncell_pop = Ncell_pop, Nround = Nround, Nrep = Nrep,
                         wt = cfg[["effect"]][["wt.effect"]],
                         messages = TRUE)
  } else{
    stop(paste('Invalid assay mode', type_exp))
  }

  # generate sequencing count
  sequencing <- generate_sequencing_counts(counts = counts,
                                           depth = cfg[["seq"]][["seq.depth"]],
                                           dispersion = cfg[["seq"]][["seq.disp"]] * cfg[["seq"]][["seq.shrink"]],
                                           Nrep = Nrep,
                                           Nround = Nround)

  return(list(counts = counts,
              sequencing = sequencing))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Generate row name for counts/effects
# Warning: Nmut + 1 mutants in total (one is ctrl)
# Example: ctrl_pos001, mutant01_pos001
mutant_index <- function(Npos, Nmut){
  mut_labels <- formatC(seq(1:Nmut), width = nchar(as.character(Nmut)), flag = "0")
  mut_labels <- c('ctrl', paste0('mutant', mut_labels))

  pos_labels <- formatC(seq(1:Npos), width = nchar(as.character(Npos)),
                        flag = "0")
  pos_labels <- paste0('pos', pos_labels, '_')
  index <- paste0(rep(pos_labels, each = Nmut + 1), rep(mut_labels, Nmut + 1))

  return(index)
}

# Generate column name for counts/effects
# Example: rep1.c_0
rep_index <- function(Nround, Nrep){
  index <- paste0(rep(paste0('rep', seq(1:Nrep)), each = Nround + 1),
                  rep(paste0('.c_', seq(0, Nround)), Nrep))
  return(index)
}


# A multinomial distribution generator for starting cell population
starting_count_distribution <- function(Ncell_pop, Npos, Nmut, disp_start){
  Nvar <-  Npos * (Nmut + 1)
  gen <- function() {
    p <- as.numeric(compositions::rDirichlet.rcomp(n = 1, alpha = rep(1/Nvar, Nvar) * disp_start * Nvar))
    distln <- stats::rmultinom(n = 1, size = Ncell_pop, prob = p)
    return(distln)
  }

  return(gen)
}

# Initialize cell count data.frame and generate starting cell population
generate_starting_counts <- function(dist_start,
                                     Nround, Nrep, Npos, Nmut,
                                     replicate_mode){
  counts_df <- data.frame(matrix(0, nrow = Npos * (Nmut + 1),
                                 ncol = (Nround + 1) * Nrep))
  colnames(counts_df) <- rep_index(Nround, Nrep)
  row.names(counts_df) <- mutant_index(Npos, Nmut)

  if (replicate_mode == 'bio'){
    # independent draws for each replicate
    counts_df[, grep("c_0$", colnames(counts_df))] <-
      replicate(Nrep, as.numeric(dist_start()), simplify = FALSE)
  } else if(replicate_mode == 'tech'){
    # same draw for all replicates
    counts_df[, grep("c_0$", colnames(counts_df))] <- as.numeric(dist_start())
  } else{
    stop(paste("Invalid replicate mode", replicate_mode))
  }

  return(counts_df)
}


# A normal distribution generator for expected effect
expected_effect_normal <- function(loc, scale) {
  gen <- function(size) {
    distn <- stats::rnorm(size, mean = loc, sd = scale)
    return(distn)
  }
  return(gen)
}
# zero generator for expected effect
expected_effect_zero <- function() {
  gen <- function(size) {
    return(rep(0, size))
  }
  return(gen)
}

# A list of distributions for ctrl and each variant group
expected_effect_distribution_list <- function(var_para,
                                              var_shrink_normal,
                                              effect_scale,
                                              null_group){
  syn_mean <- var_para[var_para$var_group == 'ctrl',]$mean

  expected_effect <- list()

  expected_effect <- lapply(1:nrow(var_para), function(i){
    row <- var_para[i, ]
    if (row$var_group == null_group || row$var_group == "ctrl") {
      expected_effect_zero()
    } else {
      expected_effect_normal(loc = (row$mean - syn_mean) * effect_scale,
                             scale = var_shrink_normal * row$sd * effect_scale^2)
    }
  })

  names(expected_effect) <- var_para$var_group
  return(expected_effect)
}

# Generate expected effect based on mutant distribution
generate_expected_effect <- function(expected_effect_list, mut_para, mut_count,
                                     Gvar, Npos, Gmut, Nmut, var_names, pos_flag){

  # from alpha to percentage
  # first compute percentage of (mut group i, var group j) at each position
  # then normalize percentage of var group within each mut group at each position
  # mut_count is the total number of variants in each mutant group
  alpha <- mut_para$alpha
  mx_perc <- matrix(as.vector(t(compositions::rDirichlet.rcomp(Npos, alpha))),
                    nrow = Gmut * Npos, ncol = Gvar - 1, byrow = TRUE)
  df_perc <- as.data.frame(t(apply(mx_perc, 1, function(x) x/sum(x))))
  df_perc$total <- rep(mut_count, Npos)

  # from percentage to count
  # generate number of variants in each var group from each mut group at each position
  mx_count <- t(apply(df_perc, 1,
                      function(row) rmultinom(1, row[Gvar], prob = row[1:(Gvar - 1)])))
  if (pos_flag) {
    perc <- colSums(mx_count) / sum(mx_count)
    pgroup_count <- stats::rmultinom(n = 1, size = Npos, prob = perc)
    pgroup_count <- rep(1:(Gvar-1), times = pgroup_count)

    for (i in 1:Npos) {
      mx_count[(i*4-3):(i*4), pgroup_count[i]] <- mut_count
      mx_count[(i*4-3):(i*4), -pgroup_count[i]] <- 0
    }
  }

  # from count to variant group list
  var_labels <- mut_para$var_group
  expected_value <- matrix(0, nrow = Npos, ncol = Nmut + 1)
  expected_Gvar <- matrix("", nrow = Npos, ncol = Nmut + 1)
  for(i in 0:(Npos-1)){
    keys = c('ctrl')
    for(j in 1:Gmut) {
      key <- rep(var_names, mx_count[i*Gmut+j,])
      key <- sample(key)
      keys <- c(keys, key)
    }
    expected_Gvar[i + 1,] <- keys
    for(k in 1:(Nmut + 1)) {
      expected_value[i + 1, k] <- expected_effect_list[[keys[k]]](size = 1)
    }
  }

  result <- list()
  result[['expected_effect']] <- c(apply(expected_value, 1, c)) # score
  result[['expected_Gvar']] <- c(apply(expected_Gvar, 1, c)) # var label

  return(result)
}

# Generate true effects (helper function)
# probability being selected/growth rate depending on simulation type
calc_true_effect_from_expected <- function(expected_effect, type_exp, wt, Nround){

  if(type_exp == 'binding'){
    true_effect <- exp(expected_effect/Nround + log(wt))
    if(sum(true_effect > 1) > 0){
      print(paste(sum(true_effect > 1), "variants have binding effect more than 1. Shrink to 1."))
      true_effect[true_effect > 1] <- 1
    }
  }
  else if(type_exp == 'growth'){
    if(wt < 0){
      wt_effect <- -1
      true_effect <- wt_effect * (1 + (expected_effect/log(2)/Nround/wt))
      if(sum(true_effect > 0) > 0){
        warning(sum(true_effect > 0), " variants have growth effect larger than 0. ",
                "Maximum true effect is ", max(true_effect),
                call. = FALSE, immediate. = TRUE)
        # true_effect[true_effect > 0] <- 0
      }
    }
    else if(wt > 0){
      wt_effect <- 1
      true_effect = wt_effect * (1 + (expected_effect/log(2)/Nround/wt))
      if(sum(true_effect < 0) > 0){
        warning(sum(true_effect < 0), "variants have growth effect smaller than 0.",
                "Minimum true effect is ", min(true_effect),
                call. = FALSE, immediate. = TRUE)
        # true_effect[true_effect < 0] <- 0
      }
    }
    else{
      stop(paste("Invalid wild type doubling rate", wt))
    }
  }
  else{
    stop(paste("Invalid assay mode ", type_exp))
  }

  return(true_effect)
}

# Genrerate true effects: dataframe
generate_true_effect <- function(expected_effect,
                                 type_exp,
                                 wt_effect,
                                 Nround, Npos, Nmut, Nrep){

  idx <- mutant_index(Npos, Nmut)
  true_effect <-
    calc_true_effect_from_expected(expected_effect = expected_effect,
                                   type_exp = type_exp,
                                   wt = wt_effect,
                                   Nround = Nround)

  df <- data.frame(matrix(true_effect, nrow = length(true_effect),
                          ncol = Nrep + 1, byrow = FALSE))
  colnames(df) <- c('effect', paste0('rep', seq(1:Nrep)))
  row.names(df) <- idx

  return(df)
}

# randomly add noise to effects
add_effect_noise <- function(effects, listn, noise_freq, type_exp, wt, Nround){
  if(is.null(noise_freq)){
    stop('not applicable to noise-free experiment')
  }

  affected_variants <- sample(row.names(effects),
                              size=round(nrow(effects)* noise_freq),
                              replace = FALSE)
  indices <- round(seq(1, length(affected_variants), length.out = ncol(effects)),
                   digits = 0)

  for(i in 1:(ncol(effects) - 1)){
    new_effect <- listn(size = indices[i + 1] - indices[i])
    new_effect <- calc_true_effect_from_expected(new_effect, type_exp, wt, Nround)
    effects[affected_variants[indices[i]:(indices[i + 1] - 1)], i + 1] <- new_effect
  }

  return(effects)
}

# Resample cells to initial cell population
resample_counts <- function(count, depth, dispersion=0){

  if (dispersion > 0){
    alpha <- nrow(count) * dispersion * count / sum(count)
    var_zero <- subset(alpha, alpha[, 1] == 0)
    var_nonzero <- subset(alpha, alpha[, 1] != 0)
    p <- c(compositions::rDirichlet.rcomp(n = 1, alpha = var_nonzero[,1]))
    p <- data.frame(p, row.names = row.names(var_nonzero))
    colnames(p) <- colnames(var_zero)
    p <- rbind(p, var_zero)
    p <- p[rownames(count), 1]
  } else if(dispersion == 0){
    p <- count / sum(count)
    p <- p[, 1]
  } else{
    stop(paste("Invalid dispersion value", dispersion))
  }

  sample_counts <- stats::rmultinom(n = 1, size = length(p) * as.integer(depth), p)
  return(unname(sample_counts))
}

# TODO: NOT TESTED!!!
# Run binding screen
run_binding <- function(counts, effects, Ncell_pop, Nround, Nrep, messages = TRUE){
  for(rep in 1:Nrep){
    for(round in 1:Nround){
      for(row in 1:nrow(counts)){
        Nprev <- counts[row, paste0('rep', rep, '.c_', round - 1)]
        p <- effects[row, paste0('rep', rep)]
        counts[row, paste0('rep', rep, '.c_', round)] <- stats::rbinom(1, Nprev, p)
      }

      # grow cells to original population
      counts[paste0('rep', rep, '.c_', round)] <-
        resample_counts(counts[paste0('rep', rep, '.c_', round)], depth = Ncell_pop / nrow(counts))

      if(messages){
        print(paste('Finished rep', rep, "round", round))
      }
    }
  }

  return(counts)
}

# Run growth screen
run_growth <- function(counts, effects, Ncell_pop, Nround, Nrep, wt, messages = TRUE){

  wt_effect <- sign(wt)
  if (wt == 0) {
    stop(paste("Invalid wild type doubling rate", wt))
  }
  t <- wt * log(2) / wt_effect

  for(rep in 1:Nrep){
    for(round in 1:Nround){
      for(row in 1:nrow(counts)){
        c <- counts[row, paste0('rep', rep, '.c_', round - 1)]
        if(c > 0){
          counts[row, paste0('rep', rep, '.c_', round)] <-
            stats::rpois(1, c * exp(t * effects[row, paste0('rep', rep)]))
          # if(wt > 0){
          #   counts[row, paste0('rep',rep,'.c_',round)] <-
          #     np.random.negative_binomial(c, exp(-1 * t * effects[row, paste0('rep',rep)])) + c
          # }
          # if(wt < 0){
          #   counts[row, paste0('rep',rep,'.c_',round)] <-
          #     rpois(1, c * exp(t * effects[row, paste0('rep',rep)]))
          # }
        } else{
          counts[row, paste0('rep', rep, '.c_', round)] <- 0
        }
      }

      # if exceeding initial cell population
      if(sum(counts[paste0('rep', rep, '.c_', round)]) > Ncell_pop){
        counts[paste0('rep', rep, '.c_', round)] <-
          as.numeric(resample_counts(counts[paste0('rep', rep, '.c_', round)],
                          depth = Ncell_pop / nrow(counts)))
      } else if(sum(counts[paste0('rep', rep, '.c_', round)] == 0) > 0.8 * nrow(counts)){
        stop(paste('Too many rounds of cell growth. 80% variants die at round', round))
      }

      if(messages){
        print(paste('Finished rep', rep, "round", round))
      }
    }
  }

  return(counts)
}

# Generate sequencing count with dispersion
generate_sequencing_counts <- function(counts, depth, dispersion, Nrep, Nround){

  seq_counts <- data.frame(matrix(0, nrow = nrow(counts), ncol = ncol(counts)))
  colnames(seq_counts) <- colnames(counts)
  row.names(seq_counts) <- row.names(counts)

  for(rep in 1:Nrep){
    for(round in 0:Nround){
      seq_counts[paste0('rep',rep,'.c_',round)] =
        as.numeric(resample_counts(counts[paste0('rep', rep, '.c_', round)], depth, dispersion))
    }
  }

  return(seq_counts)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal: output
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_tsv <- function(effects_list, counts_list, save.dir) {

  sub.save.dir <- file.path(save.dir, "tsv")
  if (!dir.exists(sub.save.dir)) {
    dir.create(sub.save.dir, recursive = TRUE)
  }

  # output effects
  df.effects.output <- data.frame(expected_effects = effects_list$expected_effect,
                                  expected_group = effects_list$expected_group,
                                  effects_list$effects)
  df.effects.output$variants <- rownames(df.effects.output)
  rownames(df.effects.output) <- NULL
  df.effects.output <- df.effects.output %>% dplyr::relocate(.data$variants)
  readr::write_tsv(df.effects.output, file = file.path(sub.save.dir, "effects.tsv"))

  # output counts
  counts <- counts_list$counts
  sequencing <- counts_list$sequencing
  readr::write_tsv(counts, file = file.path(sub.save.dir, "cell_counts.tsv"))
  readr::write_tsv(sequencing, file = file.path(sub.save.dir, "sequencing_counts.tsv"))

}

output_rosace <- function(effects_list, counts_list, Nrep, Nround, save.dir) {

  # saving directory
  sub.save.dir <- file.path(save.dir, "rosace")
  if (!dir.exists(sub.save.dir)) {
    dir.create(sub.save.dir, recursive = TRUE)
  }

  # create rosace object
  for (i in 1:Nrep) {
    idx_start <- (Nround + 1) * i - Nround
    idx_end <- (Nround + 1) * i
    assay <- CreateAssayObject(counts = as.matrix(counts_list$sequencing[idx_start:idx_end]),
                               var.names = rownames(counts_list$sequencing),
                               key = "simulation", rep = i, type = "growth")
    if (i == 1) {
      rosace <- CreateRosaceObject(object = assay)
    } else {
      rosace <- AddAssayData(object = rosace, assay = assay)
    }
  }

  # add ground truth data
  effects <- data.frame(variants = rownames(effects_list$effects),
                        expected_effects = effects_list$expected_effect,
                        expected_group = effects_list$expected_group)
  rosace@misc <- list(effects = effects)

  # process var.data
  rosace@var.data <- rosace@var.data %>%
    tidyr::separate(.data$variants, into = c("position", "mutation"), remove = FALSE) %>%
    mutate(position = as.numeric(substr(.data$position, 4, nchar(.data$position))))

  # save rosace object
  saveRDS(rosace, file = file.path(sub.save.dir, "rosace.rds"))

}

output_enrich2 <- function(counts_list, Nrep, Nround, mode.sim, save.dir) {

  # saving directory
  sub.save.dir <- file.path(save.dir, "enrich2")
  if (!dir.exists(sub.save.dir)) {
    dir.create(sub.save.dir, recursive = TRUE)
  }
  data.dir <- file.path(sub.save.dir, "data")
  if (!dir.exists(data.dir)) {
    dir.create(data.dir, recursive = TRUE)
  }
  result.dir <- file.path(sub.save.dir, "results")
  if (!dir.exists(result.dir)) {
    dir.create(result.dir, recursive = TRUE)
  }

  # add wild-type info
  ctrl_label <- endsWith(rownames(counts_list$sequencing), 'ctrl')

  # output counts
  for (i in 1:Nrep) {
    for (j in 0:Nround) {
      count <- counts_list$sequencing[(i - 1) * (Nround + 1) + j + 1]
      colnames(count) <- "count"
      count[nrow(count) + 1, 1] <- sum(count[[1]][ctrl_label])
      rownames(count)[nrow(count)] <- "_wt"
      utils::write.table(count, quote = FALSE, sep="\t",
                         file = file.path(data.dir, paste("count_rep", i, "_c", j, ".tsv", sep = "")))
    }
  }

  # output json file
  experiment <- list('name' = 'simulation', 'output directory' = result.dir, 'conditions' = list())
  condition <- list('name' = mode.sim, selections = list())
  for (i in 1:Nrep) {
    selection <- list('name' = paste("rep", i, sep = ""), 'libraries' = list())
    for (j in 0:Nround) {
      seqlib = list('counts file' = file.path(data.dir, paste("count_rep", i, "_c", j, ".tsv", sep = "")),
                    'identifiers' = list(),
                    'name' = paste("rep", i, "_c", j, sep = ""),
                    'report filtered read' = FALSE,
                    'timepoint' = j)
      selection$libraries[[j + 1]] <- seqlib
    }
    condition$selections[[i]] <- selection
  }
  experiment$conditions[[1]] <- condition

  # save json file
  jsonData <- rjson::toJSON(experiment, indent = 2)
  write(jsonData, file = file.path(sub.save.dir, "config.json"))
}


