#' @importFrom readr read_tsv
#' @importFrom jsonlite fromJSON
#' @importFrom rjson toJSON
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Write the config list in json format
#'
#' @param cfg A list of config
#' @param file File name
#' @param savedir Directory to save the file
#'
#' @return None
#'
#' @export
#'
WriteConfig <- function(cfg, file, savedir) {
  jsonData <- rjson::toJSON(cfg, indent = 2)
  write(jsonData, file = paste0(file.path(savedir, file),".json"))
}

#' Read a full config file
#' (with null.var.group, mut.var.alpha, var.dist as data frame)
#'
#' @param file Config file path
#'
#' @return A list of config
#'
#' @rdname ReadConfig
#' @export
#'
ReadFullConfig <- function(file) {
  cfg <- jsonlite::fromJSON(file)
  return(cfg)
}

#' Read a path config file
#' (with null.var.group, mut.var.alpha, var.dist as file path)
#'
#' @param file Config file path
#'
#' @return A list of config
#'
#' @rdname ReadConfig
#' @export
#'
ReadPathConfig <- function(file) {
  cfg <- jsonlite::fromJSON(file)
  effect <- cfg$effect
  cfg$effect <-
    listEffectConfig(null.var.group = effect$null.var.group,
                     mut.var.alpha = effect$mut.var.alpha, # need to read
                     var.dist = effect$var.dist, # need to read
                     mut.perc = effect$mut.perc, # need to read
                     rounds = effect$rounds,
                     wt.effect = effect$wt.effects,
                     var.shrink = effect$var.shrink)
  return(cfg)
}


#' Create the config list for simulation from a Rosette object
#'
#' @param object A Rosette object
#' @param n.sim Number of simulations
#' @param save.sim Directory to save the simulation data
#' @param type.sim Type of simulation: growth or binding
#' @param n.rep Number of replicates
#' @param null.var.group Name of the null variant group
#' @param wt.effect Wild-type effect (binding) or doubling rate (growth)
#' @param n.round Number of rounds
#' @param n.mut Number of mutants (currently has to be the same as in the Rosette object)
#' optional, can be inferred from the Rosette object.
#' @param n.pos Number of positions
#' optional, can be inferred from the Rosette object.
#' @param pop.size Cell population size
#' optional, by default 1000 cells per variant.
#' @param mode.sim Simulation mode: "clean" (default) or "reperror"
#' @param mode.rep Replicate mode: "bio" (default) or "tech"
#' @param seq.shrink Shrinkage factor for sequencing dispersion (times seq.disp)
#' @param seq.depth Sequencing depth, by default 200
#' @param lib.shrink Shrinkage factor for library dispersion
#' @param var.shrink Shrinkage factor for variance (var.dist)
#' @param pos.flag whether the simulation is position-centric
#' @param ... Other parameters for noise.
#' Reperror: noise.freq, noise.mult
#'
#' @return A list of config
#'
#' @export
#'
# (1): input
# (2): optional. missing okay. can be inferred from rosette.
# (3): optional. have default value.
# (4): no input. directy from rosette.
CreateConfig <-
  function(object,
           n.sim, save.sim, type.sim, # sim (1)
           n.rep, # exp (1)
           null.var.group, wt.effect, # effect (1)
           n.round, n.mut, n.pos, # exp: optional from rosette (2)
           pop.size, # pop: optional from rosette (2)
           mode.sim = "clean", # sim (3)
           mode.rep = "bio", # exp (3)
           seq.shrink = 2, seq.depth = 200,  # seq (3)
           lib.shrink = 2, # pop (3)
           var.shrink = 0.8, pos.flag = TRUE, # effect (3)
           ... # noise: noise.freq, noise.mult (1)
           # seq.disp # seq: from rosette (4)
           # lib.disp # pop: from rosette (4)
           # mut.var.alpha, var.dist, mut.perc, rounds # effect: from rosette (4)
           ) {


    # sim config
    sim <- listSimConfig(n.sim = n.sim, save.sim = save.sim, type.sim = type.sim, # sim (1)
                         mode.sim = mode.sim) # sim (3)

    # exp config
    if (missing(n.round)) {
      n.round <- object@rounds
    }
    if (missing(n.mut)) {
      n.mut <- nrow(object@mut.label)
    }
    if (missing(n.pos)) {
      n.pos <- length(unique(object@score.df$pos))
    }
    exp <- listExpConfig(n.rep = n.rep, # exp (1)
                         mode.rep = mode.rep, # exp (3)
                         n.round = n.round, n.mut = n.mut, n.pos = n.pos) # exp: optional from rosette (2)

    # pop config
    if (missing(pop.size)) {
      pop.size <- (n.mut + 1) * n.pos * 200
    }
    lib.disp <- object@disp.start
    pop <- listPopConfig(pop.size = pop.size, # pop: optional from rosette (2)
                         lib.disp = lib.disp, # pop (4)
                         lib.shrink = lib.shrink) # pop (3)

    # seq config
    if (length(object@disp) == 0) {
      stop("Object has empty slot 'disp'. Run 'AddDisp' first.")
    }
    seq.disp <- object@disp
    seq <- listSeqConfig(seq.disp = seq.disp,  # seq: from rosette (4)
                         seq.shrink = seq.shrink, seq.depth = seq.depth) # seq (3)

    # effect config
    output <- GenerateOutput(object, save.files = FALSE)
    if (length(object@rounds) == 0) {
      stop("Object initialization error. Slot 'rounds' is empty.")
    }

    mut.var.alpha <- output$mut.var.alpha
    var.dist <- output$var.dist
    mut.perc <- output$mut.perc
    rounds <- object@rounds
    effect <-
      listEffectConfig(null.var.group = null.var.group, wt.effect = wt.effect, # effect (1)
                       mut.var.alpha = mut.var.alpha, var.dist = var.dist,
                       mut.perc = mut.perc, rounds = rounds, # effect: from rosette (4)
                       var.shrink = var.shrink, pos.flag = pos.flag) # effect (3)

    # noise config
    noise <- listNoiseConfig(mode.sim = mode.sim, ...)

    cfg <- listConfig(list.sim = sim,
                      list.exp = exp,
                      list.seq = seq,
                      list.pop = pop,
                      list.effect = effect,
                      list.noise = noise)

    return(cfg)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

listConfig <- function(list.sim, list.exp, list.seq,
                       list.pop, list.effect, list.noise) {
  return(list(sim = list.sim,
              exp = list.exp,
              seq = list.seq,
              pop = list.pop,
              effect = list.effect,
              noise = list.noise))
}

listSimConfig <- function(n.sim, save.sim, type.sim,
                          mode.sim = "clean") {

  return(list(n.sim = n.sim, # number of simulation
              save.sim = save.sim, # saving directory
              type.sim = type.sim, # type of experiment: binding or growth
              mode.sim = mode.sim))  # mode: clean, (reperror or jackpot)

}

listExpConfig <- function(n.rep,
                          mode.rep = "bio",
                          n.round = 3,
                          n.mut = 20, n.pos = 500) {

  return(list(n.rep = n.rep, # number of replicates
              mode.rep = mode.rep, # replicate mode
              n.round = n.round, # number of rounds
              n.mut = n.mut, # number of mutants
              n.pos = n.pos)) # number of positions
}

listSeqConfig <- function(seq.disp, seq.shrink = 1.5, seq.depth = 200) {

  if (seq.shrink < 1) {
    stop("Seq config: sequencing dispersion shrinkage has to be no less than 1.")
  }

  return(list(seq.disp = seq.disp, # sequencing dispersion ***import from rosette
              seq.shrink = seq.shrink, # shrink dispersion
              seq.depth = seq.depth)) # sequencing depth
}

listPopConfig <- function(pop.size, lib.disp, lib.shrink) {

  return(list(pop.size = pop.size,
              lib.disp = lib.disp,
              lib.shrink = lib.shrink)) # population size

}

listEffectConfig <- function(null.var.group,
                             mut.var.alpha, var.dist, mut.perc, rounds,
                             wt.effect, var.shrink = 0.8, pos.flag = TRUE) {

  if (!"data.frame" %in% class(mut.var.alpha)) {
    mut.var.alpha <- readr::read_tsv(mut.var.alpha)
  }
  if (!'data.frame' %in% class(var.dist)) {
    var.dist <- readr::read_tsv(var.dist)
  }
  if (!'data.frame' %in% class(mut.perc)) {
    mut.perc <- readr::read_tsv(mut.perc)
  }

  if (!null.var.group %in% var.dist$var_group) {
    stop("Effect config: invalid 'null.var.group' entry.
         Check var group distribution and choose.")
  }

  if (var.shrink > 1) {
    stop("Effect config: var.shrink has to be no more than 1.")
  } else if (var.shrink < 0.5) {
    warnings("Effect config: var.shink is lower than 0.5.")
  }

  return(list(null.var.group = null.var.group, # null variant group
              mut.var.alpha = mut.var.alpha, # DATAFRAME ***import from rosette
              var.dist = var.dist, # DATAFRAME ***import from rosette
              mut.perc = mut.perc, # DATAFRAME ***import from rosette
              rounds = rounds, # max rounds of exp ***import from rosette
              wt.effect = wt.effect, # wt effect
              var.shrink = var.shrink, # var.shrink
              pos.flag = pos.flag))
}

listNoiseConfig <- function(mode.sim, noise.freq = 0.1, noise.mult = 50) {

  if (mode.sim == "clean") {
    return(list())
  } else if (mode.sim == "reperror" && !missing(noise.freq)) {

    if (noise.freq > 1 || noise.freq <= 0) {
      stop("Noise config: noise frequence has to be (0, 1].")
    }

    return(list(noise.freq = noise.freq))
  } else if (mode.sim == "jackpot" && !missing(noise.freq) && !missing(noise.mult)) {

    if (noise.freq > 1 || noise.freq <= 0) {
      stop("Noise config: noise frequence has to be (0, 1].")
    }
    if (noise.mult <= 1) {
      stop("Noise config: jackpot noise multiplication has to be > 1.")
    }

    return(list(noise.freq = noise.freq,
                noise.mult = noise.mult))
  } else {
    stop("Noise config input error.
         If mode is reperror input 'noise.freq'.
         If mode is jackpot input both 'noise.freq' and 'noise.mult'.")
  }

}




