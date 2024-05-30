# rosace

<p align="left">
  <img src="man/figures/rosace_logo.png" width="150">
</p>

## Update

v1.1: May 30, 2024
- Added options to take out nonsense/stop mutations from position-level estimates
- Fixed a bug in position-level lfsr computation 

## Overview

__rosace__ is an R package for analyzing growth-based deep mutational scanning screen data.  

## Installation

__rosace__ uses [cmdstanr](https://mc-stan.org/cmdstanr/) to run inference. Please ensure that __cmdstanr__ is properly installed before installing __rosace__. Below is a concise installation command; for complete details, please refer to the official website. 
```{r eval=FALSE}
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# use cmdstanr to install CmdStan, this requires a working C++ toolchain and compiler
library(cmdstanr)
install_cmdstan(cores = 4)
```

To install __rosace__ start [R](https://www.r-project.org) and first install devtools by typing
```{r eval=FALSE}
install.packages("devtools")
```

and install __rosace__ by typing
```{r eval=FALSE}
devtools::install_github("pimentellab/rosace")
```

If you prefer to use Docker, we also provide a Docker [image](https://hub.docker.com/r/cbmacdo/rosace-docker) for rosace. You can pull the image in the command line with
```sh
docker pull cbmacdo/rosace-docker
```

See the full [Installation Instructions](vignettes/installation_instructions.Rmd) for further details and alternative installation options.

## Getting started

```{r eval=FALSE}
library("rosace")
```

We recommend starting with the [vignette](https://pimentellab.github.io/rosace/doc/intro_rosace.html). A vignette for the simulation module [Rosette](vignettes/rosette_simulation.Rmd) is also avaliable.

## Further help

You may submit a bug report here on GitHub as an issue or you could send an email to roserao@ucla.edu.

## Citing rosace

Please cite the following publication if you use __rosace__: Rao, J., Xin, R., Macdonald, C. et al. Rosace: a robust deep mutational scanning analysis framework employing position and mean-variance shrinkage. Genome Biol 25, 138 (2024). https://doi.org/10.1186/s13059-024-03279-7
