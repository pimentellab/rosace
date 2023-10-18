# rosace

<p align="left">
  <img src="man/figures/rosace_logo.png" width="150">
</p>

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

See the full [Installation Instructions](vignettes/issues_with_installation.Rmd) for further details and alternative installation options.

## Getting started

```{r}
library("rosace")
```

We recommend starting with the [vignette](vignettes/intro_rosace.Rmd). 

## Further help

You may submit a bug report here on GitHub as an issue or you could send an email to roserao@ucla.edu.

## Citing rosace

Please cite the following publication if you use __rosace__: xxxxxx