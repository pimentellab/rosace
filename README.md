# rosace

__rosace__ is an R package for analyzing growth-based deep mutational scanning screen data.  

# Installation

__rosace__ uses [cmdstanr](https://mc-stan.org/cmdstanr/) to run inference. Please ensure that __cmdstanr__ is properly installed before installing __rosace__. Below is a concise installation command; for complete details, please refer to the official website and also check out potential issues reported in the issues_with_installation.Rmd under vignettes (will be mounted to a webpage soon). 
```{r eval=FALSE}
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# use cmdstanr to install CmdStan, this requires a working C++ toolchain and compiler
library(cmdstanr)
install_cmdstan(cores = 2)
```

To install __rosace__ start [R](https://www.r-project.org) and first install devtools by typing
```{r eval=FALSE}
install.packages("devtools")
```

and install __rosace__ by typing
```{r eval=FALSE}
devtools::install_github("pimentellab/rosace")
```

If downloading from local type
```{r eval=FALSE}
devtools::install(".") # change to the directory of rosace
```

Next load __rosace__ with
```{r}
library("rosace")
```

# Documentation

We recommend starting with the vignette...

# Further help

