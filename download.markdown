---
layout: page
title: Download
permalink: /download/
order: 2
---

To install `rosace` start R and first install `cmdstanr` by typing:

```
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# use cmdstanr to install CmdStan, this requires a working C++ toolchain and compiler
library(cmdstanr)
install_cmdstan(cores = 4)
```

Then install devtools by typing

```
install.packages("devtools")
```

and install `rosace` by typing

```
devtools::install_github("pimentellab/rosace")
```

Next load `rosace` with

```
library("rosace")
```

