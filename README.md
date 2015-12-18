# What is wasabi?

Wasabi allows you to easily prepare [Sailfish](https://github.com/kingsfordgroup/sailfish) and [Salmon](https://github.com/COMBINE-lab/salmon) output for downstream analysis.  
Currently, its main purpose it to prepare output for downstream analysis with [sleuth](http://pachterlab.github.io/sleuth/).

# How to use wasabi

First, install the wasabi package.  This is easy with `devtools`:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("devtools")    # only if devtools not yet installed
biocLite("COMBINE-lab/wasabi")
```

You can load the library with:

```r
library(wasabi)
```

Now that wasabi is installed, we can use it to convert the Sailfish / Salmon output into sleuth-compatible format.
Imagine we have some samples (Sailfish quant directories) sitting in a data directory:

```
data/samp1   data/samp2   data/samp3   data/samp4
```

First, we create a simple vector containing these directories (the `>` below is an R prompt):

```r
> sfdirs <- file.path("data", c("samp1", "samp2", "samp3"))
```

Now, we simply run the `prepare_fish_for_sleuth` function:

```r
> prepare_fish_for_sleuth(sfidrs)
```

The function will write some status messages to the console and, when it's done, each directory will now contain 
and `abundance.h5` file in a sleuth-compatible format.  From this point forward, you can simply run sleuth normally.
