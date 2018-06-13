**Important note for Windows users**:  If you are running Wasabi under R or R studio on Windows, there is a known issue (bug, I would say) that prevents normal operation.  Salmon and Sailfish output important extra information like bootstrap samples and the fragment length distribution to a subdirectory of the quantification directory named `aux`.  Windows, unfortunately, [forbids `aux` as a folder name](https://blog.onetechnical.com/2006/11/16/forbidden-file-and-folder-names-on-windows/).  If you are lucky enough to be reading this *before* running Sailfish/Salmon, you can pass an alternative folder name to the `--auxDir` command line option to avoid this issue.  Otherwise, we suggest the following work-around.  First, rename the `aux` subdirectory of each Sailfish/Salmon quantification directory to be a folder name that is not forbidden under Windows e.g. `aux2`.  Then, open the `cmd_info.json` file in the quantification directory, and add the following line:

```
"auxDir" : "aux2",
```

Now, things should work normally under windows.  Salmon versions (>=0.7.0) use `aux_info` as the default folder to avoid this issue.  Future releases of Sailfish (>0.10.1) will use a different default name for the auxiliary directory to avoid this issue as well.

# What is wasabi?

[![Join the chat at https://gitter.im/COMBINE-lab/wasabi](https://badges.gitter.im/COMBINE-lab/wasabi.svg)](https://gitter.im/COMBINE-lab/wasabi?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Wasabi allows you to easily prepare [Sailfish](https://github.com/kingsfordgroup/sailfish) and [Salmon](https://github.com/COMBINE-lab/salmon) output for downstream analysis.  
Currently, its main purpose it to prepare output for downstream analysis with [sleuth](http://pachterlab.github.io/sleuth/).

# How to use wasabi


## Installation 

First, you need to install the wasabi package.  There are two main ways to accomplish this:

#### Installation with devtools
  With `devtools`, it's easy to install `wasabi`:
  ```r
  source("http://bioconductor.org/biocLite.R")
  biocLite("devtools")    # only if devtools not yet installed
  biocLite("COMBINE-lab/wasabi")
  ```
    
#### Installation with bioconda
  Alternatively, you can use the [conda](http://conda.pydata.org/miniconda.html) package manager, along with the [bioconda](https://bioconda.github.io/) channel to install `wasabi`:
  ```
  conda create -n wasabi r-wasabi
  ```

## Loading and using wasabi

Once wasabi is installed, you can load the library with:

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
