#' Convert sailfish/salmon results for one or more samples to kallisto HDF5
#'
#' @param fish_dirs a character vector of length greater than one where each
#' string points to a sailfish/salmon output directory
#' @export
prepare_fish_for_sleuth <- function(fish_dirs, force=FALSE, fallback_mu=200, fallback_sd=80, fallback_num_reads=-1) {
  
  len_mu = length(fallback_mu)
  len_sd = length(fallback_sd)
  len_nr = length(fallback_num_reads)
  len_dir = length(fish_dirs)
  
  # from : http://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
  compare <- function(v) all(sapply( as.list(v[-1]), FUN=function(z) {identical(z, v[1])}))
  same_len = compare(c(1, len_mu, len_sd, len_nr))
  
  # We'll assume, for the time being, that all paths are from 
  # the same version of the software
  testdir <- fish_dirs[1]
  
  # The default aux directory
  aux_dir = "aux"
  if (file.exists(file.path(testdir, "cmd_info.json"))) {
    cmd_info <- rjson::fromJSON(file=file.path(testdir, "cmd_info.json"))
    # If the user provided a different auxiliary directory, use that instead
    if (is.element("auxDir", names(cmd_info))) {
      aux_dir = cmd_info$auxDir
    }
  }
  # The path to the aux directory if it exists
  auxPath <- file.path(testdir, aux_dir)

  ## If we're dealing with the new format files 
  if (file.exists(file.path(auxPath, "meta_info.json"))) {
    if (!same_len) {
      same_len = compare(c(len_dir, len_nr))
    }
    if (same_len) {
      mapply(fish_to_hdf5, fish_dirs, force=force, fallback_num_reads=fallback_num_reads)
    }
  } else {
  ## We're dealing with the old format files
    same_len = compare(c(1, len_mu, len_sd, len_nr))
    if (!same_len) {
      same_len = compare(c(len_dir, len_mu, len_sd, len_nr))
    }
    if (same_len) {
      mapply(fish_to_hdf5_old, fish_dirs, force=force, 
           fallback_mu=fallback_mu, fallback_sd=fallback_sd, fallback_num_reads=fallback_num_reads)
    }
  }
  
  if (!same_len) {
    msg <- paste("ERROR: If you are passing in fallback arguments (mu, sd, and num_reads) there must",
                 "be either a single value for each to apply to all results, or it must be a",
                 "vector of the same length as fish_dirs");
    message(msg)
    FALSE
  } else {
    TRUE
  }
}

#' Convert sailfish results in new format
#' SF ver >= 0.9.0, Salmon ver >= 0.6.0 for one sample to
#' kallisto HDF5
#'
#' @param fish_dir path to a sailfish output directory
#' @param force if TRUE re-create the h5 file even if it exists
fish_to_hdf5 <- function(fish_dir, force, fallback_num_reads) {
  h5file <- file.path(fish_dir, 'abundance.h5')
  if (!force && file.exists(h5file)) {
    print(paste("Skipping conversion: abundance.h5 already in ", fish_dir))
    return()
  }
  # If we're forcing it, then we have to remove the file now
  # or h5createFile will complain later
  if (file.exists(h5file)) {
    file.remove(h5file)
  }

  # load quantification data
  quant <- fread(file.path(fish_dir, 'quant.sf'))
  setnames(quant, c('target_id', 'length', 'eff_length', 'tpm', 'est_counts'))
  if (!is.character(quant$target_id)) {
    quant$target_id <- as.character(quant$target_id)
  }

  aux_dir = "aux"
  cmd_info <- rjson::fromJSON(file=file.path(fish_dir, "cmd_info.json"))
  # If the user provided a different auxiliary directory, use that instead
  if (is.element("auxDir", names(cmd_info))) {
    aux_dir = cmd_info$auxDir
  }
  auxPath <- file.path(fish_dir, aux_dir)
  
  # get all of the meta info
  minfo <- rjson::fromJSON(file=file.path(auxPath, "meta_info.json"))
  sampType <- NULL
  # check if we have explicitly recorded the type of posterior sample
  # (salmon >= 0.7.3)
  if (is.element("samp_type", names(minfo))) {
    sampType <- minfo$samp_type
  }
  
  
  # load bootstrap data if it exists
  knownSampleTypes <- c("gibbs", "bootstrap")
  numBoot <- minfo$num_bootstraps
  if (numBoot > 0) {
    bootCon <- gzcon(file(file.path(auxPath, 'bootstrap', 'bootstraps.gz'), "rb"))
    #knownSampleType <- ifelse(is.null(sampType), FALSE, sampType %in% knownSampleTypes)
    #if (knownSampleType) {
    #  boots <- readBin(bootCon, "double", n = minfo$num_targets * minfo$num_bootstraps)
    #} else {
    #}
    
    ##
    # Gibbs samples *used* to be integers, and bootstraps were doubles
    # Now, however, both types of samples are doubles.  The code below 
    # tries to load doubles first, but falls back to integers if it fails.
    ##   
    expected.n <- minfo$num_targets * minfo$num_bootstraps
    boots <- tryCatch({
      boots.in <- readBin(bootCon, "double", n = expected.n)
      stopifnot(length(boots.in) == expected.n)
      boots.in
    }, error=function(...) {
      # close and re-open the connection to reset the file
      close(bootCon)
      bootCon <- gzcon(file(file.path(auxPath, 'bootstrap', 'bootstraps.gz'), "rb"))
      readBin(bootCon, "integer", n = expected.n)
    })
    close(bootCon)
    boots <- as.double(boots)
    # rows are transcripts, columns are bootstraps
    dim(boots) <- c(minfo$num_targets, minfo$num_bootstraps)
  }

  # load stats
  numProcessed <- minfo$num_processed

  # If we're dealing with salmon, this field should be defined
  if (is.element("mapping_type", names(minfo))) {
    # In alignment mode (since unaligned reads may not be included in the BAM file)
    # we prefer to have the user supply the total number of reads in the input library
    if (minfo$mapping_type == "alignment") {
      if (fallback_num_reads == -1) {
        msg <- paste("Since salmon was run in alignment mode, it is recommended you provide",
                     "the total number of input reads via the fallback_num_reads argument to",
                     "prepare_fish_for_sleuth()")
        message(msg)
      } else {
        msg <- paste("Fish dir =", fish_dir, ":: fallback # reads =", fallback_num_reads)
        message(msg)
        numProcessed <- fallback_num_reads
      }
    } else if(minfo$mapping_type == "mapping" && fallback_num_reads != -1) {
      msg <- paste("It is generally NOT recommended to provide a fallback # of reads",
                   "when salmon was run in mapping-based mode, as it should already",
                   "accurately record this information. Please make sure you really",
                   "mean to do this!")
      message(msg)
      msg <- paste("Fish dir =", fish_dir, ":: fallback # reads =", fallback_num_reads)
      message(msg)
      numProcessed <- fallback_num_reads
    }
  }
  # Otherwise, this is either sailfish or salmon's mapping-mode and the 
  # numProcessed estimate should be OK.
  
  # build the hdf5
  rhdf5::h5createFile(h5file)

  # counts are at root
  rhdf5::h5write(quant$est_counts, h5file, '/est_counts')

  # aux group has metadata about the run and targets
  rhdf5::h5createGroup(h5file, 'aux')
  rhdf5::h5write(numProcessed, h5file, 'aux/num_processed')
  rhdf5::h5write(numBoot, h5file, 'aux/num_bootstrap')
  rhdf5::h5write(quant$length, h5file, 'aux/lengths')
  rhdf5::h5write(quant$eff_length, h5file, 'aux/eff_lengths')
  rhdf5::h5write(quant$target_id, h5file, 'aux/ids')
  rhdf5::h5write('10', h5file, 'aux/index_version')
  rhdf5::h5write('sailfish', h5file, 'aux/kallisto_version')
  rhdf5::h5write(timestamp(prefix="", suffix=""), h5file, "aux/start_time")

  # bootstrap group has (.. wait for it ..) bootstrap data
  if (numBoot > 0) {
    rhdf5::h5createGroup(h5file, 'bootstrap')
    sapply(0:(numBoot-1), function(i) {
      bootid <- paste('bs', i, sep='')
      rhdf5::h5write(unlist(boots[,i+1]),
              h5file, paste('bootstrap', bootid, sep='/'))
    })
  }

  bootCon <- gzcon(file(file.path(auxPath, 'fld.gz'), "rb"))
  fld <- readBin(bootCon, "int", n=minfo$frag_dist_length)
  close(bootCon)
  rhdf5::h5write(fld, h5file, 'aux/fld')

  bObsCon <- gzcon(file(file.path(auxPath, 'observed_bias.gz'), "rb"))
  bObs <- readBin(bObsCon, "int", n=minfo$num_bias_bins)
  close(bObsCon)
  rhdf5::h5write(bObs, h5file, 'aux/bias_observed')

  bExpCon <- gzcon(file(file.path(auxPath, 'expected_bias.gz'), "rb"))
  bExp <- readBin(bObsCon, "double", n=minfo$num_bias_bins)
  close(bExpCon)
  rhdf5::h5write(bExp, h5file, 'aux/bias_normalized')


  if(utils::packageVersion('rhdf5') < "2.23.0") {
    rhdf5::H5close()
  } else {
    rhdf5::h5closeAll()
  }
  print(paste("Successfully converted sailfish / salmon results in", fish_dir, "to kallisto HDF5 format"))
}

#' Convert sailfish results in the older format
#' sf ver <= 0.8.0, salmon ver <= 0.5.1 for one sample
#' to kallisto HDF5
#'
#' @param fish_dir path to a sailfish output directory
#' @param force if TRUE re-create the h5 file even if it exists
fish_to_hdf5_old <- function(fish_dir, force, fallback_mu, fallback_sd, fallback_num_reads) {
  h5file <- file.path(fish_dir, 'abundance.h5')
  if (!force && file.exists(h5file)) {
    print(paste("Skipping conversion: abundance.h5 already in ", fish_dir))
    return()
  }
  # If we're forcing it, then we have to remove the file now
  # or h5createFile will complain later
  if (file.exists(h5file)) {
    file.remove(h5file)
  }

  # load quantification data
  quant <- fread(file.path(fish_dir, 'quant.sf'))
  setnames(quant, c('target_id', 'length', 'tpm', 'est_counts'))
  # NOTE: this can re-order the rows!
  if (!is.character(quant$target_id)) {
    quant$target_id <- as.character(quant$target_id)
  }
  setkey(quant, 'target_id')
  
  # load bootstrap data if it exists
  bootspath <- file.path(fish_dir, 'quant_bootstraps.sf')
  numBoot <- 0
  if (file.exists(bootspath)) {
    boots <- fread(bootspath)
    target_ids <- as.character(names(boots))
    boots <- data.table(t(boots))
    setnames(boots, sapply(0:(ncol(boots)-1), function(i) paste('bs', i, sep='')))
    numBoot <- ncol(boots)
    boots[, target_id:=target_ids]
    setkey(boots, 'target_id')
    quant <- merge(quant, boots)
  }

  # load stats
  stats_file <- file.path(fish_dir, 'stats.tsv')
  ##
  # If the stats.tsv file exists, use that to get the 
  # number of observed fragments and effective lengths
  ##
  if (file.exists(stats_file)) {
    stats_tbl <- fread(stats_file)
    stats <- stats_tbl$V2
    names(stats) <- stats_tbl$V1
    stats_tbl <- stats_tbl[-1]
    setnames(stats_tbl, c('target_id', 'eff_length'))
    if (!is.character(quant$target_id)) {
      stats_tbl$target_id <- as.character(stats_tbl$target_id)
    }
    setkey(stats_tbl, 'target_id')
    quant <- merge(quant, stats_tbl)
    numProcessed <- stats[['numObservedFragments']]
  } else {
  ##
  # Otherwise, use hte provided value for the number of observed
  # fragments and compute the effective lengths given the provided
  # fragment length mean and standard deviation.
  ##
    max_len <- 1000
    print('Found no stats.tsv file when parsing old fish directory')
    print(sprintf('Generating fragment length distribution mu = %f, sd = %f, max len = %d', fallback_mu, fallback_sd, max_len))
    norm_counts <- get_norm_fl_counts(fallback_mu, fallback_sd, max_len)
    correction_factors <- get_eff_length_correction_factors(norm_counts, max_len)
    quant$eff_length <- get_eff_lengths(quant$length, correction_factors, max_len)
    if (fallback_num_reads < 0) {
      numProcessed <- round(sum(quant$est_counts))
    } else {
      msg <- paste("Fish dir =", fish_dir, ":: fallback # reads =", fallback_num_reads)
      message(msg)
      numProcessed <- fallback_num_reads
    }
    print(sprintf('Setting number of processed reads to %d', numProcessed))
  }
  

  # build the hdf5
  rhdf5::h5createFile(h5file)

  # counts are at root
  rhdf5::h5write(quant$est_counts, h5file, 'est_counts')

  # aux group has metadata about the run and targets
  rhdf5::h5createGroup(h5file, 'aux')
  rhdf5::h5write(numProcessed, h5file, 'aux/num_processed')
  rhdf5::h5write(numBoot, h5file, 'aux/num_bootstrap')
  rhdf5::h5write(quant$length, h5file, 'aux/lengths')
  rhdf5::h5write(quant$eff_length, h5file, 'aux/eff_lengths')
  rhdf5::h5write(quant$target_id, h5file, 'aux/ids')
  rhdf5::h5write('10', h5file, 'aux/index_version')
  rhdf5::h5write('sailfish', h5file, 'aux/kallisto_version')
  rhdf5::h5write(timestamp(prefix="", suffix=""), h5file, "aux/start_time")


  # bootstrap group has (.. wait for it ..) bootstrap data
  if (numBoot > 0) {
    rhdf5::h5createGroup(h5file, 'bootstrap')
    sapply(0:(numBoot-1), function(i) {
      bootid <- paste('bs', i, sep='')
      rhdf5::h5write(unlist(quant[, bootid, with=FALSE]),
              h5file, paste('bootstrap', bootid, sep='/'))
    })
  }

  if(utils::packageVersion('rhdf5') < "2.23.0") {
    rhdf5::H5close()
  } else {
    rhdf5::h5closeAll()
  }
  print(paste("Successfully converted sailfish / salmon results in", fish_dir, "to kallisto HDF5 format"))
}


#' Produce a collection of counts consistent with a truncated normal
#' distribution from 1 to max_len with the given mean and standard deviation
#' 
get_norm_fl_counts <- function(mean = mean, std = std, max_len) {
  dist <- vector(mode="numeric", max_len)
  totCount <- 10000
  
  # Evaluate the function at the point p
  kernel <- function(p) {
    invStd = 1.0 / std
    x = invStd * (p - mean)
    exp(-0.5 * x * x) * invStd
  }
  
  totMass = sum(unlist(lapply(1:max_len, kernel)))
  
  currDensity = 0.0
  if (totMass > 0.0) {
    for (i in 1:max_len) {
      currDensity = kernel(i)
      dist[i] = round(currDensity * totCount / totMass)
    }
  }
  return(dist)
}

#' Given a collection of counts for each fragment length, produce
#' a list of correction factors that should be applied to obtain the 
#' effective length of each transcript from its un-normalized length
#'
get_eff_length_correction_factors <- function(counts, max_len) {
  correctionFactors <- vector(mode="numeric", length=max_len)
  vals <- vector(mode="numeric", length=max_len)
  multiplicities <- vector(mode="numeric", length=max_len)
  
  multiplicities[1] = counts[1]
  
  for (i in 2:max_len) {
    ci = i
    pi = i-1
    v = counts[ci]
    vals[ci] = (v * ci) + vals[pi]
    multiplicities[ci] = v + multiplicities[pi]
    if (multiplicities[ci] > 0) {
      correctionFactors[ci] = vals[ci] / multiplicities[ci]
    }
  }
  return(correctionFactors)
}

#' Given a list of transcript lengths and correction factors, return
#' a list of effective lengths.
get_eff_lengths <- function(lengths, correction_factors, max_len) {
  
  eff_length <- function(orig_len) {
    correction_factor <- 
      if (orig_len >= max_len) { 
        correction_factors[max_len] 
      } else { 
        correction_factors[orig_len] 
      }
    
    eff_len <- orig_len - correction_factor + 1.0
    if (eff_len < 1.0) { eff_len <- orig_len }
    
    return(eff_len)
  }
  
  unlist(lapply(lengths, eff_length))
}


