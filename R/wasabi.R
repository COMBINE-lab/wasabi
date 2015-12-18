#' Convert sailfish/salmon results for one or more samples to kallisto HDF5
#'
#' @param fish_dirs a character vector of length greater than one where each
#' string points to a sailfish/salmon output directory
#' @export
prepare_fish <- function(fish_dirs, force=FALSE) {
  testdir <- fish_dirs[1]
  if (file.exists(file.path(testdir, "aux", "meta_info.json"))) {
    sapply(fish_dirs, fish_to_hdf5, force=force)
  } else {
    sapply(fish_dirs, fish_to_hdf5_old, force=force)
  }
  fish_dirs
}

#' Convert sailfish results in new format
#' SF ver >= 0.9.0, Salmon ver >= 0.6.0 for one sample to
#' kallisto HDF5
#'
#' @param fish_dir path to a sailfish output directory
#' @param force if TRUE re-create the h5 file even if it exists
fish_to_hdf5 <- function(fish_dir, force) {
  library(data.table)
  library(rjson)

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

  # get all of the meta info
  minfo <- rjson::fromJSON(file=file.path(fish_dir, "aux", "meta_info.json"))

  # load bootstrap data if it exists
  auxPath <- file.path(fish_dir, 'aux')
  numBoot <- minfo$num_bootstraps
  if (numBoot > 0) {
    bootCon <- gzcon(file(file.path(auxPath, 'bootstrap', 'bootstraps.gz'), "rb"))
    boots <- readBin(bootCon, "double", n = minfo$num_targets * minfo$num_bootstraps)
    close(bootCon)

    # rows are transcripts, columns are bootstraps
    dim(boots) <- c(minfo$num_targets, minfo$num_bootstraps)
  }

  # load stats
  numProcessed <- minfo$num_processed

  # build the hdf5
  library(rhdf5)
  h5createFile(h5file)

  # counts are at root
  h5write(quant$est_counts, h5file, '/est_counts')

  # aux group has metadata about the run and targets
  h5createGroup(h5file, 'aux')
  h5write(numProcessed, h5file, 'aux/num_processed')
  h5write(numBoot, h5file, 'aux/num_bootstrap')
  h5write(quant$length, h5file, 'aux/lengths')
  h5write(quant$eff_length, h5file, 'aux/eff_lengths')
  h5write(quant$target_id, h5file, 'aux/ids')
  h5write('10', h5file, 'aux/index_version')
  h5write('sailfish', h5file, 'aux/kallisto_version')
  h5write(timestamp(prefix="", suffix=""), h5file, "aux/start_time")

  # bootstrap group has (.. wait for it ..) bootstrap data
  if (numBoot > 0) {
    h5createGroup(h5file, 'bootstrap')
    sapply(0:(numBoot-1), function(i) {
      bootid <- paste('bs', i, sep='')
      h5write(unlist(boots[,i+1]),
              h5file, paste('bootstrap', bootid, sep='/'))
    })
  }

  bootCon <- gzcon(file(file.path(auxPath, 'fld.gz'), "rb"))
  fld <- readBin(bootCon, "int", n=minfo$frag_dist_length)
  close(bootCon)
  h5write(fld, h5file, 'aux/fld')

  bObsCon <- gzcon(file(file.path(auxPath, 'observed_bias.gz'), "rb"))
  bObs <- readBin(bObsCon, "int", n=minfo$num_bias_bins)
  close(bObsCon)
  h5write(bObs, h5file, 'aux/bias_observed')

  bExpCon <- gzcon(file(file.path(auxPath, 'expected_bias.gz'), "rb"))
  bExp <- readBin(bObsCon, "double", n=minfo$num_bias_bins)
  close(bExpCon)
  h5write(bExp, h5file, 'aux/bias_normalized')

  H5close()
  print(paste("Successfully converted sailfish / salmon results in", fish_dir, "to kallisto HDF5 format"))
}

#' Convert sailfish results in the older format
#' sf ver <= 0.8.0, salmon ver <= 0.5.1 for one sample
#' to kallisto HDF5
#'
#' @param fish_dir path to a sailfish output directory
#' @param force if TRUE re-create the h5 file even if it exists
fish_to_hdf5_old <- function(fish_dir, force) {
  library(data.table)

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
  setkey(quant, 'target_id')


  # load bootstrap data if it exists
  bootspath <- file.path(fish_dir, 'quant_bootstraps.sf')
  numBoot <- 0
  if (file.exists(bootspath)) {
    boots <- fread(bootspath)
    target_ids <- names(boots)
    boots <- data.table(t(boots))
    setnames(boots, sapply(0:(ncol(boots)-1), function(i) paste('bs', i, sep='')))
    numBoot <- ncol(boots)
    boots[, target_id:=target_ids]
    setkey(boots, 'target_id')
    quant <- merge(quant, boots)
  }

  # load stats
  stats_tbl <- fread(file.path(fish_dir, 'stats.tsv'))
  stats <- stats_tbl$V2
  names(stats) <- stats_tbl$V1
  stats_tbl <- stats_tbl[-1]
  setnames(stats_tbl, c('target_id', 'eff_length'))
  setkey(stats_tbl, 'target_id')
  quant <- merge(quant, stats_tbl)

  numProcessed <- stats[['numObservedFragments']]

  # build the hdf5
  library(rhdf5)
  h5createFile(h5file)

  # counts are at root
  h5write(quant$est_counts, h5file, 'est_counts')

  # aux group has metadata about the run and targets
  h5createGroup(h5file, 'aux')
  h5write(numProcessed, h5file, 'aux/num_processed')
  h5write(numBoot, h5file, 'aux/num_bootstrap')
  h5write(quant$length, h5file, 'aux/lengths')
  h5write(quant$eff_length, h5file, 'aux/eff_lengths')
  h5write(quant$target_id, h5file, 'aux/ids')
  h5write('10', h5file, 'aux/index_version')
  h5write('sailfish', h5file, 'aux/kallisto_version')
  h5write(timestamp(prefix="", suffix=""), h5file, "aux/start_time")


  # bootstrap group has (.. wait for it ..) bootstrap data
  if (numBoot > 0) {
    h5createGroup(h5file, 'bootstrap')
    sapply(0:(numBoot-1), function(i) {
      bootid <- paste('bs', i, sep='')
      h5write(unlist(quant[, bootid, with=FALSE]),
              h5file, paste('bootstrap', bootid, sep='/'))
    })
  }

  H5close()
  print(paste("Successfully converted sailfish / salmon results in", fish_dir, "to kallisto HDF5 format"))
}
