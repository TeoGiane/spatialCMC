#' Builds the spatialCMC executable
#'
#' After the build, if no error has occurred, it saves the path into the \code{RSPATIALCMC_EXE} environment variable.
#' Such variable is defined only when this package is loaded in the R session.
#'
#' @param nproc Number of processes to use for parallel compilation. Thanks to \code{parallel} package,
#' this parameter defaults to half of the available processes (through \code{\link[parallel]{detectCores}} function)
#' @param build_dirname Name for the build directory of \code{spatialCMC}. Default is "build".
#' @return \code{TRUE} if build is successful,it raises errors otherwise
#'
#' @export
build_spatialCMC <- function(nproc = ceiling(parallel::detectCores()/2), build_dirname = "build") {

  # Check input types
  if(!is.numeric(nproc)) { stop("nproc must be a number") }
  if(!is.character(build_dirname)) { stop("build_dirname must be a string") }

  # Get .Renviron file from package
  renviron = system.file("RspatialCMC.Renviron", package = "RspatialCMC")

  # Set spatialcmc_home folder from RSPATIALCMC_HOME
  readRenviron(renviron)
  home_dir = Sys.getenv("RSPATIALCMC_HOME")
  if(home_dir == ""){
    stop("RSPATIALCMC_HOME environment variable is not set")
  }
  spatialcmc_home = dirname(dirname(home_dir))

  # Create build/ subdirectory
  build_dir = sprintf("%s/%s", spatialcmc_home, build_dirname)
  dir.create(build_dir, showWarnings = F)

  # Configure spatialCMC
  cat("*** Configuring spatialCMC ***\n")
  CONFIGURE = 'cmake .. -G "Unix Makefiles"'
  errlog <- withr::with_dir(build_dir, system(CONFIGURE, ignore.stderr = TRUE))
  if(errlog != 0L){
    errmsg <- "Something went wrong during configure: command '%s' exit with status %d"
    stop(sprintf(errmsg, CONFIGURE, errlog))
  }
  cat("\n")

  # Make run_cmc executable
  cat("*** Building run_cmc executable ***\n")
  BUILD = sprintf("make run_cmc -j%d", nproc)
  errlog <- withr::with_dir(build_dir, system(BUILD))
  if (errlog != 0L) {
    errmsg <- "Something went wrong during build: command '%s' exit with status %d"
    stop(sprintf(errmsg, BUILD, errlog))
  }
  cat("\n")

  # Make run_mcmc executable
  cat("*** Building run_mcmc executable ***\n")
  BUILD = sprintf("make run_mcmc -j%d", nproc)
  errlog <- withr::with_dir(build_dir, system(BUILD))
  if (errlog != 0L) {
    errmsg <- "Something went wrong during build: command '%s' exit with status %d"
    stop(sprintf(errmsg, BUILD, errlog))
  }
  cat("\n")

  # Make run_poisreg_mcmc executable
  cat("*** Building run_poisreg_mcmc executable ***\n")
  BUILD = sprintf("make run_poisreg_mcmc -j%d", nproc)
  errlog <- withr::with_dir(build_dir, system(BUILD))
  if (errlog != 0L) {
    errmsg <- "Something went wrong during build: command '%s' exit with status %d"
    stop(sprintf(errmsg, BUILD, errlog))
  }
  cat("\n")

  # Make run_poisreg_cmc executable
  cat("*** Building run_poisreg_cmc executable ***\n")
  BUILD = sprintf("make run_poisreg_cmc -j%d", nproc)
  errlog <- withr::with_dir(build_dir, system(BUILD))
  if (errlog != 0L) {
    errmsg <- "Something went wrong during build: command '%s' exit with status %d"
    stop(sprintf(errmsg, BUILD, errlog))
  }
  cat("\n")

  # Set RSPATIALCMC_EXE environment variable
  cat("*** Setting CMC_EXE environment variable ***\n")
  write(x = sprintf("CMC_EXE=%s/run_cmc", build_dir), file = renviron, append = TRUE)

  # Set MCMC_EXE environment variable
  cat("*** Setting MCMC_EXE environment variable ***\n")
  write(x = sprintf("MCMC_EXE=%s/run_mcmc", build_dir), file = renviron, append = TRUE)

  # Set POISREG_MCMC_EXE environment variable
  cat("*** Setting POISREG_MCMC_EXE environment variable ***\n")
  write(x = sprintf("POISREG_MCMC_EXE=%s/run_poisreg_mcmc", build_dir), file = renviron, append = TRUE)

  # Set POISREG_CMC_EXE environment variable
  cat("*** Setting POISREG_CMC_EXE environment variable ***\n")
  write(x = sprintf("POISREG_CMC_EXE=%s/run_poisreg_cmc", build_dir), file = renviron, append = TRUE)

  # Set TBB_PATH environment variable (for Windows)
  cat("*** Setting TBB_PATH environment variable ***\n")
  tbb_path = sprintf('%s/lib/bayesmix-src/lib/_deps/math-src/lib/tbb/', spatialcmc_home)
  write(x = sprintf('TBB_PATH=%s', tbb_path), file = renviron, append = TRUE)
  cat("\n")

  # Parse .Renviron file to get environment variables
  readRenviron(renviron)
  cat("Successfully installed spatialCMC\n")
}
