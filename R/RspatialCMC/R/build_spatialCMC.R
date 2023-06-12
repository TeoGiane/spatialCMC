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

  # Set spatialcmc_home folder from RSPATIALCMC_HOME
  home_dir = Sys.getenv("RSPATIALCMC_HOME")
  spatialcmc_home = dirname(dirname(home_dir))

  # Build spatialCMC library
  build_dir = sprintf("%s/%s", spatialcmc_home, build_dirname)
  cat("*** Configuring spatialCMC ***\n")

  CONFIGURE = sprintf("mkdir -p %s && cd %s && cmake ..", build_dir, build_dir)
  tryCatch({
    system(CONFIGURE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  },
  error = function(e) {
    ERR_MSG <- "Something went wrong during configure. Failed with error:\n%s\n"
    stop(sprintf(ERR_MSG, as.character(e)))
  })

  # errlog = system(CONFIGURE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  # if(errlog != 0) { stop(sprintf("Something went wrong during configure: exit with status %d", errlog)) }

  # Make run_cmc executable
  cat("*** Building run_cmc executable ***\n")
  BUILD = sprintf("cd %s && make run_cmc -j%d", build_dir, nproc)
  errlog = system(BUILD)
  if(errlog != 0) { stop(sprintf("Something went wrong during build: exit with status %d", errlog)) }

  # Make run_mcmc executable
  cat("*** Building run_mcmc executable ***\n")
  BUILD = sprintf("cd %s && make run_mcmc -j%d", build_dir, nproc)
  errlog = system(BUILD)
  if(errlog != 0) { stop(sprintf("Something went wrong during build: exit with status %d", errlog)) }

  # Create .renviron file
  renviron = system.file("RspatialCMC.Renviron", package = "RspatialCMC")

  # Set SPATIALCMC_EXE environment variable
  cat("*** Setting RSPATIALCMC_EXE environment variable ***\n")
  write(x = sprintf("RSPATIALCMC_EXE=%s/run_cmc", build_dir), file = renviron, append = TRUE)
  # Set MCMC_EXE environment variable
  cat("*** Setting MCMC_EXE environment variable ***\n")
  write(x = sprintf("MCMC_EXE=%s/run_mcmc", build_dir), file = renviron, append = TRUE)

  # Parse .renviron file
  readRenviron(renviron)

  cat("spatialCMC successfully installed and environment variables set\n")
  # return(TRUE)
}
