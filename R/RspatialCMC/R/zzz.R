# Parse internal renviron file to set variables
.onAttach <- function(...) {
  readRenviron(system.file("RspatialCMC.Renviron", package = "RspatialCMC"))
  import_protobuf_messages()
}

# Unset variables on detaching
.onDetach <- function(...) {
  RProtoBuf::resetDescriptorPool()
  Sys.unsetenv("RSPATIALCMC_HOME")
  Sys.unsetenv("CMC_EXE")
  Sys.unsetenv("MCMC_EXE")
  Sys.unsetenv("POISREG_CMC_EXE")
  Sys.unsetenv("POISREG_MCMC_EXE")
  Sys.unsetenv("TBB_PATH")
}
