#' Import Protocol Buffers Descriptors for spatialCMC
#'
#' This utility loads in the workspace the protocol buffers descriptors defined
#' in the \code{bnp-downscaling} library, via \code{RProtoBuf} package. These
#' descriptors can be used to deserialize the output of \code{\link{run_mcmc}}
#' function.
#'
#' @return NULL
#'
#' @export
import_protobuf_messages <- function() {

  # Deduce BUILD_DIR from RSPATIALCMC_EXE variable in RspatialCMC.Renviron
  BUILD_DIR <- dirname(Sys.getenv("RSPATIALCMC_EXE"))
  SPATIALCMC_HOME <- dirname(dirname(Sys.getenv("RSPATIALCMC_HOME")))

  # Deduce protocol buffer proto paths
  bayesmix_protoPath <- sprintf("%s/_deps/bayesmix-src/src/proto", BUILD_DIR)
  spatialcmc_protoPath <- sprintf("%s/src/proto", SPATIALCMC_HOME)

  # Import protocol buffer message descripts via RProtoBuf
  RProtoBuf::readProtoFiles2(protoPath = bayesmix_protoPath)
  RProtoBuf::readProtoFiles2(protoPath = spatialcmc_protoPath)
}
