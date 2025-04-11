#' Run Spatial Markov Chain Monte Carlo sampler
#'
#' In this light version, this function calls the \code{run_mcmc} executable from a subprocess via \code{\link[base]{system}} command.
#'
#' @param data A numeric vector or matrix of shape (n_samples, n_dim). These are the observations on which to fit the model.
#' @param offset A numeric vector containing the fixed offsets used in the Poisson regression.
#' @param cov_matrix A matrix with n_samples rows and p columns, each row is the vector of covariates associated to i-th observation
#' @param geometry A \code{sfc} object of size \code{n_samples}. The spatial structure to consider.
#' @param shard_allocation A numeric vector of size \code{n_samples}.
#' This represents the shard in which the corresponding datum will be assigned.
#' @param hier_params A text string containing the hyperparameters of the hierarchy or a file name where the hyperparameters are stored.
#' A protobuf message of the corresponding type will be created and populated with the parameters.
#' @param mix_params A text string containing the hyperparameters of the mixing or a file name where the hyperparameters are stored.
#' A protobuf message of the corresponding type will be created and populated with the parameters.
#' @param algo_params A text string containing the hyperparameters of the algorithm or a file name where the hyperparameters are stored.
#' @param out_dir A string. If not \code{NULL}, is the folder where to store the output.
#' If \code{NULL}, a temporary directory will be created and destroyed after the sampling is finished.
#'
#' @return A list whose elements are Google Protocol Buffer Messages of type
#' \code{spatialcmc::PoisRegAlgorithmState}. Such object store a generic iteration of the MCMC chain.
#'
#' @export
pois_reg_cmc <- function(data, offset, cov_matrix, geometry, shard_allocation,
                         hier_params, mix_params, algo_params, out_dir = NULL) {

  # Check input types
  if (!is.numeric(data)) { stop("'data' parameter must be a numeric vector or matrix") }
  if (!is.numeric(offset)) { stop("'offset' parameter must be a numeric vector or matrix") }
  if (!is.numeric(cov_matrix)) { stop("'cov_matrix' parameter must be a numeric vector or matrix") }
  if (!is(geometry, "sfc")) { stop("'geometry' parameter must be an 'sfc' object") }
  if (!is.numeric(shard_allocation)) { stop("'shard_allocation' parameter must be a numeric vector") }
  if (!is.character(hier_params)) { stop("'hier_params' parameter must be a string") }
  if (!is.character(mix_params)) { stop("'mix_params' parameter must be a string") }
  if (!is.character(algo_params)) { stop("'algo_params' parameter must be a string") }
  if (!is.character(out_dir) & !is.null(out_dir)) { stop("'out_dir' parameter must be a string, NULL otherwise") }

  # Check compatibility between matrices
  if((length(data) != nrow(cov_matrix)) || (length(data) != length(offset)) || (length(data) != length(geometry))) {
    stop("Input data does not match in size")
  }

  # Get POISREG_MCMC_EXE and check if set
  POISREG_CMC_EXE = Sys.getenv("POISREG_CMC_EXE")
  if(POISREG_CMC_EXE == ""){
    stop("POISREG_CMC_EXE environment variable not set")
  }

  # Set-up template for run_poisreg_mcmc command
  params = paste('--data-file %s --hier-cov-file %s',
                 '--adj-matrix-file %s --shard-assignment-file %s',
                 '--algo-params-file %s --hier-prior-file %s',
                 '--mix-prior-file %s --chain-file %s')

  # Set run_mcmc command template
  RUN_CMD = paste(POISREG_CMC_EXE, params)

  # Use temporary directory if out_dir is not set
  if(is.null(out_dir)) {
    rngstr <- paste0(sample(letters, 7, replace=T), collapse = "")
    out_dir = sprintf("%s/Rtmp-%s", dirname(POISREG_CMC_EXE), rngstr); dir.create(out_dir, showWarnings = F)
    remove_out_dir = TRUE
  } else {
    remove_out_dir = FALSE
  }

  # Compute adjacency matrix from geometry
  adj_matrix <- spdep::nb2mat(spdep::poly2nb(geometry, queen = F), style = "B", zero.policy = T)

  # Add offsets to covariance matrix
  hier_cov <- cbind(offset, cov_matrix)

  # Prepare files for data and outcomes
  data_file = paste0(out_dir,'/data.csv'); file.create(data_file)
  hier_cov_file = paste0(out_dir, '/hier_cov.csv'); file.create(hier_cov_file)
  adj_matrix_file = paste0(out_dir, '/adj_matrix.csv'); file.create(adj_matrix_file)
  shard_allocation_file = paste0(out_dir,'/shard_allocation.csv'); file.create(shard_allocation_file)
  chain_name = sprintf('/chain_%s.recordio', format(Sys.time(), "%Y%m%d-%H%M"))
  chain_file = paste0(out_dir, chain_name); file.create(chain_file)

  # Prepare protobuf configuration files
  hier_params_file = maybe_print_to_file(hier_params, "hier_params", out_dir)
  mix_params_file = maybe_print_to_file(mix_params, "mix_params", out_dir)
  algo_params_file = maybe_print_to_file(algo_params, "algo_params", out_dir)

  # Set-up NULL filenames for arg-parse
  write.table(data, file = data_file, sep = ",", col.names = F, row.names = F)
  write.table(hier_cov, file = hier_cov_file, sep = ",", row.names = F, col.names = F)
  write.table(adj_matrix, file = adj_matrix_file, sep = ",", row.names = F, col.names = F)
  write.table(shard_allocation, file = shard_allocation_file, sep = ",", row.names = F, col.names = F)

  # Resolve run_mcmc command
  CMD = sprintf(RUN_CMD, data_file, hier_cov_file, adj_matrix_file, shard_allocation_file,
                algo_params_file, hier_params_file, mix_params_file, chain_file)

  # Execute run_poisreg_mcmc
  errlog <- system(CMD)
  if(errlog != 0L){
    unlink(out_dir, recursive = TRUE)
    errmsg <- "Something went wrong: 'run_poisreg_cmc()' exit with status %d"
    stop(sprintf(errmsg, errlog))
  }

  # Manage return object - Serialized MCMC chain
  tryCatch({
    RspatialCMC::import_protobuf_messages()
    chain <- RProtoBuf::read(spatialcmc.PoisRegMCMC, chain_file)
    chain <- lapply(chain$state, function(s){ RProtoBuf::serialize(s, NULL) })
  },
  error = function(e) {
    # Print error
    ERR_MSG <- "Parsing %s failed with error:\n%s"
    message(sprintf(ERR_MSG, chain_file, as.character(e)))
    # Clean up and return NULL
    unlink(out_dir, recursive = TRUE)
    return(NULL)
  })

  # Clean temporary files and return
  if(remove_out_dir) {
    unlink(out_dir, recursive = TRUE)
  }

  # Return deserialized chain
  return(chain)
}
