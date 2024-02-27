#' Run Spatial Markov Chain Monte Carlo sampler
#'
#' In this light version, this function calls the \code{run_mcmc} executable from a subprocess via \code{\link[base]{system}} command.
#'
#' @param data A numeric vector or matrix of shape (n_samples, n_dim). These are the observations on which to fit the model.
#' @param geometry A \code{sfc} object of size \code{n_samples}. Default value is \code{NULL}.
#' @param hier_type A text string containing the enum label of the hierarchy to use in the algorithm.
#' @param hier_params A text string containing the hyperparameters of the hierarchy or a file name where the hyperparameters are stored.
#' A protobuf message of the corresponding type will be created and populated with the parameters.
#' @param mix_type A text string containing the enum label of the mixing to use in the algorithm.
#' @param mix_params A text string containing the hyperparameters of the mixing or a file name where the hyperparameters are stored.
#' A protobuf message of the corresponding type will be created and populated with the parameters.
#' @param algo_params A text string containing the hyperparameters of the algorithm or a file name where the hyperparameters are stored.
#' @param out_dir A string. If not \code{NULL}, is the folder where to store the output.
#' If \code{NULL}, a temporary directory will be created and destroyed after the sampling is finished.
#'
#' @return A list whose elements are Google Protocol Buffer Messages of type
#' \code{bayesmix::AlgorithmState}. Such object store a generic iteration of the MCMC chain.
#'
#' @export
run_mcmc <- function(data, geometry = NULL, hier_type, hier_params,
                     mix_type, mix_params, algo_params, out_dir = NULL) {

  # Check input types
  if (!is.numeric(data)) { stop("'data' parameter must be a numeric vector or matrix") }
  if (!is(geometry, "sfc") & !is.null(geometry)) { stop("'geometry' parameter must be an 'sfc' object, NULL otherwise") }
  if (!is.character(hier_type)){ stop("'hier_type' parameter must parameter must be a string") }
  if (!is.character(hier_params)) { stop("'hier_params' parameter must be a string") }
  if (!is.character(mix_type)){ stop("'mix_type' parameter must parameter must be a string") }
  if (!is.character(mix_params)) { stop("'mix_params' parameter must be a string") }
  if (!is.character(algo_params)) { stop("'algo_params' parameter must be a string") }
  if (!is.character(out_dir) & !is.null(out_dir)) { stop("'out_dir' parameter must be a string, NULL otherwise") }

  # Get MCMC_EXE and check if set
  MCMC_EXE = Sys.getenv("MCMC_EXE")
  if(MCMC_EXE == ""){
    stop("MCMC_EXE environment variable not set")
  }

  # Set-up template for run_mcmc command
  params = paste('--data-file %s --adj-matrix-file %s',
                 '--algo-params-file %s --hier-type %s',
                 '--hier-prior-file %s --mix-type %s',
                 '--mix-prior-file %s --chain-file %s')

  # Set run_mcmc command template
  RUN_CMD = paste(MCMC_EXE, params)

  # Use temporary directory if out_dir is not set
  if(is.null(out_dir)) {
    rngstr <- paste0(sample(letters, 7, replace=T), collapse = "")
    out_dir = sprintf("%s/Rtmp-%s", dirname(MCMC_EXE), rngstr); dir.create(out_dir, showWarnings = F)
    remove_out_dir = TRUE
  } else {
    remove_out_dir = FALSE
  }

  # Compute adjacency matrix from geometry
  if (!is.null(geometry)) {
    adj_matrix <- spdep::nb2mat(spdep::poly2nb(geometry, queen = F), style = "B", zero.policy = T)
  }

  # Prepare files for data and outcomes
  data_file = paste0(out_dir,'/data.csv'); file.create(data_file)
  adj_matrix_file = paste0(out_dir, '/adj_matrix.csv'); file.create(adj_matrix_file)
  chain_name <- sprintf('/chain_%s.recordio', format(Sys.time(), "%Y%m%d-%H%M"))
  chain_file = paste0(out_dir, chain_name); file.create(chain_file)

  # Prepare protobuf configuration files
  hier_params_file = maybe_print_to_file(hier_params, "hier_params", out_dir)
  mix_params_file = maybe_print_to_file(mix_params, "mix_params", out_dir)
  algo_params_file = maybe_print_to_file(algo_params, "algo_params", out_dir)

  # Set-up NULL filenames for arg-parse
  EMPTYSTR = '\\"\\"'
  write.table(data, file = data_file, sep = ",", col.names = F, row.names = F)
  if (is.null(geometry)) {
    adj_matrix_file <- EMPTYSTR
  } else {
    write.table(adj_matrix, file = adj_matrix_file, sep = ",", row.names = F, col.names = F)
  }

  # Resolve run_mcmc command
  CMD = sprintf(RUN_CMD, data_file, adj_matrix_file, algo_params_file,
                hier_type, hier_params_file, mix_type, mix_params_file, chain_file)

  # Execute run_mcmc
  errlog <- system(CMD)
  if(errlog != 0L){
    unlink(out_dir, recursive = TRUE)
    errmsg <- "Something went wrong: 'run_mcmc()' exit with status %d"
    stop(sprintf(errmsg, errlog))
  }

  # Manage return object - Serialized MCMC chain
  tryCatch({
    RspatialCMC::import_protobuf_messages()
    chain <- RProtoBuf::read(spatialcmc.MCMCChain, chain_file)
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
