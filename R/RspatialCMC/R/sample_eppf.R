#' Sample from the EPPF of a (spatial) Product Partition Model
#'
#' In this light version, this function calls the \code{run_mcmc} executable from a subprocess via \code{\link[base]{system}} command.
#'
#' @param n A double representing the number of data.
#' @param geometry A \code{sfc} object of size \code{n}. Default value is \code{NULL}.
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
sample_eppf <- function(n, geometry = NULL, mix_type,
                        mix_params, algo_params, out_dir = NULL) {

  # Define fictional dataset, Empty hierarchy with no prior
  data <- numeric(n)
  hier_type <- "Empty"
  hier_params <- ""

  # Call run_mcmc function
  return(run_mcmc(data, geometry, hier_type, hier_params,
                  mix_type, mix_params, algo_params, out_dir))
}
