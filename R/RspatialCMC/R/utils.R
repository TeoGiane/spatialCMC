# If maybe_proto is a file, returns the file name. If maybe_proto is a string representing a message,
# prints the message to a file and returns the file name.
maybe_print_to_file <- function(maybe_proto, proto_name = NULL, out_dir = NULL) {
  if(file.exists(maybe_proto)){
    return(maybe_proto)
  }
  proto_file = sprintf("%s/%s.asciipb", out_dir, proto_name)
  write(maybe_proto, file = proto_file)
  return(proto_file)
}

#' Reads an MCMC chain from a protobuf file and returns it as a list of serialized messages.
#'
#' @param filename The path to the protobuf file containing the MCMC chain.
#' @return A list of serialized protobuf messages representing the MCMC chain.
#' 
#' @export
read_mcmc_chain <- function(filename) {
  RspatialCMC::import_protobuf_messages()
  chain <- sapply(RProtoBuf::read(spatialcmc.MCMCChain, filename)$state,
                  function(x) {RProtoBuf::serialize(x, NULL)})
  return(chain)
}

#' Convert a vector of values into sequential numeric labels
#'
#' This function accepts a vector of character, factor, numeric, or logical
#' values and maps each unique non-NA value (in order of first appearance)
#' to a sequential integer starting at `start`. NAs in the input are preserved
#' as `NA` in the output. The function returns a list containing the new
#' integer vector and the vector of unique values (including `NA` if present)
#' in order of appearance; the latter is useful for plot labels.
#'
#' @param x A vector (character, factor, numeric, logical) to relabel.
#' @param start Integer scalar giving the first value for the sequence (default 0).
#' @return A list with components:
#'   - `new_values`: integer vector of sequential labels (start..start+k-1). NAs are preserved as `NA`.
#'   - `old_labels`: vector of the unique values from `x` in order of appearance (may include `NA`).
#' @export
relabel <- function(x, start = 0L) {
  if (!is.numeric(start) || length(start) != 1 || is.na(start) || start != as.integer(start) || start < 0) {
    stop("`start` must be a single integer >= 0")
  }
  start <- as.integer(start)

  if (is.null(x)) return(list(new = integer(0), values = vector(mode = "list", length = 0)))

  if (is.matrix(x) || is.data.frame(x)) x <- as.vector(x)

  # For factors use character representation so labels are meaningful
  if (is.factor(x)) {
    x_vals <- as.character(x)
  } else {
    x_vals <- x
  }

  unique_vals <- unique(x_vals)

  # Prepare mapping for non-NA values only; keep NA in `values` but do not
  # assign a numeric label to NA (preserve as NA in `new`).
  unique_non_na <- unique_vals[!is.na(unique_vals)]
  n_non_na <- length(unique_non_na)

  new <- rep(NA_integer_, length(x_vals))
  if (n_non_na > 0) {
    # match each non-NA element to the index within unique_non_na, then shift by start-1
    non_na_idx <- which(!is.na(x_vals))
    if (length(non_na_idx) > 0) {
      matched <- match(x_vals[non_na_idx], unique_non_na)
      new[non_na_idx] <- as.integer(matched + start - 1L)
    }
  }

  return(list(new_values = new, old_labels = unique_vals))
}
