#' @useDynLib sumgp
#' @importFrom Rcpp sourceCpp
# ================================


# Function to create the the function K that will be used
# in a Gaussian process (Andrew's Version)
kernel_function <- function(squared_distance_matrix, nu, phi) {

  # Calculating the square matrix
  kernel_matrix <- exp(-squared_distance_matrix / (2 * phi^2)) / nu

  # Case nu = 0
  if(nu == 0 || nu > 1e13){
    kernel_matrix <- matrix(0, nrow = dim(squared_distance_matrix)[1],
                            ncol = dim(squared_distance_matrix)[2])
  }
  # Getting the kernel matrix
  return(kernel_matrix)
}
