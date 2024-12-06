#' @title Parallelized Bag of Little Bootstraps (BLB) Algorithm
#' @description Implements a parallelized version of the BLB algorithm.
#' @param data A numeric vector representing the input data.
#' @param beta A numeric value between 0 and 1, controlling the subset size.
#' @param s Number of bootstrap subsets.
#' @param r Monte Carlo samples per subset.
#' @param u Function to compute point estimates.
#' @param xi Function to calculate the statistic (e.g., confidence intervals).
#' @return A numeric vector representing averaged statistics from bootstrap subsets.
#' @import parallel foreach doParallel
#' @export
BLBpar <- function(data, beta, s, r, u, xi) {
  point_estimate <- u(data)
  n <- length(data)
  b <- round(n^beta)

  # Use multi-cores to conduct parallel computing
  ncores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  # Ensure cluster is stopped on exit
  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()  # Reset to sequential backend
  }, add = TRUE)

  xi_star <- foreach::foreach(j = 1:s) %dopar% {
    indices <- sample(x = 1:n, b, replace = FALSE)
    u_star <- numeric(r)
    for (k in 1:r) {
      resample <- sample(indices, n, replace = TRUE)
      u_star[k] <- u(data[resample])
    }
    xi(u_star - point_estimate, point_estimate)
  }

  # Compute average values of xi computed for different subsets
  average <- numeric(2)
  for (i in 1:s) {
    average <- average + xi_star[[i]]
  }
  average <- average / s
  return(average)
}
