#' Bag of Little Bootstraps (BLB) Algorithm
#'
#' @description This function implements the Bag of Little Bootstraps (BLB) algorithm, which is designed for scalable resampling of large datasets.
#' The algorithm resamples subsets of data, computes point estimates, and then calculates the uncertainty based on these resamples.
#' It is particularly useful for estimating confidence intervals for large datasets while minimizing computational overhead.
#'
#' @param data A numeric vector representing the data from which the resamples are drawn.
#' @param beta A numeric value between 0 and 1 that controls the size of each bootstrap sample subset. The subset size is determined as `b = n^beta`, where `n` is the length of the data.
#' @param s An integer specifying the number of bootstrap subsets to generate.
#' @param r An integer specifying the number of Monte Carlo samples to draw for each subset.
#' @param u A function that calculates a point estimate from the data (e.g., mean, median, etc.).
#' @param xi A function that computes a statistic, such as a confidence interval, based on the differences between point estimates for different resamples.
#'
#' @return A numeric vector of length 2, representing the average values of the statistic `xi` computed for different subsets.
#'         This could, for example, represent the lower and upper bounds of a confidence interval.
#'
#' @examples
#' # Example usage:
#' data <- rnorm(1000)  # Generate a sample of 1000 random numbers
#' beta <- 0.5          # Subset size is n^beta, where n is length of data
#' s <- 100             # Number of subsets to generate
#' r <- 500             # Number of resampling iterations per subset
#' u <- mean            # Point estimate function (mean)
#' xi <- function(resamples, estimate) {
#'   quantile(resamples, c(0.025, 0.975))  # 95% confidence interval
#' }
#'
#' result <- BLB(data, beta, s, r, u, xi)
#' print(result)
#'
#' @export
BLB <- function(data, beta, s, r, u, xi) {
  point_estimate <- u(data)
  n <- length(data)
  b <- round(n^beta)
  xi_star <- vector("list", length = s)
  for (j in 1:s) {
    indices <- sample(x = 1:n, b, replace = FALSE)
    u_star <- numeric(r)
    for (k in 1:r) {
      resample <- sample(indices, n, replace = TRUE)
      u_star[k] <- u(data[resample])
    }
    xi_star[[j]] <- xi(u_star - point_estimate, point_estimate)
  }

  # Compute average values of xi for different subsets
  average <- numeric(2)
  for (i in 1:s) {
    average <- average + xi_star[[i]]
  }
  average <- average / s
  return(average)
}
