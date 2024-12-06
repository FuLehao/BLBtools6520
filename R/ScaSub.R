#' Title: Scalable Subsampling Algorithm
#'
#' @description
#' This function implements a scalable subsampling algorithm, which divides the dataset into overlapping or non-overlapping blocks
#' of size `b` (controlled by `beta`) and computes statistics on each block. The function is particularly useful for estimating
#' confidence intervals for large datasets with reduced computational costs compared to full bootstrapping methods.
#'
#' @param data A numeric vector representing the dataset to be subsampled.
#' @param beta A numeric value between 0 and 1 that determines the size of each subsample block as `b = n^beta`, where `n` is the length of the dataset.
#' @param c1 A numeric value greater than 0 that controls the overlap between blocks. The stride length for block selection is `h = c1 * b`.
#' @param u A function that computes a point estimate (e.g., mean, median) from the data.
#' @param xi A function that calculates a statistic (e.g., confidence interval) based on the distribution of the root distances.
#'
#' @return A numeric value or vector representing the result of applying the statistic function `xi` to the root distances.
#'         For example, it may return confidence intervals.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' data <- rnorm(1000)  # Generate a dataset of 1000 random numbers
#' beta <- 0.5          # Subsample block size is n^beta, where n is the data length
#' c1 <- 0.2            # Stride length as a fraction of block size
#' u <- mean            # Point estimate function (mean)
#' xi <- function(root_dist, estimate) {
#'   quantile(root_dist, c(0.025, 0.975))  # 95% confidence interval
#' }
#'
#' result <- ScaSub(data, beta, c1, u, xi)
#' print(result)

ScaSub <- function(data, beta, c1, u, xi) {
  point_estimate <- u(data)
  n <- length(data)
  b <- round(n^beta)
  h <- round(c1 * b)
  q <- floor((n - b) / h) + 1
  u_star <- numeric(q)
  for (i in 1:q) {
    u_star[i] <- u(data[((i - 1) * h + 1):((i - 1) * h + b)])
  }
  root_dist <- sqrt(b / n) * (u_star - point_estimate)
  return(xi(root_dist, point_estimate))
}
