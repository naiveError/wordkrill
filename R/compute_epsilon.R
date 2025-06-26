#' Calculates epsilon so that the estimated target variables (mean, variance
#' and covariance) lie within an appropriate range around the target values.
#'
#' @param dtm a quanteda DTM object
#' @param alpha control parameters to determine the strictness of the epsilon
#' environment
#' @param norm logical value indicating normality distribution of the thetas is
#' assumed or not
#' @param K dimension parameter for the model to be estimated
#' @return An object of class numeric, containing:
#'
#' \item{epsilon}{estimated optimal epsilon for the model to be estimated}
#'
#' @author Benjamin Riesch
#' @references Riesch (2025) 'Wordkrill: Extending Wordfish into the
#' multidimensional political space'. https://doi.org/10.48550/arXiv.2506.20275
#'
#' @importFrom stats cov qchisq qnorm rnorm var cor
#'
#' @export
compute_epsilon <- function(dtm, alpha = 0.05, norm = T, K=2) {
  tY <- as.matrix(dtm)
  n <- nrow(tY)  # number of documents
  epsilon <- NULL

  # mean
  z <- qnorm(1 - alpha / 2)
  epsilon_mean <- z / sqrt(797)
  epsilon <- epsilon_mean

  if(norm==T){
    # variance
    df <- n - 1
    chi_lower <- qchisq(1 - alpha / 2, df)
    chi_upper <- qchisq(alpha / 2, df, lower.tail = FALSE)

    var_lower <- df / chi_upper
    var_upper <- df / chi_lower
    epsilon_var <- max(abs(var_lower - 1), abs(var_upper - 1))
    epsilon <- c(epsilon, epsilon_var)
    if(K>1){
      z <- qnorm(1 - alpha / 2)
      epsilon_cov <- z / sqrt(n-1)
      epsilon <- epsilon <- c(epsilon, epsilon_cov)
    }else{}

  }else{}

  # Return the most conservative epsilon
  return(max(epsilon))
}
