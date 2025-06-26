#' Get initial values for a Wordkrill model with K dimensions.
#'
#' This function is only called by model fitting routines and does therefore
#' not take a dtm classes. dtm is assumed to be in document by feature form.
#'
#' In the poisson form of the model incidental parameters (alpha) are set to
#' logged rowmeans divided by the first entry of rowmeans. Intercept (psi)
#' values are set to log(colmeans). These are subtracted from a the data matrix,
#' which is logged and decomposed by SVD. Word slope (beta) and document
#' position (theta) are estimated by rescaling SVD output.
#'
#' @param dtm a document by word matrix of counts
#' @param K number of dimensions of the corresponding wordkrill model
#' @param ort logical value indicating whether the initial values are
#' orthogonalised
#'
#' @return List with elements:
#' \item{alpha}{starting values of alpha parameters}
#' \item{psi}{starting values of psi parameters}
#' \item{theta_1...theta_K}{starting values of theta_1 to theta_K}
#' \item{theta_1...theta_K}{starting values of theta_1 to theta_K}
#' \item{beta_1...beta_K}{starting values of beta_1 to beta_K}
#'
#' @author Benjamin Riesch
#'
#' This is a multidimensional extension of the algorithm proposed by
#' Slapin and Proksch for obtaining initial values.
#' @references Slapin and Proksch (2008) 'A Scaling Model for Estimating
#' Time-Series Party Positions from Texts.' American Journal of Political
#' Science 52(3):705-772.
#' @references Lowe (2020) 'Austin: Do Things with Words'.
#' https://conjugateprior.github.io/austin
#' @references Riesch (2025) 'Wordkrill: Extending Wordfish into the
#' multidimensional political space'. https://doi.org/10.48550/arXiv.2506.20275
#'
#' @export
initialize_urkrill <- function(dtm, K, ort=F){
  tY      <- as.matrix(dtm)
  D       <- nrow(tY)
  V       <- ncol(tY)

  numword <- rep(1:V, each=D)
  numdoc  <- rep(1:D, V)

  dat           <- matrix(1, nrow=V*D, ncol=3)
  dat[,1]       <- as.vector(as.matrix(tY))
  dat[,2]       <- as.vector(numword)
  dat[,3]       <- as.vector(numdoc)
  dat           <- data.frame(dat)
  colnames(dat) <- c("y", "word", "doc")
  dat$word      <- factor(dat$word)
  dat$doc       <- factor(dat$doc)

  b0     <- log(colMeans(tY))
  rmeans <- rowMeans(tY)
  alpha  <- log(rmeans/rmeans[1])

  ystar  <- log(dat$y + 0.1) - alpha[dat$doc] - b0[dat$word]
  ystarm <- matrix(ystar, D, V, byrow=FALSE)
  res    <- svd(ystarm, nu=K)  # Adjusted for K dimensions

  theta_list <- list()
  beta_list  <- list()

  for (i in 1:K) {
    th <- as.vector(res$u[, i]) - res$u[1, i]
    b <- as.numeric(res$v[, i] * res$d[i])
    beta_list[[i]] <- b * sd(th)

    th <- (th-mean(th))/sd(th)

    # Gram-Schmidt-Orthogonalisierung
    if (i > 1 & ort==T) {
      for (j in 1:(i - 1)) {
        proj <- sum(th * theta_list[[j]]) / sum(theta_list[[j]]^2)
        th <- th - proj * theta_list[[j]]
      }
      th <- (th - mean(th)) / sd(th) # centralisation
    }

    theta_list[[i]] <- th
  }

  ## Name the lists appropriately
  names(theta_list) <- paste0("theta_", 1:K)
  names(beta_list)  <- paste0("beta_", 1:K)

  params <- list(
    alpha = as.numeric(alpha),
    psi = as.numeric(b0),
    beta = beta_list,
    theta =  theta_list
  )

  return(params)
}
