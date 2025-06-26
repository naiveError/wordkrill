#' Compute Bootstrap Standard Errors
#' 
#' Computes bootstrap standard errors for document positions from a fitted
#' Wordkrill model.
#' 
#' This function computes a parametric bootstrap by resampling counts from the
#' fitted feature counts, refitting the model, and storing the document 
#' positions.
#' The standard deviations for each resampled document position are returned.
#' 
#' @param model a fitted Wordkrill model
#' @param B how many replications
#' @param verbose give progress updates
#' @param K dimension of the analysed Wordkrill model
#' @return list with standard errors for document positions
#' @author Benjamin Riesch
#' @importFrom stats rpois sd
#' @importFrom methods is
#' @export 
bootstrap_wordkrill <- function(model, B = 500, K = NULL, verbose = TRUE) {
  if (is.null(K)) stop("Please provide the number of dimensions K.")
  if (verbose) cat("Starting parametric bootstrap with", B, "replications...\n")
  
  I <- ncol(model$data)
  J <- nrow(model$data)
  
  # Create controls
  alpha <- model$alpha
  psi <- model$psi
  theta <- vector("list", K)
  beta <- vector("list", K)
  
  for (d in 1:K) {
    theta[[d]] <- model[[paste0("theta_", d)]]
    beta[[d]]  <- model[[paste0("beta_", d)]]
  }
  
  names(theta) <- paste0("theta_", 1:K)
  names(beta)  <- paste0("beta_", 1:K)
  
  start_params <- list(
    alpha = as.numeric(alpha),
    psi = as.numeric(psi),
    theta = theta,
    beta = beta
  )
  
  # reconstruct lambda matrix
  log_lambda <- outer(model$alpha, rep(1, J)) + outer(rep(1, I), model$psi)
  for (k in 1:K) {
    theta_k <- model[[paste0("theta_", k)]]
    beta_k <- model[[paste0("beta_", k)]]
    log_lambda <- log_lambda + outer(theta_k, beta_k)
  }
  lambda <- exp(log_lambda)
  
  # storage for bootstrap results
  boot_thetas <- vector("list", K)
  for (k in 1:K) {
    boot_thetas[[k]] <- matrix(NA, nrow = I, ncol = B)
  }
  
  for (b in 1:B) {
    if (verbose && b %% 10 == 0) cat("Bootstrap sample", b, "\n")
    
    # simulate new DFM
    sim_dtm <- matrix(rpois(I * J, lambda), nrow = I, ncol = J)
    
    # Fit model on simulated data
    fit <- try(wordkrill(sim_dtm, K = K, 
                                dir = model$dir, 
                                epsilon = model$epsilon, 
                                control = list(tol = 1e-08, 
                                               sigma = model$sigma, 
                                               startparams = NULL, 
                                               conv.check = 'll'),
                                verbose = verbose), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      for (k in 1:K) {
        boot_thetas[[k]][, b] <- fit[[paste0("theta_", k)]]
      }
    } else {
      if (verbose) cat("Fit failed on bootstrap sample", b, "\n")
    }
  }
  
  names(boot_thetas) <- paste0("theta_", 1:K)
  
  # compute standard errors
  theta_se <- lapply(boot_thetas, function(mat) {
    apply(mat, 1, sd)
  })
  
  return(theta_se)
}