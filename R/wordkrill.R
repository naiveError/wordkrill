#' K-dimensional Wordkrill Model with Conditional Maximum Likelihood
#'
#' This function estimates a Wordkrill aka a K-dimensional Wordfish model,
#' based on a Document-Term Matrix (DTM). It outputs document positions
#' (`theta_1`,...,`theta_K`) and word parameters (`beta_1`,..., `beta_K`),
#' besides the usual wordfish parameters estimated with conditional maximum
#' likelihood.
#'
#' @param dtm A quanteda DTM object
#' @param K integer to select the dimensionality of the Wordkrill model
#' @param dir set global identification by forcing \code{theta_k[dir[1]]} <
#' \code{theta_k[dir[2]]}
#' @param epsilon numeric to  determine the strictness of the epsilon
#' environment
#' @param control list of estimation options (tol: tolerance for likelihood
#' procedure; sigma: regularization parameter for betas in poisson form;
#' startparams: initial values of the parameters to be estimated; conv.check:
#' selecting the conversion check)
#' @param verbose logical, produce a running commentary
#' @param asympSE logical, compute asymptotcia standard errors
#'
#' @return An object of class wordkrill. This is a list containing:
#' \item{dir}{global identification of the dimension}
#' \item{theta_1 ... theta_K}{document positions dim 1 to K}
#' \item{alpha}{document fixed effects}
#' \item{beta_1 ... beta_K}{word slope parameters dim 1 to K}
#' \item{psi}{word fixed effects}
#' \item{docs}{names of the documents}
#' \item{features}{names of features}
#' \item{sigma}{regularization parameter for betas in Poisson form}
#' \item{ll}{final log likelihood}
#' \item{theta_SE_1 ... theta_SE_K}{standard errors for document position dim 1
#' to K}
#' \item{data}{the original data}
#'
#' @author Benjamin Riesch
#'
#' @references Slapin and Proksch (2008) 'A Scaling Model for Estimating
#' Time-Series Party Positions from Texts.' American Journal of Political
#' Science 52(3): 705-772.
#' @references Lowe (2020) 'Austin: Do Things with Words'.
#' https://conjugateprior.github.io/austin
#' @references Riesch (2025) 'Wordkrill: Extending Wordfish into the
#' multidimensional political space'. https://doi.org/10.48550/arXiv.2506.20275
#'
#' @examples
#'
#' simulation <- sim_wordkrill(docs = 10, vocab =20, K=2, doclen =5000,
#'                             dist = "normal")
#'
#' test_model <- wordkrill(dtm=simulation$Y, dir=c(1,2), K=2,
#'                                epsilon = 0.1,
#'                                control=list(tol=1e-08,
#'                                             sigma=3,
#'                                             startparams=NULL,
#'                                             conv.check='ll'),
#'                                verbose=TRUE, asympSE=TRUE)
#'
#' print(summary(object = test_model, top_n = 15))
#'
#' @importFrom alabama constrOptim.nl
#' @importFrom quanteda convert featnames docnames as.dfm
#' @importFrom methods is
#' @importFrom stats dmultinom cor median optim optimize var rnorm cov
#' @importFrom utils flush.console
#' @importFrom numDeriv hessian
#' @export
wordkrill <- function(dtm, K = NULL,
                             dir = c(1, 2),
                             epsilon = 1,
                             control = list(tol = 1e-08,
                                            sigma = 3,
                                            startparams = NULL,
                                            conv.check = c('ll', 'cor')),
                             verbose=T,
                             asympSE=F) {
  thecall <- match.call()

  tY <- as.matrix(dtm)
  D <- nrow(tY)  # number of documents
  V <- ncol(tY)  # number of words

  # possible reasons for aborting the algorithm
  stopifnot("Select integer for the dimension K." = is.null(K)==F)
  stopifnot("The number of selected dimensions K exceeds the number of fed-in
            documents D. It is suggested to keep K significantly smaller than D
            to maintain a well-posed model." = K <= D)
  stopifnot("Direction improperly selected." = dir[1]<=D)
  stopifnot("Direction improperly selected." = dir[2]<=D)
  stopifnot("Epsilon cannot be negative." = all(epsilon >= 0))

  # enforce control values
  tol <- control$tol %||% 1e-05
  sigma <- control$sigma %||% 3
  conv.check <- control$conv.check %||% 'll'

  optimfail <- "Warning: optim failed to converge"

  # enforce initial parameters
  if (is.null(control$startparams)) {
    theta_list <- list() # empty list
    beta_list  <- list() # empty list

    for (i in 1:K) {
      th <- rnorm(D, 0, 1) # fill with normal distributed initial values
      th <- (th-mean(th))/sd(th) # centralising the theta

      # Gram-Schmidt orthogonalisierung
      if (i > 1) {
        for (j in 1:(i - 1)) {
          proj <- sum(th * theta_list[[j]]) / sum(theta_list[[j]]^2)
          th <- th - proj * theta_list[[j]]
        }
        # centralisation
        th <- (th - mean(th)) / sd(th)
      }

      theta_list[[i]] <- th
      beta_list[[i]]  <- rnorm(V, 0, 1) # fill with normal distributed initial values
    }

    ## Name the lists appropriately
    names(theta_list) <- paste0("theta_", 1:K)
    names(beta_list)  <- paste0("beta_", 1:K)

    params <- list(
      theta = theta_list,
      beta = beta_list,
      psi = rnorm(V, 0, 1),
      alpha = rnorm(D, 0, 1)
    )
  } else {
    params <- control$startparams
    stopifnot((length(params$beta) == K) &&
                (length(params$psi) == ncol(tY)) &&
                (length(params$alpha) == nrow(tY)) &&
                (length(params$theta) == K))
  }

  ### Constraint
  constraint <- function(p, y, beta, psi) {

    constraints <- NULL
    theta <- matrix(p[1:(D * K)], nrow = D, ncol = K, byrow = F)

    # Mean constraints
    mean_vals <- colMeans(theta)
    constraints <- constraints <- c(epsilon - abs(mean_vals))

    # Variance constraints
    for (k in 1:K) {
      constraints <- c(constraints, epsilon - abs(var(theta[, k]) - 1))
    }

    if(K>1){
      # Covariance constraints (off-diagonal)
      Theta_Cov <- cov(theta)  # K × K covariance between dimensions

      # Extract off-diagonal elements
      for (i in 1:(K - 1)) {
        for (j in (i + 1):K) {
          constraints <- c(constraints, epsilon - abs(Theta_Cov[i, j]))
        }
      }}else{}

    return(constraints)
  }

  ### Log-likelihood function (overall LL)
  LL <- function(params, tY) {
    logexpected <- matrix(0, ncol = V, nrow = D)
    for (i in 1:K) {
      logexpected <- logexpected + params$theta[[i]] %*% t(params$beta[[i]])
    }

    logexpected <- sweep(logexpected, 1, params$alpha, "+")  # add alpha per row
    logexpected <- sweep(logexpected, 2, params$psi, "+")    # add psi per column

    # sum over all entries
    LL <- sum(sum(tY * logexpected - exp(logexpected))) -
      0.5 * sum(sapply(params$beta, function(b) sum(b^2))) / sigma^2
    return(LL)
  }

  ### Log-Likelihood for Psi and Beta
  LL.psi.beta <- function(p, y, theta, alpha) {

    # Reshape vector p into beta (V × K) and psi (V × 1)
    beta <- matrix(p[1:(V * K)], nrow = V, ncol = K, byrow = TRUE)
    psi <- p[(V * K + 1):(V * K + V)]

    # Compute eta: D × V matrix
    eta <- theta %*% t(beta)  # D × V
    eta <- sweep(eta, 1, alpha, "+")  # add alpha across rows
    eta <- sweep(eta, 2, psi, "+")   # add psi across columns

    # Compute log-likelihood
    ll <- sum(sum(y * eta - exp(eta))) - 0.5 * sum(beta^2) / sigma^2

    return(-ll)
  }

  ### Gradient for LL.psi.beta
  DLL.psi.beta <- function(p, y, theta, alpha) {
    # Reshape vector p into beta (V × K) and psi (V × 1)
    beta <- matrix(p[1:(V * K)], nrow = V, ncol = K, byrow = TRUE)
    psi <- p[(V * K + 1):(V * K + V)]

    # Compute eta: D × V matrix
    eta <- theta %*% t(beta)
    eta <- sweep(eta, 1, alpha, "+")
    eta <- sweep(eta, 2, psi, "+")

    # Compute mu (expected counts)
    mu <- exp(eta)

    # Compute gradients
    grad_beta <- t(y - mu) %*% theta - beta / sigma^2  # V × K
    grad_psi <- colSums(y - mu)  # V × 1

    # Flatten gradients into a long vector
    dll <- c(as.vector(grad_beta), grad_psi)

    return(-dll)
  }

  LL.alpha.theta <- function(p, y, beta, psi) {
    # Reshape p into (D × K) theta matrix and (D × 1) alpha vector
    theta <- matrix(p[1:(D * K)], nrow = D, ncol = K, byrow = TRUE)
    alpha <- c(0, p[(D * K + 2):length(p)])  # Ensure alpha[1] = 0

    # Compute eta: (D × V) matrix
    eta <- matrix(0, nrow = D, ncol = V)
    for (i in 1:K) {
      eta <- eta + theta[, i] %*% t(beta[[i]])
    }
    eta <- sweep(eta, 1, alpha, "+")
    eta <- sweep(eta, 2, psi, "+")

    # Compute log-likelihood
    ll <- sum(y * eta - exp(eta))
    return(-ll)  # Negative because we minimize
  }

  DLL.alpha.theta <- function(p, y, beta, psi) {
    # Reshape p into (D × K) theta matrix and (D × 1) alpha vector
    theta <- matrix(p[1:(D * K)], nrow = D, ncol = K, byrow = TRUE)
    alpha <- c(0, p[(D * K + 2):length(p)])  # Ensure alpha[1] = 0

    # Compute eta: (D × V) matrix
    eta <- matrix(0, nrow = D, ncol = V)
    for (i in 1:K) {
      eta <- eta + theta[, i] %*% t(beta[[i]])
    }
    eta <- sweep(eta, 1, alpha, "+")
    eta <- sweep(eta, 2, psi, "+")

    mu <- exp(eta)  # Expected counts

    # Compute gradient w.r.t theta (D × K)
    grad_theta <- matrix(0, nrow = D, ncol = K)
    for (i in 1:K) {
      grad_theta[, i] <- rowSums((y - mu) * matrix(rep(beta[[i]], each = D), nrow = D))
    }

    # Compute gradient w.r.t alpha (D × 1)
    grad_alpha <- rowSums(y - mu)
    grad_alpha[1] <- 0  # Ensure alpha[1] = 0

    grad_theta <- t(grad_theta)

    # Flatten gradients and combine
    dll <- c(as.vector(grad_theta), grad_alpha)  # Exclude alpha_1
    return(-dll)  # Negative because we minimize
  }

  if (conv.check=='ll'){
    ll <- LL(params, tY)
    if (verbose) cat("0\tLL:", ll, "\n")
  }

  iter <- 0
  diff <- Inf

  while (diff > tol) {
    if (conv.check == 'cor') {
      old.theta <- params$theta
    } else {
      old.ll <- LL(params, tY)
    }

    if (verbose) cat(iter, "\t")

    # constrained optimisation of theta and alpha
    resa <- constrOptim.nl(
      par = c(unlist(params$theta), params$alpha),
      fn = LL.alpha.theta,
      gr = DLL.alpha.theta,
      hin = constraint,
      heq = NULL,
      y = tY,
      beta = params$beta,
      psi = params$psi,
      control.outer = list(trace = 0, itmax = 100, eps = tol),
      control.optim = list(maxit = 500, factr = 1e7, pgtol = tol, trace = 0)
    )

    if (!is.null(resa$convergence) && resa$convergence != 0)
      cat("Warning: constrOptim.nl failed to converge on alpha/theta step\n")

    for (d in 1:K) {
      params$theta[[d]] <- resa$par[(1+(d-1)*D):(d*D)] # getting thetas of the K dimensions
    }
    params$alpha <- resa$par[(K*D+1):(length(resa$par))] # getting alpha

    # global identification
    for (i in 1:K) {
      if (params$theta[[i]][dir[1]] < params$theta[[i]][dir[2]])
        params$theta[[i]] <- params$theta[[i]]*(-1)
    }

    # Estimate beta and psi
    resa <- optim(par=c(unlist(params$beta), params$psi),
                  fn=LL.psi.beta,
                  gr=DLL.psi.beta,
                  y=tY,
                  theta=do.call(cbind, params$theta),
                  alpha=params$alpha,
                  method=c('BFGS'),
                  control = list(maxit = 500, factr = 1e4, pgtol = tol, trace = 0))

    for (d in 1:K) {
      params$beta[[d]] <- resa$par[(1+(d-1)*V):(d*V)] # getting betas of the K dimensions
    }
    params$psi <-  resa$par[(K*V+1):(length(resa$par))] # getting psi

    if (resa$convergence!=0)
      cat(optimfail, "while estimating betas and psi \n")

    # Convergence check
    if (conv.check == 'll') {
      ll <- LL(params, tY)
      diff <- abs((ll - old.ll) / old.ll)
      if (verbose) {
        cat("LL:", ll, "\n")
        flush.console()
      }
    } else {
      diff <- 1 - Reduce(`*`, sapply(1:K, function(d) cor(params$theta[[d]], old.theta[[d]])))
      if (verbose) {
        cat("1-cor:", diff, "\n")
        flush.console()
      }
    }
    iter <- iter + 1
  }

  model <- list(dir=dir,
                alpha=params$alpha,
                psi=params$psi,
                sigma=sigma,
                docs=docnames(dtm),
                features=featnames(dtm),
                ll=ll,
                sigma = sigma,
                epsilon = epsilon,
                K=K,
                data=t(tY),
                call = thecall)

  for (d in 1:K) {
    model[[paste0("theta_", d)]] <- params$theta[[d]]
    names( model[[paste0("theta_", d)]]) <- NULL
    model[[paste0("beta_", d)]] <- params$beta[[d]]
    names(model[[paste0("beta_", d)]]) <- NULL
  }

  ### Asymptotic standard errors
  if(asympSE==T){

    # generate standard errors
    for (d in 1:K) {
      model[[paste0("theta_SE_", d)]] <- params$theta[[d]]
    }

    ## Multinomial Negative Log-Likelihood for K-dimensional theta
    mnll <- function(tt, B, P, y) {
      linpred <- c(0, as.vector(tt %*% B + t(P)))  # K-dimensional linear predictor
      return(dmultinom(y, size=sum(y), prob=exp(linpred), log=TRUE))
    }

    ## Initialize model
    new.theta <- matrix(0, nrow=D, ncol=K)  # Document positions (D x K)
    T_all <- do.call(cbind, params$theta) # D x K
    B_all <- do.call(rbind, params$beta) # K x V
    P_all <- model$psi

    for (i in 1:D) {
      ## Construct beta and psi differences (excluding first word)
      new.B <- matrix(B_all[, 2:V] - B_all[, 1], nrow=K)  # K x (V-1)
      new.P <- matrix(P_all[2:V] - P_all[1], ncol = 1) # (V-1) x 1

      ## Optimize theta (K-dimensional)
      opt_result <- optim(par=T_all[i,], fn=mnll, B=new.B, P=new.P, y=tY[i,],
                          lower = -6, upper = 6,
                          method="L-BFGS-B", hessian=TRUE)

      ## Store theta estimates
      new.theta[i, ] <- opt_result$par

      ## Compute standard errors from Hessian
      neghess <- -opt_result$hessian
      invneghess <- solve(neghess)  # Inverse Hessian (K × K covariance matrix)
      SE <- sqrt(diag(invneghess))  # Extract diagonal (SE for each theta component)

      for (d in 1:K) {
        model[[paste0("theta_SE_", d)]][i] <- SE[d]
      }

    }

    for (d in 1:K) {
      names(model[[paste0("theta_SE_", d)]]) <- NULL
    }

  }else{}

  class(model) <- c("wordkrill", class(model))
  return(model)

}
