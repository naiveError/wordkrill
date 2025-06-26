#' Summarizes a Wordkrill model
#' 
#' @param object a fitted Wordkrill model
#' @param ... additional arguments (e.g., top_n to select number of displayed
#' features)
#' 
#' @method summary wordkrill
#' @export
summary.wordkrill <- function(object, ...) {
  args <- list(...)
  top_n <- top_n <- if (!is.null(args$top_n)) args$top_n else 15
  
  m <- object
  K <- m$K
  docs <- m$docs
  terms <- m$features
  alpha <- m$alpha
  psi <- m$psi
  
  # Assemble document positions
  theta_df <- data.frame(Doc = docs)
  for (k in 1:K) {
    theta_k <- m[[paste0("theta_", k)]]
    se_k <- if (!is.null(m[[paste0("theta_SE_", k)]])) m[[paste0("theta_SE_", k)]] else rep(NA, length(theta_k))
    theta_df[[paste0("theta_", k)]] <- theta_k
    theta_df[[paste0("se_", k)]] <- se_k
  }
  
  # Prepare beta and psi for each dimension
  feature_summary <- list()
  for (k in 1:K) {
    beta_k <- m[[paste0("beta_", k)]]
    top_indices <- order(abs(beta_k), decreasing = TRUE)[1:min(top_n, length(beta_k))]
    feature_summary[[paste0("dim_", k)]] <- data.frame(
      term = terms[top_indices],
      beta = beta_k[top_indices],
      psi = psi[top_indices],
      row.names = NULL
    )
  }
  
  result <- list(
    call = m$call,
    theta = theta_df,
    features = feature_summary,
    K = K
  )
  class(result) <- "summary.wordkrill"
  return(result)
}

#' Print summarised Wordkrill model
#' 
#' @param x a summary.wordkrill object
#' @param ... additional arguments
#'
#' @method print summary.wordkrill
#' @export
print.summary.wordkrill <- function(x, ...) {
  cat("##\n## Estimated Document Positions:\n")
  print(x$theta, row.names = FALSE, digits = 4)
  
  for (k in 1:x$K) {
    cat("\n## Most Influential Features (Dimension", k, "):\n")
    print(x$features[[paste0("dim_", k)]], row.names = FALSE, digits = 4)
  }
}


#' Get document positions
#' @param model a fitted Wordkrill model 
#' @export
get_document_para <- function(model) {
  # Create summary to extract theta + SE
  summary_model <- summary.wordkrill(model) # summarises wordkrill model
  theta_df <- summary_model$theta
  results <- cbind(theta_df, alpha = model$alpha) # add alpha
  return(as.data.frame(results))
}

#' Simulate data and parameters for a Wordkrill model
#'
#' Simulates data assuming K-dimensional Wordkrill model. Word counts are drawn
#' from Poisson distributions with log means determined by multiple document
#' positions and word slopes (beta), plus word intercepts (psi).
#'
#' @param docs Number of documents
#' @param vocab Number of words (must be divisible by 4)
#' @param doclen Document lengths (scalar or vector of length `docs`)
#' @param K Number of latent dimensions
#' @param dist Distribution for document positions (per dimension)
#' @param scaled Whether to standardize document positions
#' @return List with elements:
#' \item{Y}{Document-term matrix}
#' \item{theta}{Document positions (docs × K)}
#' \item{doclen}{Document lengths}
#' \item{psi}{Word intercepts}
#' \item{beta}{Word slopes (vocab × K)}
#' 
#' @importFrom stats rmultinom
#' @importFrom quanteda as.dfm
#' @export
sim_wordkrill <- function(docs = 10,
                          vocab = 20,
                          doclen = 500,
                          K = 2,
                          dist = c("spaced", "normal"),
                          scaled = TRUE) {
  
  if ((vocab %% 4) > 0)
    stop("This function assumes vocab is divisible by 4")
  
  dist <- match.arg(dist)
  
  # Generate orthogonal and standardized document positions (theta)
  raw_theta <- matrix(rnorm(docs * K), nrow = docs, ncol = K)
  
  # Orthogonalize using QR decomposition
  qr_decomp <- qr(raw_theta)
  theta_ortho <- qr.Q(qr_decomp)[, 1:K]  # Take first K orthogonal columns
  
  # Scale each column to mean 0, sd 1
  theta <- scale(theta_ortho)
  
  theta <- theta[order(theta[, 1]), , drop = FALSE]  # sort by 1st dim
  
  # Generate word intercepts (psi) and word slopes (beta): vocab × K
  psi <- rep(0, vocab)  # intercepts
  beta <- matrix(NA, nrow = vocab, ncol = K)
  
  for (k in 1:K) {
    beta[, k] <- rep(c(-1, 1), each = vocab / 2)
  }
  
  if (length(doclen) == 1)
    doclen <- rep(doclen, docs)
  
  # Compute log-mu: vocab × docs
  logmu <- matrix(0, nrow = vocab, ncol = docs)
  for (d in 1:docs) {
    for (k in 1:K) {
      logmu[, d] <- logmu[, d] + beta[, k] * theta[d, k]
    }
    logmu[, d] <- logmu[, d] + psi
  }
  
  mu <- exp(logmu)
  
  # Sample from multinomial
  data <- matrix(0, nrow = vocab, ncol = docs)
  for (d in 1:docs) {
    data[, d] <- rmultinom(1, doclen[d], mu[, d])
  }
  
  zfill <- function(vals) {
    v <- as.character(vals)
    wid <- max(nchar(v))
    formatstr <- paste0("%0", wid, "d")
    sprintf(formatstr, vals)
  }
  
  rownames(data) <- paste0("W", zfill(1:vocab))
  colnames(data) <- paste0("D", zfill(1:docs))
  
  val <- list(
    Y = as.dfm(data),
    theta = theta,
    doclen = doclen,
    psi = psi,
    beta = beta
  )
  
  class(val) <- c('wordkrill.simdata', class(val))
  
  return(val)
}


