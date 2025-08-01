% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wordkrill.R
\name{wordkrill}
\alias{wordkrill}
\title{K-dimensional Wordkrill Model with Conditional Maximum Likelihood}
\usage{
wordkrill(
  dtm,
  K = NULL,
  dir = c(1, 2),
  epsilon = 1,
  control = list(tol = 1e-08, sigma = 3, startparams = NULL, conv.check = c("ll", "cor")),
  verbose = T,
  asympSE = F
)
}
\arguments{
\item{dtm}{A quanteda DTM object}

\item{K}{integer to select the dimensionality of the Wordkrill model}

\item{dir}{set global identification by forcing \code{theta_k[dir[1]]} <
\code{theta_k[dir[2]]}}

\item{epsilon}{numeric to  determine the strictness of the epsilon
environment}

\item{control}{list of estimation options (tol: tolerance for likelihood
procedure; sigma: regularization parameter for betas in poisson form;
startparams: initial values of the parameters to be estimated; conv.check:
selecting the conversion check)}

\item{verbose}{logical, produce a running commentary}

\item{asympSE}{logical, compute asymptotcia standard errors}
}
\value{
An object of class wordkrill. This is a list containing:
\item{dir}{global identification of the dimension}
\item{theta_1 ... theta_K}{document positions dim 1 to K}
\item{alpha}{document fixed effects}
\item{beta_1 ... beta_K}{word slope parameters dim 1 to K}
\item{psi}{word fixed effects}
\item{docs}{names of the documents}
\item{features}{names of features}
\item{sigma}{regularization parameter for betas in Poisson form}
\item{ll}{final log likelihood}
\item{theta_SE_1 ... theta_SE_K}{standard errors for document position dim 1
to K}
\item{data}{the original data}
}
\description{
This function estimates a Wordkrill aka a K-dimensional Wordfish model,
based on a Document-Term Matrix (DTM). It outputs document positions
(\code{theta_1},...,\code{theta_K}) and word parameters (\code{beta_1},..., \code{beta_K}),
besides the usual wordfish parameters estimated with conditional maximum
likelihood.
}
\examples{

simulation <- sim_wordkrill(docs = 10, vocab =20, K=2, doclen =5000,
                            dist = "normal")

test_model <- wordkrill(dtm=simulation$Y, dir=c(1,2), K=2,
                               epsilon = 0.1,
                               control=list(tol=1e-08,
                                            sigma=3,
                                            startparams=NULL,
                                            conv.check='ll'),
                               verbose=TRUE, asympSE=TRUE)

print(summary(object = test_model, top_n = 15))

}
\references{
Slapin and Proksch (2008) 'A Scaling Model for Estimating
Time-Series Party Positions from Texts.' American Journal of Political
Science 52(3): 705-772.

Lowe (2020) 'Austin: Do Things with Words'.
https://conjugateprior.github.io/austin

Riesch (2025) 'Wordkrill: Extending Wordfish into the
multidimensional political space'. https://doi.org/10.48550/arXiv.2506.20275
}
\author{
Benjamin Riesch
}
