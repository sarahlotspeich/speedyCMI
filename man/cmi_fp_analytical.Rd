% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmi_fp_analytical.R
\name{cmi_fp_analytical}
\alias{cmi_fp_analytical}
\title{Fully parametric conditional mean imputation for a right-censored covariate with conditional means using the analytical solution}
\usage{
cmi_fp_analytical(
  imputation_model,
  dist,
  data,
  nintervals = NULL,
  breaks = NULL,
  maxiter = 100,
  boots = 0,
  seed = NULL
)
}
\arguments{
\item{imputation_model}{imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.}

\item{dist}{imputation model distribution passed through to \code{survreg}. See \code{survreg} documentation for more details.}

\item{data}{dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_model}.}

\item{nintervals}{(only if \code{dist = "PWE"}) integer, number of disjoint subintervals used to split the hazard function of \code{W}. Must specify either this or \code{breaks}.}

\item{breaks}{(only if \code{dist = "PWE"}) vector, fixed subinterval boundaries used to split the hazard function of \code{W}. Must specify either this or \code{nintervals}.}

\item{boots}{(optional) numeric, for multiple imputation supply the desired number of imputations (obtained via bootstrapping) to \code{boots}. Default is \code{0}, which is single imputation.}

\item{seed}{(optional) numeric, for multiple imputation set the random seed for the bootstrapping with \code{seed}. Default is \code{NULL}, which does not reset \code{seed}.}

\item{max_iter}{(optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{max_iter = 100}.}
}
\value{
A list (or, in the case of multiple imputation, a list of \code{boots} lists) containing:
\item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
\item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
\item{aic}{Akaike information criterion (AIC) from the \code{imputation_model} fit.}
\item{bic}{Bayesian information criterion (AIC) from the \code{imputation_model} fit.}
\item{coefficients}{Vector of coefficients from the \code{imputation_model} fit.}
\item{scale}{Scale from the \code{imputation_model} fit.}
}
\description{
Fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then uses an analytical solution to compute conditional means.
}
