% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmi_fp_stabilized.R
\name{cmi_fp_stabilized}
\alias{cmi_fp_stabilized}
\title{Fully parametric imputation for a right-censored covariate with conditional means using the stabilized integral}
\usage{
cmi_fp_stabilized(
  imputation_model,
  dist,
  data,
  with_mean = FALSE,
  use_cumulative_hazard = TRUE,
  maxiter = 100,
  boots = 0,
  seed = NULL
)
}
\arguments{
\item{imputation_model}{imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.}

\item{dist}{imputation model distribution passed through to \code{survreg}. Options of \code{dist =} \code{"exponential"}, \code{"loglogistic"}, \code{"lognormal"}, and \code{"weibull"} currently accepted.}

\item{data}{Dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_model}.}

\item{with_mean}{logical, if \code{TRUE} the stabilized integral with the mean is used. Default is \code{FALSE}.}

\item{use_cumulative_hazard}{(optional) logical, if \code{use_cumulative_hazard = TRUE} and \code{with_mean = FALSE} the survival function is transformed to the cumulative hazard before integration. Default is \code{TRUE}.}

\item{boots}{(optional) numeric, for multiple imputation supply the desired number of imputations (obtained via bootstrapping) to \code{boots}. Default is \code{0}, which is single imputation.}

\item{seed}{(optional) numeric, for multiple imputation set the random seed for the bootstrapping with \code{seed}. Default is \code{NULL}, which does not reset \code{seed}.}

\item{max_iter}{(optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.}
}
\value{
A list (or, in the case of multiple imputation, a list of lists) containing:
\item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
\item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
\item{aic}{Akaike information criterion (AIC) from the \code{imputation_model} fit.}
\item{bic}{Bayesian information criterion (AIC) from the \code{imputation_model} fit.}
\item{coefficients}{Vector of coefficients from the \code{imputation_model} fit.}
\item{scale}{Scale from the \code{imputation_model} fit.}
}
\description{
Fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then computes conditional means using the stabilized integral.
}
