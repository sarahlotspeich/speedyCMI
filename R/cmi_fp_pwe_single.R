#' Single, fully parametric conditional mean imputation for a right-censored covariate (piecewise exponential distribution) with conditional means following Equation (11)
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using a piecewise exponential model to estimate the conditional survival function and then uses an analytic solution to compute conditional means, as in Equation (11) of the manuscript.
#'
#' @param imputation_formula imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_formula}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#' @param nintervals integer, number of disjoint subintervals used to split the hazard function of \code{W}. Must specify either this or \code{breaks}.
#' @param breaks vector, fixed subinterval boundaries used to split the hazard function of \code{W}. Must specify either this or \code{nintervals}.
#' @return
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @importFrom survival survreg
#' @importFrom survival Surv
#' @importFrom survival psurvreg
#' @importFrom survival survSplit
#' @importFrom formula.tools terms

cmi_fp_pwe_single = function(imputation_formula, data, maxiter = 100, nintervals = NULL, breaks = NULL) {
  ## Checks
  if (is.null(nintervals)) {
    if (is.null(breaks)) {
      stop('One of nintervals / breaks must be non-null')
    }
  }

  ## Perform checks and start setup
  setup = cmi_pwe_setup(imputation_formula, data)
  W  = data[[setup$Wname]]
  Delta = data[[setup$Deltaname]]

  ## Get terms from formula
  tt = terms(imputation_formula,
             data = data)

  ## Compute breaks based on nintervals quantiles (if specified)
  if (!(is.null(nintervals))) {
    if (is.null(breaks)) {
      breaks = rep(NA, nintervals + 1)
      breaks[1] = 0
      breaks[length(breaks)] = Inf
      ## Compute percentiles of observed event times
      Wobs = W[Delta == 1]
      Wobsqtl = quantile(x = Wobs,
                         probs = (1:(nintervals - 1)) / nintervals)
      breaks  = c(0, Wobsqtl, max(Wobs) + 100000)
    } else {
      warning('Both nintervals and breaks were specified. Ignoring nintervals and using breaks...')
    }
  }
  ## Create pseudo-data for Poisson likelihood;
  pdata = survSplit(imputation_formula,
                     data = data,
                     cut = breaks,
                     episode = "interval",
                     start = "start",
                     id = 'id')
  ## compute time-at-risk
  pdata$risktime = pdata[[setup$Wname]] - pdata[['start']]
  ## Obtain breaks assuming last cutpoint = infinity; create names for factor variable
  breaks[length(breaks)] = Inf
  intvlnames = paste0('[', breaks[-length(breaks)], ', ', breaks[-1],')')
  pdata$interval = factor(pdata$interval,
                          labels = intvlnames)

  ## Create glm formula--drop any uncessary variables
  tlabs = attr(tt, 'term.labels')
  fmla = reformulate( c('interval', tlabs, 'offset(log(risktime))'), response = setup$Deltaname, intercept = FALSE )

  ## Fit PWE model via Poisson likelihood
  fit = glm(formula = fmla,
            data = pdata,
            family = poisson)

  ## -----------------------
  ## Impute censored times
  ## -----------------------
  ## Obtain censored values
  cenindx = data[[setup$Deltaname]] == 0
  lower = data[[setup$Wname]][cenindx]
  idcen = (1:nrow(data))[cenindx]

  ## Filter pdata so that only censored individuals are included
  pdata = pdata[which(pdata$id %in% idcen), ]
  data2 = data
  data2$id = 1:nrow(data2)
  data2 = data2[which(data2$id %in% idcen), ]
  C = data2[[setup$Wname]] ## Censored times

  ## Extract design matrix (ignoring intercept terms / basline hazards)
  X = model.matrix(imputation_formula, data2)
  if ( '(Intercept)' %in% colnames(X) ) {
    intindx = which(colnames(X) == '(Intercept)')
    X = X[, -intindx]
    if (is.null(ncol(X))) {
      X = data.matrix(frame = X)
     }
  }
  ## Extract coefficients (ignoring intercept terms)
  beta = coef(fit)
  blhazindx = which(names(beta) %in% paste0('interval', intvlnames))
  logblhaz = beta[blhazindx]
  beta = beta[-blhazindx]
  lp = (X %*% beta)[, 1]
  J = length(breaks) - 1

  ## Compute survival function at each endpoint
  ## and at censored time
  cumHaz = matrix(0, nrow = nrow(X), ncol = J+1)   ## cumulative hazard at each cut point
  cumHazCen = matrix(0, nrow = nrow(X), ncol = J+1)   ## cumulative hazard at each point until censoring
  survProbs = matrix(1, nrow = nrow(X), ncol = J+1)   ## survival probabilities at each cut point
  survProbsCen = matrix(1, nrow = nrow(X), ncol = J+1)   ## survival probabilities at each cut point
  survInt = numeric(nrow(X))  ## \int_{0}^{\infty} S(t) dt
  survIntCen = numeric(nrow(X))  ## \int_{0}^{C} S(t) dt
  logSurvCen = numeric(nrow(X))  ## log S(C)
  for (j in 1:J) {
    ## Compute cumulative hazard and survival at interval j
    cumHaz[, j+1]    = cumHaz[, j] + exp( logblhaz[j] + lp + log(breaks[j+1] - breaks[j]) )
    survProbs[, j+1] = exp(-cumHaz[, j+1])
    survInt          = survInt + exp(log(survProbs[, j] - survProbs[, j+1]) - logblhaz[j])
    ## See if subject was at risk in interval j;
    ##  compute time at risk and cumulative hazard in that interval
    indx_j                    = (C > breaks[j])
    risktime_j                = pmin(C[indx_j], breaks[j+1]) - breaks[j]
    cumHazCen[indx_j, j+1]    = cumHazCen[indx_j, j] + exp( logblhaz[j] + lp[indx_j] + log(risktime_j) )
    survProbsCen[indx_j, j+1] = exp(-cumHazCen[indx_j, j+1])
    survIntCen[indx_j]        = survIntCen[indx_j] +
      exp(log(survProbsCen[indx_j, j] - survProbsCen[indx_j, j+1]) - logblhaz[j])
    logSurvCen[indx_j] = logSurvCen[indx_j] - exp( logblhaz[j] + lp[indx_j] + log(risktime_j) )
  }
  ## Impute censoring time
  imp = C + exp( log( survInt - survIntCen ) - logSurvCen )

  ## Return data set with imputed time
  data$imp  = data[[setup$Wname]]
  data$imp[data[[setup$Deltaname]] == 0] = imp

  # Return input dataset with appended column imp containing imputed values
  return_list = list(imputed_data = data,
                     code = !any(is.na(data$imp)))
  return(return_list)
}
