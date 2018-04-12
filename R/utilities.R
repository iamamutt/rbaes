#' Posterior predict data.table
#'
#' @param object stanfit object
#' @param newdata data to use for prediction
#' @param append add data columns to predictions
#' @param fn  posterior_predict or posterior_linpred
#' @param ... options passed to rstanarm::posterior_predict
#'
#' @return Posterior prediction matrix as a data.table, merging newdata if available.
#' @export
#'
#' @examples
#' #TODO:
post_pred_dtable <-
  function(object,
           newdata = NULL,
           append = TRUE,
           fn = rstanarm::posterior_predict,
           ...) {
    pred <- fn(object, newdata = newdata, ...)
    pred <- data.table::as.data.table(pred)
    pred[, .sample := 1:.N]
    pred <-
      melt.data.table(
        pred,
        id.vars = ".sample",
        variable.name = ".obs",
        value.name = ".y"
      )
    pred[, .obs := as.integer(.obs)]
    pred[, .sample := as.integer(.sample)]
    if (!is.null(newdata) & append) {
      if (!is.data.table(newdata)) {
        dt <- as.data.table(newdata)
      } else {
        dt <- data.table::copy(newdata)
      }
      dt[, .obs := 1:.N]
      pred <- merge(dt, pred, by = ".obs")
    }
    return(pred)
  }

#' Stan example model
#'
#' @param iter number of iterations for stan_glmer
#' @param chains number of chains for stan_glmer
#'
#' @return stanfit obj
#' @export
#'
#' @examples
#' stanfit <- example_stanreg()
example_stanreg <- function(iter = 1000, chains = 2) {
  rstanarm::stan_glmer(
    log_radon ~ floor + log_uranium + (1 | county),
    data = rstanarm::radon,
    iter = iter,
    chains = chains
  )
}

#' Example stan model and data
#'
#' @param n sample size per group
#' @param seed seed value to get same data samples
#' @param data_only return only the data without fitting the model
#'
#' @return list(data, pars, inits, model, fit)
#' @export
#'
#' @examples
#' stan_ex <- example_stan(seed=19851905)
example_stan <- function(n=30, seed=NULL, data_only=FALSE) {
  if (!is.null(seed))
    set.seed(seed)

  y1 <- rnorm(n, 75, 15)
  y2 <- rnorm(n, 100, 10)

  data <- list(y = c(y1, y2), x = rep(c(0,1), each=n), n = n*2)
  pars <- c("Beta", "Sigma", "log_lik")
  inits <- function()
    list(Beta = c(75, 25), Sigma = 12.75)
  model <-
    "
    data {
      int<lower=1> n;
      vector[n]    y;
      vector[n]    x;
    }
    parameters {
      vector[2]     Beta;
      real<lower=0> Sigma;
    }
    model {
      Sigma ~ student_t(7, 0, 10);
      y ~ normal(Beta[1] + Beta[2] * x, Sigma);
    }
    generated quantities {
      vector[n] log_lik;
      for (i in 1:n)
        log_lik[i] = normal_lpdf(y[i] | Beta[1] + Beta[2] * x[i], Sigma);
    }
    "

  stan_stuff <- list(
    data = data,
    pars = pars,
    inits = inits,
    model = model,
    fit = list()
  )

  if (data_only) {
    return(stan_stuff)
  } else {
    if (!requireNamespace('rstan', quietly = TRUE)) {
      stop('cannot find package rstan')
    }

    stan_stuff$fit <- rstan::stan(
      model_code = model,
      model_name = "example",
      data = data,
      iter = 2100,
      chains = 3,
      pars = pars,
      init = inits,
      verbose = FALSE,
      open_progress = -1,
      control = list(adapt_delta=0.8))

    return(stan_stuff)
  }
}

#' Alternative rstanarm print
#'
#' @param stanreg a stanreg object
#'
#' @param header header name
#'
#' @export
print_rstanarm <- function(stanreg, header = "Stan model summary") {
  txt <- capture.output(print(stanreg, digits = 3))
  est <- which(txt == "Estimates:"):length(txt)
  cat("\n")
  mejr::print_sec(header)
  cat(stringr::str_c(txt[est], collapse = "\n"))
  cat("\n\n")
  print(rstanarm::prior_summary(stanreg))
  cat("\n")
  return(invisible(NULL))
}

#' Number of MCMC samples from model
#'
#' @param object stanreg or stanfit object
#'
#' @return scalar value
#' @export
#'
#' @examples
#'
n_mcmc_samples <- function(object) {
  if ('stanreg' %in% class(object)) {
    n_post_draws(object)
  } else if ('stanfit' %in% class(object)) {
    sum(object@sim$n_save - object@sim$warmup2)
  } else {
    stop('object is not some kind of stan fitted model')
  }
}
