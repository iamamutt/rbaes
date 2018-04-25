#' Merge posterior predictions with data
#'
#' @param stanreg
#' @param newdata data to use for prediction
#' @param post_fun
#' @param var.sample
#' @param var.rows
#' @param var.yhat
#' @param ... options passed to rstanarm::posterior_predict
#'
#' @return [data.table::data.table]
#' @export
#'
#' @examples
#' #TODO:
merge_data_and_posterior <- function(stanreg, newdata = NULL,
                                     post_fun = rstanarm::posterior_predict,
                                     var.sample = ".sample", var.rows = ".obs",
                                     var.yhat = ".y", ...) {
  pred <- post_fun(stanreg, newdata = newdata, ...)
  pred <- data.table::as.data.table(pred)
  pred[, .__mcmc_sample := 1:.N]
  pred <- melt(pred,
    id.vars = ".__mcmc_sample",
    variable.name = ".__data_obs", value.name = ".__post_pred")
  pred[, .__data_obs := as.integer(.__data_obs)]
  pred[, .__mcmc_sample := as.integer(.__mcmc_sample)]

  if (is.null(newdata)) {
    newdata <- stanreg$data
  }

  if (data.table::is.data.table(newdata)) {
    newdata <- data.table::copy(newdata)
  } else {
    newdata <- data.table::as.data.table(newdata)
  }

  newdata[, .__data_obs := 1:.N]
  pred <- merge(newdata, pred, by = ".__data_obs", all.y = TRUE)
  default_vars <- c(".__mcmc_sample", ".__data_obs", ".__post_pred")
  rename_vars <- c(var.sample, var.rows, var.yhat)
  which_exist <- rename_vars %in% names(pred)

  if (any(which_exist)) {
    these_exist <- paste0(rename_vars[which_exist], collapse = ", ")
    warning(
      "New names: ", these_exist, " already exist. ",
      "Keeping default names.")
  }

  if (!all(which_exist)) {
    data.table::setnames(
      pred, default_vars[!which_exist],
      rename_vars[!which_exist])
  }

  pred
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
  rstanarm::stan_glmer(log_radon ~ floor + log_uranium + (1 | county),
    data = rstanarm::radon, iter = iter, chains = chains)
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
example_stan <- function(n = 30, seed = NULL, data_only = FALSE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  y1 <- rnorm(n, 75, 15)
  y2 <- rnorm(n, 100, 10)

  data <- list(y = c(y1, y2), x = rep(c(0, 1), each = n), n = n * 2)
  pars <- c("Beta", "Sigma", "log_lik")
  inits <- function() {
    list(Beta = c(75, 25), Sigma = 12.75)
  }
  model <- "
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
    data = data, pars = pars, inits = inits,
    model = model, fit = list())

  if (data_only) {
    return(stan_stuff)
  } else {
    if (!requireNamespace("rstan", quietly = TRUE)) {
      stop("cannot find package rstan")
    }

    stan_stuff$fit <- rstan::stan(
      model_code = model, model_name = "example",
      data = data, iter = 2100, chains = 3, pars = pars,
      init = inits, verbose = FALSE, open_progress = -1,
      control = list(adapt_delta = 0.8)
    )

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
  if ("stanreg" %in% class(object)) {
    n_post_draws(object)
  } else {
    if ("stanfit" %in% class(object)) {
      sum(object@sim$n_save - object@sim$warmup2)
    } else {
      stop("object is not some kind of stan fitted model")
    }
  }
}

set_col_order <- function(x, first_names) {
  data.table::setcolorder(x, c(first_names, names(x)[!names(x) %chin%
    first_names]))
  invisible(x)
}
