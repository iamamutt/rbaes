
#' Stan example model
#'
#' @return stanfit obj
#' @export
#'
#' @examples
#' stanfit <- example_stanreg()
example_stanreg <- function() {
  rstanarm::stan_glmer(
    log_radon ~ floor * log_uranium + (1 | county),
    data = rstanarm::radon,
    iter = 1000, chains = 2)
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


n_post_draws <- function(stanfit) {
  pss <- stanfit$stanfit@sim$n_save
  sum(pss - stanfit$stanfit@sim$warmup2)
}


ppi <- function(data, fit, fun = identity, form = NULL, prob = 0.9, n = 1)
{
  pp <- do.call(rbind, lapply(1:n, function(i) {
    rstanarm::posterior_predict(fit, data, fun = fun, re.form = form)
  }))
  pp_int <- rstantools::posterior_interval(pp, prob = prob)
  yhat <- data.table(apply(pp, 2, median), pp_int)
  names(yhat) <- c("c", "l", "r")
  return(cbind(data, yhat))
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

cdata_list <- function() {
  library(data.table)
  stanreg <- example_stanreg()
  cdata <- list(
    c1 = contrast_data(stanreg, TRUE, margin_ignore = 'floor', subset_expression = .~ log_uranium < -0.5),
    c2 = contrast_data(stanreg, TRUE, margin_ignore = 'floor', subset_expression = .~ log_uranium >= -0.5 & log_uranium <= 0.5),
    c3 = contrast_data(stanreg, TRUE, margin_ignore = 'floor', subset_expression = .~ log_uranium > 0.5))

  assign('stanreg', stanreg, .GlobalEnv)
  assign('cdata', cdata, .GlobalEnv)
  return(NULL)
}


