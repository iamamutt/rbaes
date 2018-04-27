#' Merge posterior predictions with data
#'
#' @param stanreg x
#' @param newdata data to use for prediction
#' @param post_fun x
#' @param var.sample x
#' @param var.rows x
#' @param var.yhat x
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
    newdata <- stanreg_dtbl(stanreg, model_frame = FALSE)
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

#' Scale a numeric vector
#'
#' Vectors are scaled to have a mean of 0.5 and standard deviation of 0.5.
#'
#' The function works similar to scale, in that attributes of the original mean and standard
#' deviation are saved in the object called \code{center} and \code{scale}. The adjustment is also
#' provided in the \code{adj} attribute.
#'
#' @param x numeric vector
#' @param y alternative numeric vector that contains the attributes to use to rescale x
#' @param na.rm remove NAs from the vector before finding stats
#'
#' @return scaled vector
#' @export
#'
#' @examples
#' x <- rnorm(1000, 50, 5)
#' hist(x, br=50)
#'
#' y <- unit_scale(x)
#' hist(y, br=50)
#' mean(y) # <- 0.5
#' sd(y) # <- 0.5
unit_scale <- function(x, y = NULL, na.rm = TRUE) {
  if (is.null(y)) {
    x_mean <- mean(x, na.rm = na.rm)
    s_std <- 2 * sd(x, na.rm = na.rm)
    x_adj <- 0.5
  } else {
    y_attr <- get_unit_attr(y)
    x_mean <- y_attr$m
    s_std <- y_attr$s
    x_adj <- y_attr$adj
  }

  z <- x_adj + ((x - x_mean) / s_std)

  attr(z, "center") <- x_mean
  attr(z, "scale") <- s_std
  attr(z, "adj") <- x_adj
  return(z)
}

#' Reverse unit scaling
#'
#' @param x object to be scaled back to original scale
#' @param y optional object that has center, scale, and adj attributes to use on x. If y is not
#' used, these attributes will be taken from x.
#'
#' @return numeric vector in original scale
#' @export
#'
#' @examples
#' z <- rnorm(10000, 100, 15)
#' y <- unit_scale(z)
#' rev_unit_scale(y)
rev_unit_scale <- function(x, y = NULL, na.rm = FALSE) {
  if (is.null(y)) {
    stats <- get_unit_attr(x)
  } else {
    stats <- get_unit_attr(y)
  }

  if (na.rm) {
    x <- na.omit(x)
  }

  with(stats, as.numeric((x - adj) * s + m))
}

#' Stan formatted Cholesky factored cov/cor matrix
#'
#' @param mat A cholesky factored covariance or correlation matrix
#'
#' @examples
#' # make matrix
#' L <- solve(rWishart(1, 6, diag(5))[,,1])
#'
#' # compare
#' l1 <- chol(L)
#' l2 <- stan_chol(L)
#' @export
stan_chol <- function(mat) {
  L <- chol(mat)
  l <- dim(L)
  Lp <- array(0, l)
  for (M in l[1]:1) {
    for (N in l[2]:1) {
      Lp[M, N] <- L[N, M]
    }
  }
  return(Lp)
}

#' Get model data from stanreg
#'
#' @param stanreg object of class stanreg
#' @param model_frame If `TRUE`, then extract from `$glmod$fr`, else use
#' `$data`.
#' @param get_y logical value. Include the output columns?
#'
#' @return data.table
#' @export
#' @examples
#' stanreg <- example_stanreg()
#' model_data <- stanreg_dtbl(stanreg, TRUE, TRUE)
stanreg_dtbl <- function(stanreg, model_frame = FALSE, get_y = FALSE) {
  if (model_frame) {
    dt <- copy(as.data.table(stanreg$glmod$fr))
    if (!get_y) {
      y <- sub("~", "", deparse(stanreg$formula[1:2]))
      dt[, eval(y) := NULL]
    }
  } else {
    dt <- stanreg$data
  }

  dt
}

#' @export
dot_dot_len <- function(...) {
  length(match.call(expand.dots = TRUE)) - 1L
}


