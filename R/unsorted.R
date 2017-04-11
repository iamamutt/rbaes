
#' Stan test data
#'
#' @param n number of observations
#'
#' @return list of stan options
#' @examples
#' fake_data <- stan_test_data(n=30)
stan_test_data <- function(n=100){
  y <- sort(rnorm(n, 100, 15))
  x <- sort(sample(c(0, 1), n, replace = TRUE))
  data <- list(y = y, x = x, n = n)
  pars <- c("Beta", "Sigma", "log_lik")
  inits <- function()
    list(Beta = c(88, 24), Sigma = 10)
  model <- "
    data {
    int<lower=1> n;
    vector[n] y;
    vector[n] x;
    }
    parameters {
    vector[2]     Beta;
    real<lower=0> Sigma;
    }
    model {
    Sigma ~ cauchy(0, 3);
    y ~ normal(Beta[1] + Beta[2] * x, Sigma);
    }
    generated quantities {
    vector[n] log_lik;
    for (i in 1:n)
    log_lik[i] = normal_lpdf(y[i] | Beta[1] + Beta[2] * x[i], Sigma);
    }
    "
  return(list(
    data = data,
    pars = pars,
    inits = inits,
    model = model
  ))
}

#' Stan Mejr
#'
#' Wrapper function for simplifying \code{stan}. See \code{?rstan::stan} for more info.
#'
#' Also provides additional plots and repackages data into a single object
#'
#' @param model model code or file name, both as character objects
#' @param name name of model being fitted
#' @param data data in stan format
#' @param pars named parameters. Defaults to kept parameters
#' @param samples list of sample parameters. See example for list structure.
#' @param init initial values
#' @param control control paramters for algorithm
#' @param out directory for where to save output files
#' @param parallel use multiple cores for running multiple chains
#' @param debug do a test run with options saved to global env
#' @param repackage add additional objects to the output of this function
#' @param ... other arguments from the main \code{rstan::stan} function
#'
#' @return list of objects: \code{stan_mcmc, pram, central, env, repack}
#' @examples
#' # Using fake data for example
#' stan_dat <- stan_test_data()
#'
#' # fit model
#' stan_fit <- stan_mejr(
#'      model=stan_dat$model,
#'      name="example_model",
#'      data=stan_dat$data,
#'      pars=stan_dat$pars,
#'      samples=list(n_chains = 4, n_final = 1200, n_thin = 2, n_warm = 800),
#'      init=stan_dat$init,
#'      out="~/Desktop",
#'      parallel=FALSE)
#'
#' # Run default example model
#' stan_fit <- rstan_debug <- stan_mejr(debug=TRUE)
stan_mejr <- function(
  model,
  name,
  data,
  pars,
  samples,
  init,
  out = NULL,
  parallel = FALSE,
  debug = FALSE,
  repackage = NULL,
  ...)
{
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Need to install rstan first")
  }

  ## options list ------------------------------------------------------------
  stan_opts <- list()

  # check model name
  if (missing(name))
    name <- "mejr_model"

  stan_opts[["model_name"]] <- name

  # check for samples argument
  if (missing(samples)) {
    samples <-
      list(
        n_chains = 4,
        n_final = 1200,
        n_thin = 2,
        n_warm = 800
      )
  }

  if (length(samples) < 4) {
    stop(simpleError(
      paste0(
        "Make sure to specify the following options: ",
        "n_chains, n_final, n_thin, n_warm"
      )
    ))
  }

  stan_opts[["chains"]] <- samples$n_chains
  stan_opts[["thin"]] <- samples$n_thin
  stan_opts[["warmup"]] <- samples$n_warm
  stan_opts[["iter"]] <-
    ceiling(samples$n_final * samples$n_thin /
              samples$n_chains + samples$n_warm)

  # no need for parallel if only one chain or debug
  if (parallel & (debug || samples$n_chains == 1)) {
    message("1 chain found: parallel was set to FALSE")
    parallel <- FALSE
  }

  if (parallel) {
    rstan::rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
  } else {
    rstan::rstan_options(auto_write = FALSE)
    options(mc.cores = 1)
  }

  # check for missing initializations argument
  if (missing(init)) {
    init <- "random"
  }
  stan_opts[["init"]] <- init
  stan_opts[["enable_random_init"]] <- TRUE

  # Check for missing data argument
  if (missing(data)) {
    data <- stan_test_data()$data
    warning(simpleWarning("No data suppled. Used example data."))
  }
  stan_opts[["data"]] <- data

  # check for missing model parameters
  if (missing(pars)) {
    pars <- NA
  }
  stan_opts[["pars"]] <- pars

  # check for missing model file/obj
  if (missing(model)) {
    warning(simpleWarning("No model supplied. Used example file."))
    model <- stan_test_data()$model
  }

  if (class(model) == "stanmodel") {
    stan_opts[["object"]] <- model
    stan_opts[["model_name"]] <- NULL
    stan_method <- rstan::sampling
  } else {
    stan_method <- rstan::stan
    nlines_model <- length(readLines(textConnection(model)))
    if (nlines_model > 1) {
      stan_opts[["model_code"]] <- model
    } else {
      stan_opts[["file"]] <- model
    }
  }

  # combine arguments
  op <- list(...)
  toset <- !(names(stan_opts) %in% names(op))
  stan_opts <- c(stan_opts[toset], op)

  # set debug parameters and save opts to global
  if (debug) {
    if (is.list(stan_opts$init)) {
      stan_opts$init <- list(stan_opts$init[1][[1]])
    }
    stan_opts$chains = 1
    stan_opts$thin = 1
    stan_opts$warmup = 25
    stan_opts$iter = 50
    args <- c("stan_opts")
    for (i in args) {
      assign(bquote(.(i)), get(i, environment()), .GlobalEnv)
    }
  }

  ## start -------------------------------------------------------------------

  message("\n\nModel compiling and sampling initiated")
  startDate <- date()
  stan_fitted <- do.call(stan_method, stan_opts)
  message("\n\nModel sampling finished...")
  options(mc.cores = parallel::detectCores())

  ## Output list -------------------------------------------------------------

  message("\n\nPacking output")
  central <- tryCatch(
    stan_point_est(stan_fitted),
    error = function(cond) {
      message("Point estimation caused an error")
      message(cond)
      return(NULL)
    },
    warning = function(cond) {
      message("Point estimation caused a warning")
      message(cond)
      return(NULL)
    },
    finally = {
      NULL
    }
  )

  rstan_pack <- list(
    stan_mcmc = stan_fitted,
    central = central,
    env = stan_opts,
    repack = repackage
  )

  ## print results -----------------------------------------------------------

  old_opts <- options()[c("width", "max.print")]
  print_results <- function(go = NULL) {
    if (!is.null(go)) {
      message("\n\nSaving contents to directory")
      options(list(max.print = 1e8, width = 1000))

      fout <- function(filename, ...) {
        file.path(out, filename, ...)
      }

      sink(file = fout(paste0("results-", name, ".txt")),
           type = "output")
      print_sec("Runtime")
      print(startDate)
      cat("\n")
      cat("\n")
      print(date())
      cat("\n\n")
      print(rstan::get_elapsed_time(stan_fitted))
      print_sec(paste("Stan Model:", name))
      print(stan_fitted,
            digits = 4,
            probs = c(0.025, 0.975))
      sink()
      save(rstan_pack,
           file = fout(paste0("saved_stan-", name, ".Rdata")))
    }
    return(invisible(NULL))
  }

  tryCatch(
    print_results(out),
    error = function(cond) {
      message("Printing results caused an error")
      message(cond)
      return(NULL)
    },
    warning = function(cond) {
      message("Printing results caused a warning")
      message(cond)
      return(NULL)
    },
    finally = {
      options(old_opts)
    }
  )

  ## return object -----------------------------------------------------------
  return(rstan_pack)
}


#' Chain convergence stats
#'
#' This is to see which chains should be removed, if any.
#'
#' If you get an error while plotting and using RStudio, try to make plot window bigger.
#'
#' @param stanmodel A fitted stan model
#' @param pars Names of parameters to calculate convergence statistics, defaults to all parameters, including lp__
#' @param view Plot each chain's log-prob and its distribution
#' @examples
#' rstan_pack <- stan_mejr(parallel=TRUE)
#' chain_convergence(rstan_pack$stan_mcmc)
stan_chain_convergence <- function(stanmodel, pars, view=FALSE) {
  requireNamespace("rstan", quietly = TRUE)

  if (missing(pars)) {
    pars <- stanmodel@sim$pars_oi
    pars <- pars[!pars %in% c("log_lik")]
  }

  x <- rstan::extract(stanmodel, pars=pars, permuted=FALSE, inc_warmup=TRUE)
  w <- stanmodel@sim$warmup2
  d <- dim(x)
  l <- d[2]
  dn <- dimnames(x)
  lp <- rstan::extract(stanmodel, pars="lp__", permuted=FALSE, inc_warmup=FALSE)
  dat <- list()

  nc <- ceiling(sqrt(l))
  nr <- ceiling((l)/nc)
  par(mfrow = c(nc, nr))

  for (i in 1:l) {
    y <- array(x[, i, ], c(d[1], 1, d[3]))
    dimnames(y) <- list(iterations=NULL,
                        chains=dn$chains[i],
                        parameters=dn$parameters)
    s <- rstan::monitor(y, warmup = w[i], print=FALSE)
    R_hat <- s[, "Rhat"]
    nonNaN <- !is.nan(R_hat)
    r <- mean(R_hat[nonNaN], na.omit=TRUE)
    n <- mean(s[nonNaN, "n_eff"])
    se <- mean(s[nonNaN, "se_mean"])
    lpv <- lp[,i,1]
    dat[[i]] <- data.frame(chain=i,
                           pars_avg_se=se,
                           pars_n_eff=n,
                           pars_rhat=r,
                           lp_var=var(lpv))

    if (view) {
      plot(lpv, type="l", main=paste("LP: chain", i))
      plot(density(lpv, adjust=0.75), main=paste("LP: chain", i))
    }
  }
  return(do.call(rbind, dat))
}

#' STAN point estimates
#'
#' @param stan_obj Stan fitted object
#' @param ... additional options passed to function \code{hdiq}
#'
#' @return list object with same parameter names
stan_point_est <- function(stan_obj, ...) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Need to install rstan first")
  }
  p <- rstan::extract(stan_obj, permuted=TRUE)
  pnames <- names(p)
  central <- list()


  for (i in seq_len(length(p))) {
    # i <- 1
    # temp_pram <- array(rnorm(100), c(5,5,4))

    temp_pram <- p[[i]]
    d <- dim(temp_pram)
    dl <- length(d)

    if (dl == 1) {
      # 1d
      y <- hdiq(temp_pram, ..., warn = FALSE)$mid
    } else if (dl == 2) {
      # 2d
      yl <- lapply(1:d[2], function(ii) {
        # ii=2
        tp2 <- temp_pram[, ii]
        mp <- hdiq(tp2, ..., warn = FALSE)$mid
        return(mp)
      })
      y <- c(yl, recursive = TRUE)
    } else if (dl == 3) {
      # 3d
      y <- array(0.0, c(d[2], d[3]))
      for (m in 1:d[2]) {
        for (n in 1:d[3]) {
          tp3 <- temp_pram[, m, n]
          mp <- hdiq(tp3, ..., warn = FALSE)$mid
          y[m, n] <- mp
        }
      }
    } else if (dl == 4) {
      y <- array(0.0, c(d[2:4]))
      for (k in 1:d[2]) {
        for (m in 1:d[3]) {
          for (n in 1:d[4]) {
            tp4 <- temp_pram[, k, m, n]
            mp <- hdiq(tp4, ..., warn = FALSE)$mid
            y[k,m,n] <- mp
          }
        }
      }
    } else
      y <- NA

    central[[pnames[i]]] <- y
  }
  return(central)
}


#' Drop unwanted chains from stanfit object
#'
#' @param object object of class 'stanfit'
#' @param chain_ids an integer vector chain numbers to drop
#'
#' @return new instance of stanfit class
#'
#' @examples
#' stan_drop_chains(myFittedObject, c(2,4))
stan_drop_chains <- function(object, chain_ids) {
  if (!is(object, "stanfit"))
    stop("object must be a fitted rstan object")

  if (missing(chain_ids))
    stop("please select chain numbers to drop for argument 'chain_ids'")

  ids <- unlist(lapply(object@stan_args, function(l) l$chain_id))

  if (length(ids) == 1)
    return(object)

  drop <- ids %in% chain_ids

  if (!any(drop))
    stop("chain_ids not found in stanfit object")

  for (chain in seq_along(drop)) {
    if (drop[chain]) {
      object@sim$samples[ids[chain]] <- NULL
      object@inits[ids[chain]] <- NULL
      object@sim$permutation[ids[chain]] <- NULL
      object@stan_args[ids[chain]] <- NULL
    }
  }

  object@sim$chains <- length(which(!drop))
  object@sim$n_save <- object@sim$n_save[!drop]
  object@sim$warmup2 <- object@sim$warmup2[!drop]

  for (chain in seq_len(object@sim$chains)) {
    object@stan_args[[chain]]$chain_id <- chain
  }

  new_obj <- new(
    "stanfit",
    model_name = object@model_name,
    model_pars = object@model_pars,
    par_dims = object@par_dims,
    mode = 0L,
    sim = object@sim,
    inits = object@inits,
    stan_args = object@stan_args,
    stanmodel = object@stanmodel,
    date = date(),
    .MISC = new.env(parent = emptyenv())
  )

  return(new_obj)
}


multi_merge <- function(data_list, setkeys = FALSE, ...) {
  Reduce(function(x, y) {
    if (setkeys) {
      if (data.table::is.data.table(x)) data.table::setkey(x)
      if (data.table::is.data.table(y)) data.table::setkey(y)
    }
    merge(x, y, ...)
  }, data_list)
}


# effects_list <- list(
#     main = list(
#         X = array(runif(25), c(5, 5)),
#         B = array(runif(100), c(20, 5))
#     )
# )
stan_yhat_i <- function(fx_list) {
  batch_items <- names(fx_list)

  # checks on list names, get X
  if (all(c("X", "B") %in% batch_items)) {
    X <- fx_list$X
  } else if (all(c("B", "formula", "data") %in% batch_items)) {
    X <- model.matrix(fx_list$formula, data = fx_list$data)
  } else {
    stop("Need X & B or B, formula & data in effects_list")
  }

  # get B
  B <- fx_list$B

  # check conformity
  if (ncol(X) != ncol(B)) {
    stop("columns in X must equal columns in B")
  }

  # make index table for data.table joins
  I <- data.table::data.table(`__.samp key` = 1:nrow(B), key = '__.samp key')

  # do calculations
  Y <- I[I, X %*% B[.I, ], by = .EACHI]

  # set names, data index, return Y
  data.table::setnames(Y, "V1", fx_list$valname)
  Y[,  '__.I key' := 1:.N, by = '__.samp key']
  data.table::setkey(Y, '__.I key')
  return(Y)
}

#' Predict values from X matrix and posterior samples
#'
#' @param effects_list A list of named lists that each contain an X matrix and the parameters values B
#' @param link Logical value deciding whether to link each batch in a single data.table or keep as a list
#' @param link_fun The linking function after summing predicted values. Defaults to the identity function
#'
#' @return A list or data.table
#' @examples
#' # Two batches of effects that will be summed, and
#' # the exponential link function used afterwards.
#' # Each batch has an nObservations x nParameters model matrix and
#' #  an nSamples x nParameters posterior matrix
#' effects_list <- list(
#'     main_fx  = list(X = model.matrix(~X1*X2, data = test_data),
#'                     B = posterior_samples$Beta),
#'     other_fx = list(X = model.matrix(~0 + grp, data = test_data),
#'                     B = posterior_samples$Gamma)
#'
#' predicted_values <- stan_yhat(effects_list, TRUE, exp)
#'
#' # Numerical example:
#' effects_list <- list(
#'     fx = list(
#'         X = array(runif(25), c(5, 5)),
#'         B = array(runif(1000), c(200, 5))
#'     )
#' )
#' Y <- stan_yhat(effects_list, TRUE, sigmoid)
stan_yhat <- function(effects_list, link = FALSE, link_fun = identity) {

  # checked named list
  fx_names <- names(effects_list)
  if (length(unique(fx_names)) != length(fx_names) || is.null(fx_names)) {
    fx_names <- paste0(fx_names, 1:length(fx_names))
  }

  null_names <- fx_names == ""
  if (any(null_names)) {
    fx_names[null_names] <- paste0('myfx_mejr', 1:sum(null_names))
  }

  for (i in seq_len(length(effects_list))) {
    effects_list[[i]]$valname <- fx_names[i]
  }

  # process each batch
  yhat <- lapply(effects_list, stan_yhat_i)

  # override names
  src_names <- c('__.samp key', '__.I key')
  new_names <- c("sample_i", "data_i")

  if (link) {
    # merge batches
    n_fx <- c(lapply(yhat, nrow), recursive = TRUE)
    if (!all(n_fx[1] == n_fx)) stop('Cannot link if nrows are not the same.')
    yhat <- multi_merge(yhat, by = src_names)
    valname <- names(yhat)[!names(yhat) %in% src_names]
    setnames(yhat, src_names, new_names)
    yhat[, yhat := link_fun(rowSums(.SD[, valname, with = FALSE]))]
    setcolorder(yhat, c(new_names, 'yhat', valname))
    setkey(yhat, NULL)
  } else {
    # keep as list
    lapply(yhat, function(i) {
      valname <- names(i)[!names(i) %in% src_names]
      i[, yhat := link_fun(get(valname))]
      setnames(i, src_names, new_names)
      setcolorder(i, c(new_names, 'yhat', valname))
      setkey(i, NULL)
    })
  }

  return(yhat)
}

alpha_override = function() {
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1)),
                  fill = ggplot2::guide_legend(override.aes = list(alpha = 1)))
}


color_10 <- function(n = 2)
{
  set <- c(
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf"
  )

  rep_len(set, n)
}


color_override <- function(n = 4, values, fill = FALSE) {
  if (missing(values))
    values <- color_10(n)

  if (fill) {
    return(ggplot2::scale_fill_manual(values = values))
  } else {
    return(ggplot2::scale_color_manual(values = values))
  }
}


#' RStan plots
#'
#' Saves a series of plots from a fitted rstan object
#'
#' @param stan_obj The fitted object from \code{rstan::stan}
#' @param pars A character vector of parameter names to plot. Defaults to working directory.
#' @param out Path for where to save the plots
#' @param label a character string to attach to file names
#' @param inc_warmup Include warmup for all plots. Defaults to \code{FALSE}
#'
#' @return nothing. Prints plots to disk
#'
#' @examples
#' stan_fit <- stan_mejr()
#' stan_plots_mejr(stan_fit$stan_mcmc)
stan_plots_mejr <- function(
  stan_obj,
  pars,
  out,
  label,
  inc_warmup = FALSE,
  max_pars = 9) {

  if (missing(pars)) {
    pars <- stan_obj@sim$pars_oi
    pars <- pars[!pars %in% c("log_lik")]
  }

  if (missing(label)) {
    label <- "mejr_model"
  }

  use_colors <- color_10(stan_obj@sim$chains)

  if (length(use_colors) > 1) {
    alpha_lvl <- 0.33
  } else {
    alpha_lvl <- 1
  }

  old_warn <- getOption("warn")
  options(list(warn = -1))

  save_plots <- function(p) {
    if (any(grepl("\\[[0-9]+", p))) {
      n_plots <- length(p)
      n_pages <- ceiling(n_plots / max_pars)
    } else {
      full_p <- stan_obj@sim$fnames_oi
      p_idx <- p == sub("\\[[0-9,]+\\]", "", full_p)
      if (any(p_idx)) {
        p <- full_p[p_idx]
      } else {
        return(invisible(NULL))
      }
      n_plots <- length(p)
      n_pages <- ceiling(n_plots / max_pars)
    }

    ncol <- ceiling(sqrt(max_pars))
    plot_list <- list()

    for (p0 in seq(1, max_pars*n_pages, max_pars)) {
      p1 <- min(c(p0 + (max_pars - 1), n_plots))
      par_set <- p[p0:p1]

      pts <- suppressMessages(rstan::stan_plot(
        stan_obj,
        pars = par_set,
        inc_warmup = inc_warmup,
        fill_color = use_colors[1],
        ci_level = 2/3, outer_level = 0.95))

      trc <- suppressMessages(rstan::stan_trace(
        stan_obj,
        pars = par_set,
        alpha = alpha_lvl,
        inc_warmup = inc_warmup,
        ncol = ncol
      )+ alpha_override()+ color_override(values = use_colors, fill = FALSE))

      den <- suppressMessages(rstan::stan_dens(
        stan_obj,
        pars = par_set,
        alpha = alpha_lvl,
        separate_chains = TRUE
      )+ alpha_override()+ color_override(values = use_colors, fill = TRUE))

      acr <- rstan::stan_ac(
        stan_obj,
        pars = par_set,
        lags = 6,
        ncol = ncol,
        partial = FALSE,
        color = "white")

      plot_list <- c(plot_list, list(pts, trc, den, acr))
    }

    return(plot_list)
  }

  plot_list <- list()

  if (is.list(pars)) {
    for (i in 1:length(pars))
      plot_list <- c(plot_list, save_plots(pars[[i]]))
  } else {
    for (p in pars)
      plot_list <- c(plot_list, save_plots(p))
  }

  options(list(warn = old_warn))

  if (length(plot_list) == 0)
    return(NULL)

  message("busy arranging grobs...")
  plts <- gridExtra::marrangeGrob(plot_list, ncol = 2, nrow = 2)

  if (missing(out)) {
    return(plts)
  } else {
    message("writing pdf...")
    graphics.off()
    pdf(
      file = file.path(
        out,
        paste0(label, "-stan_plots.pdf")
      ),
      width = 16.5354,
      height = 11.6929,
      onefile = TRUE)
    grid::grid.draw(plts)
    dev.off()
  }


  return(invisible(NULL))
}
