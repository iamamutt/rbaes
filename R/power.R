
ss_one_sample <- function(n = NULL, min_n = 30, max_n = 100, increments = 10)
{
  if (is.null(n)) {
    n <- seq(min_n, max_n, increments)
  }
  return(sapply(n, function(x) list(list(n = x))))
}

ss_two_sample <- function(n1 = NULL, n2 = NULL, min_n = 30, max_n = 100, increments = 10)
{
  # number of attempts at reaching power given sample size at certain increments
  if (is.null(n1))
    n1 <- seq(min_n, max_n, increments)
  if (is.null(n2))
    n2 <- seq(min_n, max_n, increments)
  n1_l <- length(n1)
  n2_l <- length(n2)

  # if number of sample size attempts is not equal for both groups then repeat for the one with less tries
  if (n1_l != n2_l) {
    if (n1_l < n2_l) {
      n1 <- c(n1, rep(rev(n1), length.out = n2_l)[(n1_l + 1):n2_l])
    } else {
      n2 <- c(n2, rep(rev(n2), length.out = n1_l)[(n2_l + 1):n1_l])
    }
  }
  ss_attempts <- sapply(1:max(c(n1_l, n2_l)), function(x) list(list(n1 = n1[x], n2 = n2[x])))
  return(ss_attempts)
}

refit_stan_template <- function()
{
  list(data = NULL, update_stan_args = NULL)
}

add_test_item <- function(test_name, passed, check=FALSE, test_set='pwr_test', store_value=NULL)
{
  if (!is.null(store_value) & is.recursive(store_value)) {
    stop('store_value must be an atomic type as seen in is.atomic')
  }

  if (!is.logical(passed)) {
    stop('passed must be a logical value indicating if the test passed or not')
  }

  if (!is.logical(check)) {
    stop('check must be a logical value indicating if the test should be evaluated for stopping')
  }

  slot <- data.table(set=test_set, test=test_name, passed=passed, check=check, value=store_value)

  return(slot)
}

bayesian_power <- function(stanfit, draw_fun, make_fun, test_fun,
                           sample_sizes = NULL, stan_args = NULL,
                           power_goal = 0.8, max_draws = 250, min_draws = 25,
                           seed = 19851905, print_progress = 25, fileprint=NULL)
{
  set.seed(seed)
  mcmc_samples <- sum(unlist(lapply(stanfit@sim$permutation, length)))
  n_passes <- length(sample_sizes)

  if (is.null(min_draws) | is.na(min_draws)) min_draws <- max_draws
  if (min_draws > max_draws) stop('min_draws must be less than max_draws')

  draw_order <- sample(mcmc_samples, max_draws, ifelse(max_draws > mcmc_samples, TRUE, FALSE))
  pwr_record <- power_record()

  for (pass in 1:n_passes) {
    # seed is reset to append simulation values from last pass
    set.seed(seed)
    sample_size <- sample_sizes[[pass]]

    for (draw in 1:max_draws) {
      pars <- draw_fun(stanfit, draw_order[draw])
      stan_data <- make_fun(pars, sample_size)
      refitted <- power_stan(stan_data, stan_args)
      test_results <- test_fun(refitted)

      # save record of results from current draw
      pwr_record <- power_record(
        record=pwr_record, draw=draw, pass=pass, power_goal=power_goal,
        sample_size=sample_size, test_results=test_results, print_progress=print_progress)

      # print updated status given print_progress
      if (!is.null(print_progress)) {
        if (pwr_record[[pass]]$counter %% print_progress == 0)
          print_power_status(record, fileprint)
      }

      # power reached early, return from function
      if (pwr_record[[pass]]$lower_reached) {
        print_power_status(pwr_record, fileprint)
        message(sprintf('Power stoping point reached early after %d draws.'), draw)
        return(pwr_record)
      }

      # start checking if need to quit and move on to next sample size
      if (draw > min_draws) {
        upper_less_than_goal <- !pwr_record[[pass]]$upper_reached
        will_fail <- pwr_record[[pass]]$p_goal_met < .05
        if (upper_less_than_goal | will_fail) {
          print_power_status(pwr_record, fileprint)
          message(sprintf('Likelihood of reaching goal is <5%% after %d draws. Trying next.', draw))
          break
        }
      }
    }

    # lower not reached but center power value reached
    if (pwr_record[[pass]]$center_reached) {
      print_power_status(pwr_record, fileprint)
      message('Power stoping point reached.')
      return(pwr_record)
    }
  }

  # all attempts failed
  print_power_status(pwr_record, fileprint)
  warning('Power was never reached!')
  return(pwr_record)
}

power_stan <- function(data, stan_args = NULL)
{
  arg_updates <- data$update_stan_args

  if (is.null(stan_args)) {
    stan_args <- arg_updates
  } else {
    if (!is.null(arg_updates)) {
      which_updates <- names(arg_updates) %in% names(stan_args)
      # first overwrite stan_args from args_updates
      if (any(which_updates)) {
        to_update <- names(arg_updates)[which_updates]
        stan_args[names(stan_args) %in% to_update] <- arg_updates[which_updates]
      }
      # then append any additional args from args_updates
      if (any(!which_updates)) {
        stan_args <- c(stan_args, arg_updates[!which_updates])
      }
    }
  }

  if (is.null(stan_args))
    stop("No stan arguments found in stan_args or data$update_stan_args. Need at least `file`` or `model_code`")

  stan_args$data <- data$data
  old_opt <- rstan::rstan_options(auto_write = TRUE)
  stanfit <- do.call(rstan::stan, stan_args)
  rstan::rstan_options(auto_write = old_opt)

  return(stanfit)
}

power_record_blank <- function()
{
  list(list(pass = 0, draw = 0, counter = 0, total = NA, sample_size = NA, power_goal = NA, data = list()))
}

power_record <- function(record = NULL, draw, pass, power_goal, sample_size, test_results, print_progress = NULL)
{
  if (is.null(record)) {
    record <- list()
    return(record)
  }

  r <- length(record)

  if (r > 0){
    counter <- record[[r]]$counter
  } else {
    counter <- 0
  }

  if (draw == 1) {
    record <- c(record, power_record_blank())
    r <- length(record)
    record[[r]]$counter <- counter
    record[[r]]$pass <- pass
    record[[r]]$sample_size <- sample_size
    record[[r]]$power_goal <- power_goal
  }

  if (missing(test_results)) {
    stop('results argument empty')
  }

  if (is.data.table(test_results)) {
    # test function returned a single test
    test_results[, test_id := 1]
  } else {
    # test function returned a list of tests
    test_results <- data.table::rbindlist(test_results, idcol = 'test_id')
  }

  test_results[, value_id := 1:.N, .(test_id, set)]

  if (draw==1) {
    total <- test_results[, numeric(length(unique(test_id)))]
  } else {
    total <- record[[r]]$total
  }

  total <- total + test_results[, as.integer(passed[1]), test_id][[2]]
  successes <- total + 1
  failures <- 1 + draw - total

  test_results[, c('power', 'pow_lower', 'pow_upper') := compute_power_interval(
    successes[test_id], failures[test_id], mass = 0.95), test_id]

  power_goals_met <- power_results_tests(test_results, power_goal, draw)

  # update records
  record[[r]]$draw <- draw
  record[[r]]$counter <- record[[r]]$counter + 1
  record[[r]]$total <- total
  record[[r]]$data[[draw]] <- as.data.frame(test_results)
  record[[r]][names(power_goals_met)] <- power_goals_met

  return(record)
}


power_results_tests <- function(results, power_goal, draw)
{
  max_set_pwr <- results[check==TRUE, .(row=.I[which.max(pow_lower)]), .(set)]
  last_best <- results[max_set_pwr$row, .(test, power, pow_lower, pow_upper)]
  lower_reached <- last_best[, all(pow_lower >= power_goal)]
  center_reached <- last_best[, all(power >= power_goal)]
  upper_reached <- last_best[, min(pow_upper) > power_goal]

  p_goal_met <- last_best[, pbeta(power_goal,
                                  pow_upper * draw + 1,
                                  (1 - pow_upper) * draw + 1)]

  p_goal_met <- max(1 - p_goal_met)

  nlist(last_best, lower_reached, center_reached, upper_reached, p_goal_met)
}

compute_power_interval <- function(successes, failures, mass=0.95)
{
  if (mass <= 0 | mass >= 1) stop("HDI interval must have a mass between 0 and 1")

  opt_result <- optimize(beta_optimization_fn, c(0, 1 - mass),
                         a = successes, b = failures, mass = mass, tol = 1e-09)

  left_tail_p <- opt_result$minimum

  data.table(
    power=successes / (successes + failures),
    pow_lb=qbeta(left_tail_p, successes, failures),
    pow_ub=qbeta(mass + left_tail_p, successes, failures))
}

beta_optimization_fn <- function(x, a, b, mass=0.95)
{
  qbeta(x + mass, a, b) - qbeta(x, a, b)
}

print_power_status <- function(record, filename=NULL)
{
  if (!is.null(filename)) {
    plog_file <- file(filename, open='wt')
    sink(plog_file)
    sink(plog_file, type='message')
  }

  n_recs <- length(record)
  sample_size <- paste(unlist(record[[n_recs]]$sample_size), collapse=',')
  pass <- record[[n_recs]]$pass
  draw <- record[[n_recs]]$draw
  counter <- record[[n_recs]]$counter
  last_best <- record[[n_recs]]$last_best
  power_goal <- record[[n_recs]]$power_goal
  p_goal_met <- record[[n_recs]]$p_goal_met

  message(sprintf(
    paste0("Status:\n  Attempt %d, Draw %d, Iter %d, Sample Size %s\n",
           "  prob. of achieving power of %.2f is currently %.2f"),
    pass, draw, counter, sample_size, power_goal, p_goal_met))

  last_best[, message(sprintf(
    "  Power: %.2f [%.2f, %.2f]\t%s",
    power, pow_lower, pow_upper, test)), test]

  if (!is.null(filename)) {
    sink(type='message')
    sink()
  }

  return(invisible())
}

# power tests ---------------------------------------------------------------------------------


test_interval <- function(x, null_value = 0, type = c("!", "l>", "r<", "l>=", "r<=", '=='))
{
  if (is.vector(x))
    x <- matrix(x, ncol = 2)

  # make sure lower and upper are sorted
  if (!all(apply(x, 1, function(i) i[2] >= i[1])))
    stop('x must be an interval with the first value <= the second value')

  lhs <- x[, 1]
  rhs <- x[, 2]
  ttype <- type[1]

  if (ttype %in% c("!", "!=", "neq")) {
    test <- rhs < null_value | lhs > null_value
  } else if (ttype %in% c("l>", "lgt")) {
    test <- lhs > null_value
  } else if (ttype %in% c("l>=", "lgeq")) {
    test <- lhs >= null_value
  } else if (ttype %in% c("l<", "llt")) {
    test <- lhs < null_value
  } else if (ttype %in% c("l<=", "lleq")) {
    test <- lhs <= null_value
  } else if (ttype %in% c("r>", "rgt")) {
    test <- rhs > null_value
  } else if (ttype %in% c("r>=", "rgeq")) {
    test <- rhs >= null_value
  } else if (ttype %in% c("r<", "rlt")) {
    test <- rhs < null_value
  } else if (ttype %in% c("r<=", "rleq")) {
    test <- rhs <= null_value
  } else if (ttype %in% c("==", "eq")) {
    # ROPE
    if (length(null_value)!=2)
      stop('null_value must be a lower and upper region for this option')
    test <- lhs >= null_value[1] & rhs <= null_value[2]
  } else {
    stop("unknown interval test type: ", ttype)
  }
  return(test)
}


# stats ---------------------------------------------------------------------------------------

bayes_p <- function(x, null = 0)
{
  n <- length(x)
  p <- sum(x > null)/n
  return(c(leq_null = 1 - p, gt_null = p))
}

cohens_d <- function(x, sdx, y = NULL, sdy = NULL, nx = NULL, ny = NULL, one_sample_test_pt = 0)
{
  null_nx <- is.null(nx)
  null_ny <- is.null(ny)
  if (is.null(y)) {
    (x - one_sample_test_pt)/sdx
  } else {
    if (is.null(sdy))
      stop("y arg specified without specifying sdy")
    if (null_nx & null_ny) {
      # equal sample sizes
      nx <- 2
      ny <- 2
    } else if (xor(null_nx, null_ny)) {
      stop("Only one sample size specified. If sample sizes are equal, leave nx, ny = NULL")
    }
    sd_vec <- cbind(sdx, sdy)
    n_vec <- cbind(nx, ny)
    sd_p <- pooled_sd(sd_vec, n_vec)
    fx_size <- (x - y)/sd_p
  }
}

pooled_sd <- function(sd_vec, n_vec)
{
  if (is.vector(sd_vec)) {
    sd_vec <- matrix(sd_vec, ncol = length(sd_vec))
  }
  if (is.vector(n_vec)) {
    n_vec <- matrix(n_vec, ncol = length(sd_vec))
  }
  if (ncol(sd_vec) != ncol(n_vec))
    stop("SD vec must equal length of N vec")
  if (nrow(n_vec) == 1 & nrow(sd_vec) > 1)
    n_vec <- matrix(rep(n_vec, each = nrow(sd_vec)), nrow = nrow(sd_vec))
  sqrt(rowSums((n_vec - 1) * sd_vec^2)/(rowSums(n_vec) - ncol(n_vec)))
}


# example functions ---------------------------------------------------------------------------



one_sample_mvt_draw <- function(stanfit, i)
{
  post <- rstan::extract(stanfit, pars = c("mu", "Sigma", "nu"))
  pars <- list(mu = post$mu[i, ], Sigma = post$Sigma[i, , ], nu = post$nu[i])
  return(pars)
}

#' One sample multivariate-t
#'
#' Simulation function for one sample multivariate t
#'
#'
#' @param pars list of parameters needed for simulation and to create objects to pass to stan
#'
#' @return sim_fit_template
one_sample_mvt_make <- function(pars, sample_size_list = NULL)
{
  if (is.null(sample_size_list)) {
    # use sample size from data if only determining current power
    N <- pars$N
  } else {
    # since this is a one sample, it'll have a single slot `n`
    N <- sample_size_list$n
  }
  # the rest of the code uses parameter values to simulate a new stan data list
  Sigma <- pars$Sigma
  K <- ncol(Sigma)
  nu <- pars$nu
  mu <- pars$mu
  # i.i.d student-t values
  mu_t_zero <- array(rt(N * K, df = nu), c(N, K))
  # reshape t-values according to Sigma and rescale by mu
  Y <- t(apply(mu_t_zero, 1, function(x) {
    mu + x %*% chol(Sigma)
  }))
  # use the template to check for some slots not used, such as initialization values
  stan_data <- refit_stan_template()
  stan_data$data <- nlist(Y, K, N)
  inits <- function(chain_id) {
    list(mu_adjustment = c(-0.25, 0.25, 0.25), nu_minus_one = 36, sigma = c(0.7, 0.88, 0.85), chol_corr = diag(3))
  }
  stan_data$update_stan_args <- list(init = inits)
  return(stan_data)
}


one_sample_mvt_tests <- function(refitted)
{
  # extract paramters
  post <- rstan::extract(refitted, pars = c("mu", "sigma"))
  mu <- post$mu
  sigma <- post$sigma

  # mean differences
  emb_fam <- mu[, 2] - mu[, 1]
  par_fam <- mu[, 3] - mu[, 1]
  par_emb <- mu[, 3] - mu[, 2]

  # effect sizes
  fx_size_EF <- cohens_d(mu[, 2], sigma[, 2], mu[, 1], sigma[, 1])
  fx_size_PF <- cohens_d(mu[, 3], sigma[, 3], mu[, 1], sigma[, 1])

  # test values
  hdi_EF <- rbaes::hdi(emb_fam, 0.9)
  p_EF <- bayes_p(emb_fam)[1]
  hdi_D_EF <- rbaes::hdi(fx_size_EF, 0.9)
  hdi_PF <- rbaes::hdi(par_fam, 0.9)
  p_PF <- bayes_p(par_fam)[1]
  hdi_D_PF <- rbaes::hdi(fx_size_PF, 0.9)
  hdi_PE <- rbaes::hdi(par_emb, 0.9)
  p_PE <- bayes_p(par_emb)
  p_PE_test <- !as.logical(p_PE[1] < 0.05 | p_PE[2] < 0.05)

  return(
    list(
      add_test_item(test_name = "HDI E-F > 0", passed = test_interval(hdi_EF), check = TRUE, test_set = "EvF", store_value = hdi_EF),
      add_test_item("p E-F < .05", p_EF < 0.05, check = TRUE, "EvF", p_EF),
      add_test_item("HDI D > 0.1 | E,F", test_interval(hdi_D_EF, null = 0.1, type = "l>"), FALSE, "EvF", hdi_D_EF),
      add_test_item("HDI P-F > 0", test_interval(hdi_PF), check = TRUE, "PvF", hdi_PF),
      add_test_item("p P-F < .05", p_PF < 0.05, check = TRUE, "PvF", p_PF),
      add_test_item("HDI D > 0.1 | P,F", test_interval(hdi_D_PF, null = 0.1, type = "l>"), FALSE, "PvF", hdi_D_PF),
      add_test_item("HDI P-E == 0", !test_interval(hdi_PE), FALSE, "PvE", hdi_PE),
      add_test_item("p P-E > .05", p_PE_test, FALSE, "PvE", p_PE)
    )
  )
}