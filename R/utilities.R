
example_stanreg <- function() {
  d <- rstanarm::radon
  rstanarm::stan_glmer(
    log_radon ~ floor * log_uranium + (1 | county),
    data = d,
    iter = 1000, chains = 2)
}

n_post_draws <- function(stanfit) {
  pss <- stanfit$stanfit@sim$n_save
  sum(pss - stanfit$stanfit@sim$warmup2)
}


ppi <- function(
  data,
  fit,
  fun = identity,
  form = NULL,
  prob = 0.9,
  n = 1)
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
print_stanarm <- function
(
  stanreg,
  header='Stan model summary'
)
{
  txt <- capture.output(print(stanreg, digits = 3))
  est <- which(txt == 'Estimates:'):length(txt)
  cat('\n')
  mejr::print_sec(header)
  cat(stringr::str_c(txt[est], collapse = '\n'))
  cat('\n\n')
  print(rstanarm::prior_summary(stanreg))
  cat('\n')
  return(invisible(NULL))
}
