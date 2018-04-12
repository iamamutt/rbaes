
#' Highest density interval
#'
#' This is a function that will calculate the highest density interval from a posterior sample.
#'
#' The default is to calcualte the highest 95 percent interval. It can be used with any numeric
#' vector instead of having to use one of the specific MCMC classes. This function has been adapted
#' from John K. Kruschke (2011). Doing Bayesian Data Analaysis: A Tutorial with R and BUGS.
#'
#' @param x Numeric vector of a distribution of data, typically a posterior sample
#' @param prob Width of the interval from some distribution. Defaults to \code{0.95}.
#' @param warn Option to turn off multiple sample warning message Must be in the range of \code{[0,
#'   1]}.
#' @return Numeric range
#' @export
#' @examples
#' x <- qnorm(seq(1e-04, .9999, length.out=1001))
#' # x <- c(seq(0,.5,length.out=250), rep(.5, 480), seq(.5,1, length.out=250))
#' hdi_95 <- hdi(x)
#' hdi_50 <- hdi(x, .5)
#'
#' hist(x, br=50)
#' abline(v=hdi_95, col="red")
#' abline(v=hdi_50, col="green")
hdi <- function(x, prob = 0.95, warn = TRUE) {
  sorted_x <- sort(x)
  x_size <- length(sorted_x)
  window_size <- floor(prob * length(sorted_x))
  scan_index <- 1:(x_size - window_size)

  # vectorized difference between edges of cumulative distribution based on scan_length
  window_width_diff <- sorted_x[scan_index + window_size] - sorted_x[scan_index]

  # find minimum of width differences, check for multiple minima
  candidates <- which(window_width_diff == min(window_width_diff))
  n_c <- length(candidates)
  if (warn && n_c > 1) {
    warning(simpleWarning(paste0("Multiple candidate thresholds found for HDI,",
                                 "choosing the middle of possible limits.")))
  }

  # if more than one minimum, get average index
  if (length(candidates) > 1) {
    get_diff <- c(1, candidates[2:n_c] - candidates[1:(n_c - 1)])
    if (any(get_diff != 1)) {
      stop_i <- which(get_diff != 1) - 1
      candidates <- candidates[1:stop_i[1]]
    }
    min_i <- floor(mean(candidates))
  } else min_i <- candidates

  # get values based on minimum
  hdi_min <- sorted_x[min_i]
  hdi_max <- sorted_x[min_i + window_size]
  hdi_vals <- c(hdi_min, hdi_max)
  return(hdi_vals)
}

#' Posterior intervals
#'
#' Returns cutoff points from a posterior distribution
#' @return data.table
#'
#' @param x Vector of numeric values. Typically a posterior sample.
#' @param mid Central tendency estimator. Defaults to \code{"median"}. Other options include
#'   \code{c("mean", "mode")}.
#' @param widths interval widths
#' @param adj Bandwidth adjustment used only with the \code{"mode"} estimator. See \link{dmode}.
#' @param rope Region of practical equivalence. Check how much of the distribution is within rope
#'   value.
#' @param warn Turn off warning for flat intervals found (multiple possible values)
#' @param int interval type, either "hdi" or "ci"
#'
#' @examples
#' x <- rpois(1000, 15)
#' hist(x, br=50)
#' abline(v=post_int(x), col="cyan")
#'
#' post_int(x, "median")
#' post_int(x, "mean")
#' post_int(x, "mode", adj=2)
#' @export
post_int <- function(
  x,
  mid = "median",
  int = "ci",
  widths = c(.5, .90),
  adj = 1.5,
  rope = NULL,
  warn = FALSE)
{
  if (mid == 'mean') {
    m <- mean(x)
    s <- sd(x)
  } else if (mid == 'median') {
    m <- median(x)
    s <- mad(x)
  } else if (mid == 'mode') {
    m <- dmode(x, adjust = adj)
    s <- mad(x)
  } else {
    stop(sprintf('"%s" is not a measure of central tendency', mid))
  }

  if (int == 'ci') {
    fx <- function(x, p, w = NULL) {
      rstantools::posterior_interval(
        object = as.matrix(x), prob = p)
    }
  } else if (int == 'hdi') {
    fx <- function(x, p, w) {
      hdi(x = x, prob = p, warn = w)
    }
  } else {
    stop('unknown interval method')
  }

  s_int <- sum(x > m - s & x < m + s) / length(x)
  central <- data.table::data.table(
    interval = c(0, s_int, s_int),
    side = c('c', 'l', 'r'),
    type = c('mid', 'std', 'std'),
    y = c(m, m - s, m + s)
  )

  intervals <- data.table::data.table(
    widths,
    do.call(
    rbind, lapply(seq_along(widths), function(w) {
    fx(x, widths[w], warn)
  })))

  names(intervals) <- c('interval', 'l', 'r')

  intervals <- data.table::melt(
    intervals, id.vars = 'interval',
    variable.name = 'side',
    value.name = 'y')

  intervals$type <- int

  data.table::setcolorder(
    intervals,
    c('interval', 'side', 'type', 'y'))

  if (!is.null(rope)) {
    if (length(rope) != 2)
      stop("ROPE must be a lower and upper value")
    rope <- sort(rope)
    rope_int <- sum(x > rope[1] & x < rope[2]) / length(x)
    rope_p <- data.table::data.table(
      interval = c(rope_int, rope_int),
      side = c('l', 'r'),
      type = 'rope',
      y = rope)

  } else {
    rope_p <- NULL
  }

  out <- rbind(central, intervals, rope_p)
  out <- out[order(interval, type, side, y)]

  return(as.data.frame(out))
}


#' Mode from counting frequency
#'
#' Finds the most frequent value from a vector of integers
#'
#' @param x an integer vector
#'
#' @return scalar integer value
#' @export
#'
#' @examples
#' cmode(rpois(1000, 20))
cmode <- function(x) {
  x <- sort(x)
  y <- rle(x)
  y$values[which.max(y$lengths)]
}

#' Mode from density estimation
#'
#' Finds the mode using the \link{density} function and then obtains the maximum value.
#'
#' @param x Value vector. Numeric or integers.
#' @param adjust Bandwidth adjustment. See \link{density}.
#'
#' @examples
#' x <- rchisq(1000, 3)
#' hist(x, br=50)
#' abline(v = dmode(x), col = "red")
#' abline(v = median(x), col = "green")
#' abline(v = mean(x), col = "blue")
#' @export
dmode <- function(x, adjust = 1.5) {
  cut <- hdi(x, 0.99)
  d <- density(x, n = 1000, bw = "SJ", from = cut[1], to = cut[2], adjust = adjust)
  d$x[which.max(d$y)]
}


#' Pooled standard deviation
#'
#' Get pooled standard deviation from a vector of SDs and sample sizes
#'
#' @param sd_vec numeric vector of standard deviations
#' @param n_vec sample sizes corresponding to the standard deviations
#'
#' @return numeric value
#' @export
#'
#' @examples
#' pooled_sd(c(1, 5), c(100, 8))
pooled_sd <- function(sd_vec, n_vec) {
  if (is.vector(sd_vec)) {
    sd_vec <- matrix(sd_vec, ncol = length(sd_vec))
  }
  if (is.vector(n_vec)) {
    n_vec <- matrix(n_vec, ncol = length(sd_vec))
  }
  if (ncol(sd_vec) != ncol(n_vec)) {
    stop("SD vec must equal length of N vec")
  }
  if (nrow(n_vec) == 1 & nrow(sd_vec) > 1) {
    n_vec <- matrix(rep(n_vec, each = nrow(sd_vec)),
                    nrow = nrow(sd_vec))
  }
  sqrt(rowSums((n_vec - 1) * sd_vec ^ 2) /
         (rowSums(n_vec) - ncol(n_vec)))
}
