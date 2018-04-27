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
#' 1]}.
#' @return Numeric range
#' @export
#' @examples
#' x <- qnorm(seq(1e-04, .9999, length.out=1001))
#' # x <- c(seq(0,.5,length.out=250), rep(.5, 480), seq(.5,1, length.out=250))
#' hdi_95 <- hdi(x, .95)
#' hdi_50 <- hdi(x, .50)
#'
#' hist(x, br=50)
#' abline(v=hdi_95, col="red")
#' abline(v=hdi_50, col="green")
hdi <- function(x, prob = 0.95, warn = TRUE) {
  len <- length(x)

  if (len < 3) {
    warning("length of `x` < 3.", " Returning NAs")
    return(c(NA_integer_, NA_integer_))
  }

  x_sort <- sort(x)
  window_size <- as.integer(floor(prob * length(x_sort)))

  if (window_size < 3) {
    warning("window_size < 3.", " `prob` is too small.", " Returning NAs")
    return(c(NA_integer_, NA_integer_))
  }

  lower <- seq_len(len - window_size)
  upper <- window_size + lower

  # vectorized difference between edges of cumulative distribution based on
  # scan_length. Values are arranged from left to right scanning.
  window_width_diff <- x_sort[upper] - x_sort[lower]

  if (warn) {
    # find minimum of width differences, check for multiple minima
    min_i <- which(window_width_diff == min(window_width_diff))
    n_candies <- length(min_i)

    if (n_candies > 1) {
      warning(
        "Multiple candidate thresholds found for HDI.",
        " Choosing the middle of possible limits.")
      # if more than one minimum, get average index
      get_diff <- c(1, min_i[2:n_candies] - min_i[1:(n_candies - 1)])
      if (any(get_diff != 1)) {
        stop_i <- which(get_diff != 1) - 1
        min_i <- min_i[1:stop_i[1]]
      }
      min_i <- floor(mean(min_i))
    }
  } else {
    min_i <- which.min(window_width_diff)
  }

  # get values based on minimum
  c(x_sort[min_i], x_sort[upper[min_i]])
}

prop_in_range <- function(x, lower, upper) {
  sum(x > lower & x < upper) / length(x)
}

#' Posterior intervals
#'
#' Returns cutoff points from a posterior distribution
#' @return data.table
#'
#' @param x Vector of numeric values. Typically a posterior sample.
#' @param mid Central tendency estimator. Defaults to \code{"median"}. Other options include
#' \code{c("mean", "mode")}.
#' @param widths interval widths
#' @param adj Bandwidth adjustment used only with the \code{"mode"} estimator. See \link{dmode}.
#' @param rope Region of practical equivalence. Check how much of the distribution is within rope
#' value.
#' @param warn Turn off warning for flat intervals found (multiple possible values)
#' @param int interval type, either "hdi" or "ci"
#'
#' @examples
#' x <- rpois(5000, 15)
#' ints <- post_int(x)
#' hist(x, br=50)
#' abline(v=ints$c, col="cyan")
#' abline(v=ints[, c("l.wide", "r.wide")], col="magenta")
#'
#' post_int(x, "median")
#' post_int(x, "mean")
#' post_int(x, "mode", adj=2, rope = c(14, 16))
#' @export
post_int <- function(x, mid = c("median", "mean", "mode"),
                     int = c("hdi", "ci"), widths = c(.50, .95),
                     adj = 1.5, rope = NULL, warn = FALSE) {
  if (!is.vector(x)) {
    stop("`x` must be a vector.")
  }
  mid <- match.arg(mid)
  int <- match.arg(int)

  center <- NA_real_
  scale <- NA_real_

  # measure of central tendency
  switch(mid,
    "mean" = {
      center <- mean(x)
      scale <- sd(x)
    },
    "median" = {
      center <- median(x)
      scale <- mad(x)
    },
    "mode" = {
      center <- dmode(x, adjust = adj)
      scale <- mad(x)
    },
    NULL)

  # interval function as an unevaluated function call
  int_fun <- switch(
    int,
    "ci" = {
      match.call(
        posterior_interval,
        call("posterior_interval",
          object = quote(as.matrix(x)), prob = quote(w)))
    },
    "hdi" = {
      match.call(hdi, call("hdi",
        x = quote(x),
        prob = quote(w), warn = warn))
    },
    NULL
  )

  # mass within 1 sd range
  area_sd <- prop_in_range(x, center - scale, center + scale)

  # width identifiers
  widths <- sort(widths)
  width_ids <- sprintf("%.0f", widths * 100)
  width_ids[which.max(widths)] <- "wide"

  # left/right values for each width
  intervals <- data.table::rbindlist(
    Map(
      function(w, i) {
        ints <- structure(Reduce(cbind, eval(int_fun)),
          dimnames = list(NULL, c("l", "r")))
        data.table::data.table(interval = i, ints)
      },
      c(widths, area_sd),
      c(width_ids, "sd"))
  )

  # reshape to multivariate format
  intervals <- intervals %>%
    .[order(-interval), central := mid] %>%
    data.table::dcast(central ~ interval,
      value.var = c("l", "r"),
      sep = ".") %>%
    .[, c := center] %>%
    set_col_order(c("central", "c", "l.sd", "l.wide", "r.sd", "r.wide"))

  # get rope
  if (!is.null(rope)) {
    if (is.character(rope) && rope == "max") {
        rope <- c(intervals$l.wide, intervals$r.wide)
    }
    if (length(rope) != 2) {
      stop("ROPE must be a lower and upper value")
    }
    rope <- sort(rope)
    rope_mass <- prop_in_range(x, rope[1], rope[2])
    intervals[, `:=`(l.rope = rope[1], r.rope = rope[2], rope = rope_mass)]
  }

  as.data.frame(intervals)
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
  x <- trim_ends(x, 0.05)
  d <- density(x, n = 1000, bw = "SJ", adjust = adjust)
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
    n_vec <- matrix(rep(n_vec, each = nrow(sd_vec)), nrow = nrow(sd_vec))
  }
  sqrt(rowSums((n_vec - 1) * sd_vec^2) / (rowSums(n_vec) - ncol(n_vec)))
}

#' Trim extreme values at each end of a vector.
#'
#' @param x A [numeric] vector
#' @param trim Proportion of vector length to trim. Must be between 0 and 1.
#' E.g., a value 0.05 (default) trims 2.5\% off each end of a sorted vector.
#' @param na.rm omit `NA` values. May result in different size vector.
#'
#' @return A [numeric] vector in the original order of `x`, but with trimmed
#' values as `NA` if `na.rm=TRUE` or with these values removed if `FALSE`
#' (which will result in a different sized vector from the input).
#' @export
#'
#' @examples
#' x <- rgamma(10000, 1, 1)
#' range(x)
#' length(x)     # <- 10000
#' sum(is.na(x)) # <- 0
#'
#' t <- trim_ends(x, trim = 0.1)
#' range(t)
#' length(t)     # <- 9000
#' sum(is.na(t)) # <- 0
#'
#' t <- trim_ends(x, 0.1, na.rm = FALSE)
#' range(t, na.rm = TRUE)
#' length(t)     # <- 10000
#' sum(is.na(t)) # <- 1000
trim_ends <- function(x, trim = 0.05, na.rm = TRUE) {
  if (trim < 0 | trim > 1) {
    stop("trim amount must be between 0 and 1")
  }

  which_na <- is.na(x)

  len <- sum(!which_na)

  trim_size <- as.integer(floor((len * trim) / 2))

  if (trim_size > 0L) {
    real_vals <- which(!which_na)
    sort_order_real <- order(x[real_vals], na.last = TRUE)
    real_vals <- real_vals[sort_order_real]
    trim_index <- seq_len(trim_size)
    which_na[real_vals[trim_index]] <- TRUE
    which_na[real_vals[len - trim_index + 1]] <- TRUE
  }

  if (!na.rm) {
    x[which_na] <- NA
    return(x)
  }

  x[!which_na]
}

#' @export
apply_contrast <- function(contrast_coef, ...) {
  if (dot_dot_len(...) != length(contrast_coef)) {
    stop("Num. coefficients must equal num vectors.")
  }
  if (sum(contrast_coef) != 0) {
    warning("Sum of coefficients not equal to zero.")
  }
  as.vector(cbind(...) %*% contrast_coef)
}

