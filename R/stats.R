
is.atomic_vec <- function(x) {
  is.vector(x) && is.atomic(x)
}

vec_to_mat <- function(x, row_vec = FALSE) {
 if (is.matrix(x)) {
   return(x)
 }

 if (!is.atomic_vec(x)) {
   stop("x is not a vector.")
 }

 n <- length(x)

 if (n < 1) {
   return(as.matrix(get(typeof(x))()))
 }

 if (row_vec) {
   dim(x) <- c(1L, n)
 } else {
   dim(x) <- c(n, 1L)
 }
 x
}

rep_mat <- function(x, n, dim) {
  if (dim == 1) {
    matrix(rep(x, each = n), nrow = n)
  } else {
    matrix(rep(x, each = n), ncol = n)
  }
}


#' Hedges G
#' @param n_vec sample sizes corresponding to the standard deviations
#'
#' @return numeric value
#' @export
#'
#' @examples
#' pooled_sd(c(1, 5), c(100, 8))
pooled_sd <- function(sd, n) {
  sd_mat <- vec_to_mat(sd, row_vec = TRUE)
  n_mat <- vec_to_mat(n, row_vec = TRUE)

  cols_sd <- ncol(sd_mat)
  cols_n <- ncol(n_mat)

  if (cols_sd != cols_n) {
    if (cols_n == 1) {
      n_mat <- rep_mat(n_mat, cols_sd, 2)
      cols_n <- ncol(n_mat)
    } else {
      stop("SD vec must equal length of N vec")
    }
  }

  if (nrow(n_mat) == 1 & nrow(sd_mat) > 1) {
    n_mat <- rep_mat(n_mat, nrow(sd_mat), 1)
  }

  numer <- rowSums((n_mat - 1) * sd_mat^2)
  denom <- rowSums(n_mat) - cols_n
  sqrt(numer / denom)
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

