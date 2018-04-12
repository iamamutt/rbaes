
get_unit_attr <- function(x) {
    adj <- attr(x, 'adj')
    m <- attr(x, 'center')
    s <- attr(x, 'scale')

  if (any(is.null(c(adj, m, s)))) {
    stop('could not find attributes from a unit_scale object')
  }

  nlist(adj, m, s)
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
unit_scale <- function(x, y=NULL, na.rm=TRUE) {
  if (is.null(y)) {
    x_mean <- mean(x, na.rm=na.rm)
    s_std <- 2 * sd(x, na.rm=na.rm)
    x_adj <- 0.5
  } else {
    y_attr <- get_unit_attr(y)
    x_mean <- y_attr$m
    s_std <- y_attr$s
    x_adj <- y_attr$adj
  }

  z <- x_adj + ((x - x_mean) / s_std)

  attr(z, 'center') <- x_mean
  attr(z, 'scale') <- s_std
  attr(z, 'adj') <- x_adj
  return(z)
}

#' Reverse unit scaling
#'
#' @param x object to be scaled back to original scale
#' @param y optional object that has center, scale, and adj attributes to use on x. If y is not
#'   used, these attributes will be taken from x.
#'
#' @return numeric vector in original scale
#' @export
#'
#' @examples
#' z <- rnorm(10000, 100, 15)
#' y <- unit_scale(z)
#' rev_unit_scale(y)
rev_unit_scale <- function(x, y=NULL, na.rm=FALSE) {
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
