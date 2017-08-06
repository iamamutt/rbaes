
`%?n%` <- function(obj, names) {
  obj_names <- names(obj)
  if (all(names %in% obj_names)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

dots2list <- function(...) {
  eval(substitute(alist(...)))
}

argval2char <- function(...) {
  as.character(match.call())[-1L]
}

symbol2char <- function(...) {
  lapply(dots2list(...), deparse)
}

formula2char <- function(x){
  Reduce(paste, deparse(x))
}

nlist <- function(...) {
  nms <- as.character(match.call())[-1L]
  out <- list(...)
  named <- names(out)
  if (is.null(named)) { # all unnamed
    names(out) <- nms
  } else {
    which_named <- nzchar(named)
    if (!all(which_named)) { # partial named
      names(out)[!which_named] <- nms[!which_named]
    }
  }
  return(out)
}

lextract <- function(x, ...) {
  entries <- symbol2char(...)
  get_from_list <- function(l, n) {
    if (!l %?n% n)
      return(NULL)
    l[[n]]
  }
  lapply(x, function(i) {
    Reduce(get_from_list, entries, init = i, accumulate = FALSE)
  })
}

dtbl2list <- function(data, ...) {

  if (!is.data.table(data)) {
    dt <- as.data.table(data)
    dtbl <- FALSE
  } else {
    dt <- copy(data)
    dtbl <- TRUE
  }

  by_cols <- unlist(symbol2char(...))

  if (!dt %?n% by_cols) {
    stop(sprintf(
      'check that columns exist:\n  %s',
      paste(by_cols, collapse=', ')))
  }

  dt[, `__BY` := paste(unlist(.BY), collapse = '.'), by = by_cols]
  dt[, `__GRP` := .GRP, by = by_cols]

  ids <- dt[, .N, by = .(`__GRP`, `__BY`)]

  grps <- ids$`__G`
  gnames <- ids$`__BY`
  dt[, `__BY` := NULL]

  glist <- lapply(grps, function(g) {
    y <- dt[`__GRP` == g, ]
    y[, `__GRP` := NULL]
    if (!dtbl) y <- as.data.frame(y)
    return(y)
  })

  names(glist) <- gnames

  return(glist)
}

pairwise <- function(n) {
  if (n < 2)
    return(NULL)
  t(utils::combn(n, 2))
}