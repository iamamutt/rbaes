
#' all pairwise combinations
pairwise <- function(n) {
  if (n < 2) return(NULL)
  t(utils::combn(n, 2))
}


#' named list
nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))

  if (no_names) {
    has_name <- FALSE
  } else {
    has_name <- nzchar(names(out))
  }

  if (all(has_name)) return(out)

  nms <- as.character(m)[-1L]

  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }
  return(out)
}

#' print text section header
print_sec <- function(x, char=79) {
  if (missing(x)) x <- ""
  nt <- nchar(x)
  if (nt >= char) char <- nt+16
  cat(c("\n\n", rep("-", 16), x, rep("-", (char-16)-nt), "\n\n"), sep="")
}

#' convert dots to list
dots2list <- function(...) {
  eval(substitute(alist(...)))
}

#' convert symbol/equation to character
symbol2char <- function(...) {
  lapply(dots2list(...), deparse)
}

#' extract items from deep within a named list
lextract <- function(x, ...) {
  entries <- symbol2char(...)

  get_from_list <- function(l, n) {
    if (!l %?n% n)
      return(NULL)
    l[[n]]
  }

  lapply(x, function(i) {
    Reduce(get_from_list,
           entries,
           init = i,
           accumulate = FALSE)
  })
}

`%?n%` <- function(obj, names) {
  obj_names <- names(obj)
  if (all(names %in% obj_names)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
