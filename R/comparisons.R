
#' table of LOO comparisons
#'
#' @param loo_list a list of objects of the type \code{loo}
#' @param stat
#'
#' @return data.table of loo comparisons
#' @export
#'
#' @examples
#' stanreg <- example_stanreg()
#' stanreg2 <- update(stanreg, . ~ . + floor:log_uranium)
#' loo_list <- lapply(list(stanreg, stanreg2), rstanarm::loo)
#' loo_table(loo_list)
loo_table <- function(loo_list, stat = "looic")
{
  n_fit <- length(loo_list)
  n_compare <- pairwise(n_fit)
  lnames <- names(loo_list)

  if (is.null(lnames)) {
    lnames <- paste0('loo', 1:n_fit)
  }

  no_names <- !nzchar(lnames)

  if (any(no_names)) {
    lnames[no_names] <- paste0('loo', which(no_names))
  }

  ldata <- lapply(loo_list, function(i)
  {
    lp <- i$pointwise[, stat]
    nlist(lp, sum_lp = sum(lp))
  })

  loo_comp <- list()

  tbl_names <- c(
    toupper(paste(stat, "L")),
    toupper(paste(stat, "R")),
    toupper(paste(stat, "Diff.")),
    "SE DIFF.",
    "CI  5%",
    "CI 95%")

  for (p in 1:nrow(n_compare)) {
    p1 <- n_compare[p, 1]
    p2 <- n_compare[p, 2]
    mod_names <- lnames[c(p1, p2)]
    lp_data <- list(ldata[[p1]]$lp, ldata[[p2]]$lp)

    if (any(unlist(lapply(lp_data, is.na)))) {
      loo_comp[[p]] <- NULL
    } else {
      sum_data <- c(ldata[[p1]]$sum_lp, ldata[[p2]]$sum_lp)
      m_order <- order(sum_data)
      loo_diff_data <- lp_data[[m_order[2]]] - lp_data[[m_order[1]]]
      loo_diff <- sum(loo_diff_data)
      loo_se <- sqrt(var(loo_diff_data) * length(loo_diff_data))
      loo_cint <- loo_ci(loo_diff, loo_se)
      tbl <- data.table::data.table(
        C = c(better = sum_data[m_order[1]],
              worse = sum_data[m_order[2]],
              diff = loo_diff,
              se = loo_se,
              cil = loo_cint[1],
              ciu = loo_cint[2])
      )
      names(tbl) <- paste(mod_names[m_order[1]], "<", mod_names[m_order[2]])
      loo_comp[[p]] <- tbl
    }
  }

  out_table <- cbind(` ` = tbl_names, do.call(cbind, loo_comp))
  return(out_table)
}

#' Get model data from stanreg
#'
#' @param stanreg object of class stanreg
#' @param get_y logical value. Include the output columns?
#'
#' @return data.table
#' @export
#' @examples
#' stanreg <- example_stanreg()
#' model_data <- stanreg_dtbl(stanreg, TRUE)
stanreg_dtbl <- function(stanreg, get_y = TRUE) {
  dt <- copy(as.data.table(stanreg$glmod$fr))

  if (!get_y) {
    y <- sub('~', '', deparse(stanreg$formula[1:2]))
    dt[, eval(y) := NULL]
  }

  return(dt)
}

#' get data for contrasts
#'
#' @param stanreg stanreg object from rstanarm
#' @param new_group a character vector of column names (factors) that should be
#'   estimated "out of sample". If \code{NA} or \code{TRUE}, then all random
#'   effects will be out of sample, i.e., new levels in those factors.
#' @param marginalize_factors a character manually choosing which factors to
#'   marginalize over (collapsing numeric values by taking the mean for each
#'   combination of factors). If \code{NA}, then the groups used as random
#'   effects will be omitted during aggregation. The default is to marginalize
#'   over all factors in the model data set.
#' @param margin_ignore character vector of numeric variable names to ingore
#'   when marginalizing over factors
#' @param set_numerics an optional named list, where each name is the name of a
#'   column in the data and its value is the value to use to override the
#'   numeric value instead of using the margin.
#' @param subset_expression an optional formula expression to subset the data
#'   before marginalizing. Expressions must begin with a \code{.~} and have
#'   conditions after this. See examples for usage details.
#'
#' @return a data.table with the same columns used in the model.
#' @export
#'
#' @examples
#' stanreg <- example_stanreg()
#' model_data <- stanreg_dtbl(stanreg) # compare below to model data
#'
#' contrast_data(stanreg)
#'
#' contrast_data(stanreg, new_group=TRUE)
#' contrast_data(stanreg, new_group=c('county'))
#'
#' contrast_data(stanreg, margin_ignore = 'floor')
#'
#' contrast_data(stanreg, TRUE, set_numerics = list(floor = 1))
#'
#' contrast_data(stanreg, subset_expression = .~ log_uranium > 0)
#' contrast_data(stanreg, margin_ignore = c('floor'), subset_expression = 1:101)
contrast_data <- function(stanreg, new_group = NULL, marginalize_factors = NULL,
                          margin_ignore = NULL, set_numerics = NULL, subset_expression = NULL)
{
  # data from rstanarm model
  dt <- stanreg_dtbl(stanreg, get_y = FALSE)

  # data classes
  dt_classes <- sapply(dt, class)
  numerics <- c('numeric', 'integer')
  factors <-  c('factor', 'character')

  nums <- names(dt_classes[dt_classes %in% numerics])
  facs <- names(dt_classes[dt_classes %in% factors])
  grps <- names(stanreg$glmod$reTrms$cnms)

  # grps can't be numeric
  nums <- nums[!nums %in% grps]

  if (length(nums) == 0) {
    nums <- NULL
  }

  if (length(facs) == 0) {
    facs <- NULL
  }

  # find new group levels
  if (length(grps) != 0 & !is.null(new_group)) {
    new_g_true <- all(is.logical(new_group)) && all(new_group)
    new_g_na <- all(is.na(new_group))
    new_g_char <- all(is.character(new_group))

    if (new_g_true | new_g_na) {
      new_group <- grps
    } else if (new_g_char) {
      new_group <- new_group[new_group %in% grps]
    } else {
      new_group <- NULL
    }

  } else {
    new_group <- NULL
  }

  # make new group levels
  if (!is.null(new_group)) {
    dt[, eval(new_group) := lapply(.SD, as.character), .SDcols = new_group]
    for (g in new_group) {
      ran_factor <- stanreg$glmod$reTrms$flist[[g]]
      new_level <- paste0("_NEW_", g)
      add_level <- c(levels(ran_factor), new_level)
      dt[, eval(g) := new_level]
      dt[, eval(g) := lapply(.SD, factor, levels=add_level), .SDcols=g]
    }
  }

  # average numerics by factors
  if (is.null(marginalize_factors)) {
    by_fac <- facs
  } else if (is.na(marginalize_factors)) {
    by_fac <- facs[!facs %in% grps]
  } else {
    by_fac <- marginalize_factors
  }

  if (!is.null(set_numerics) & !is.null(nums)) {
    num_edit <- names(set_numerics)
    num_edit <- num_edit[num_edit %in% nums]
    for (n in seq_along(num_edit)) {
      cname <- num_edit[n]
      cval <- set_numerics[[cname]]
      dt[, eval(cname) := cval]
    }
  }

  if (!is.null(margin_ignore) & !is.null(nums)) {
    nums <- nums[!nums %in% margin_ignore]
  }

  # row subsetting
  if (!is.null(subset_expression)) {
    if (is.logical(subset_expression) && any(subset_expression)) {
      index <- which(subset_expression)
    } else if (is.numeric(subset_expression) && length(subset_expression) > 0) {
      index <- subset_expression
    } else if (class(subset_expression) == 'formula') {
      index <- get_contrast_rows(subset_expression, dt)
    } else {
      stop('subset_expression must be a formula, integer vector, or logical vector.')
    }

    dt <- dt[index, ]
  }

  if (!is.null(nums)) {
    to_numeric <- sapply(dt[, .SD, .SDcols = nums], class) == 'integer'
    if (any(to_numeric)) {
      int_nums <- nums[to_numeric]
      dt[, eval(int_nums) := lapply(.SD, as.numeric), .SDcols=int_nums]
    }
    if (is.null(by_fac)) {
      dt[, eval(nums) := lapply(.SD, mean), .SDcols = nums]
    } else {
      dt[, eval(nums) := lapply(.SD, mean), by = by_fac, .SDcols = nums]
    }
  }

  dt <- unique(dt)

  if (nrow(dt) == 0) {
    dt <- NULL
  }

  return(dt)
}


#' Posterior predictive contrast
#'
#' @param stanreg stanreg object from rstanarm.
#' @param cdata a list of data sets to use to find posterior predictions of
#'   each.
#' @param ccoef contrast coefficients corresponding to the order of the data
#'   sets in \code{cdata}.
#' @param width confidence interval width. Sent to the \code{prob} argument
#'   from \code{\link[rstanarm]{posterior_interval}}.
#' @param yfun the stat function to use to collapse the predictions into a
#'   scalar value. Corresponds to the average of the predicted responses for each
#'   data set in \code{cdata}. Default is median.
#' @param ... additional parameters passed to the function \code{\link[rstanarm]{posterior_predict}}
#'   \code{\link[rstanarm]{posterior_interval}}.
#'
#' @return list containing the confidence interval and posterior predictive
#'   contrast difference
#' @export
#'
#' @examples
#' stanreg <- example_stanreg()
#' cdata <- list(
#'  c1 = contrast_data(stanreg, TRUE, margin_ignore = 'floor', subset_expression = .~ log_uranium < -0.5),
#'  c2 = contrast_data(stanreg, TRUE, margin_ignore = 'floor', subset_expression = .~ log_uranium >= -0.5 & log_uranium <= 0.5),
#'  c3 = contrast_data(stanreg, TRUE, margin_ignore = 'floor', subset_expression = .~ log_uranium > 0.5))
#'
#' my_contrast <- pp_contrast(stanreg, cdata, c(-1, 0.5, 0.5))
#' my_contrast$ci
pp_contrast <- function
(
  stanreg,
  cdata,
  ccoef,
  width = 0.95,
  yfun = median,
  ...
)
{
  # posterior predict
  pp_list <- pp_cdata(
    stanreg,
    cdata = cdata,
    row_stat = yfun,
    ...)

  # posterior predictive contrast
  pp <- apply_contrast(pp_list, ccoef)

  # contrast interval
  ci <- rstanarm::posterior_interval(pp, prob = width)

  nlist(ci, pp)
}

#' posterior predictive contrasts for all pairwise comparisons in a list
#'
#' @inheritParams pp_contrast
#' @param pp_all posterior predictions for each list item in cdata. Defaults to
#'   \code{FALSE}
#' @return list containing CI matrix and posterior differences
#' @export
#'
#' @examples
#' stanreg <- example_stanreg()
#' cdata <- list(
#'  c1 = contrast_data(stanreg, TRUE, margin_ignore = 'floor', subset_expression = .~ log_uranium < -0.5),
#'  c2 = contrast_data(stanreg, TRUE, margin_ignore = 'floor', subset_expression = .~ log_uranium >= -0.5 & log_uranium <= 0.5),
#'  c3 = contrast_data(stanreg, TRUE, margin_ignore = 'floor', subset_expression = .~ log_uranium > 0.5))
#' pairwise_contrasts(stanreg, cdata)$ci.cont
pairwise_contrasts <- function(stanreg, cdata, width=0.95, yfun=median, pp_all=FALSE, ...)
{
  pp_list <- pp_cdata(
    stanreg,
    cdata = cdata,
    row_stat = yfun,
    ...)

  n_datasets <- length(pp_list)
  cont_pairs <- pairwise(n_datasets)
  n_pairs <- nrow(cont_pairs)
  data_names <- names(pp_list)
  cont_names <- paste(
    data_names[cont_pairs[, 1]], '-',
    data_names[cont_pairs[, 2]])

  # get posterior predictive contrasts
  diff_list <- lapply(seq_len(n_pairs), function(p) {
    apply_contrast(
      pp_list = list(
        pp_list[[cont_pairs[p, 1]]],
        pp_list[[cont_pairs[p, 2]]]),
      ccoef = c(1, -1)
    )
  })

  # posterior contrasts matrix
  pred.cont <- do.call(cbind, diff_list)
  colnames(pred.cont) <- cont_names

  # get confidence intervals
  ci.cont <- rstanarm::posterior_interval(pred.cont, prob = width)

  if (pp_all) {
    pred.each <- do.call(cbind, lapply(pp_list, function(p) {
      apply_contrast(list(p), 1)
    }))
    colnames(pred.each) <- names(pp_list)
    ci.each <- rstanarm::posterior_interval(pred.each, prob = width)
  } else {
    pred.each <- NULL
    ci.each <- NULL
  }

  nlist(ci.cont, pred.cont, ci.each, pred.each)
}

# subfunctions -----------------------------------------------------------

loo_ci <- function(x, s, p = 0.95) {
  se <- qnorm(p) * s
  ci <- c(x - se, x + se)
  names(ci) <- c("05%", "95%")
  return(ci)
}

get_contrast_rows <- function(form, data) {
  if (nchar(deparse(form[2])) != 3) {
    stop('improper contrast formula specification')
  }

  data[, eval(form[-1])][[1]]
}

apply_contrast <- function(pp_list, ccoef) {
  n_pp <- length(pp_list)

  if (length(ccoef) != n_pp) {
    stop('length coef must match length pp_list')
  }

  for (pp in seq_len(n_pp)) {
    pp_list[[pp]] <- pp_list[[pp]] * ccoef[pp]
  }

  as.matrix(apply(do.call(cbind, pp_list), 1, sum))
}

pp_cdata <- function(stanreg, cdata, row_stat = median, ...) {

  if (!'stanreg' %in% class(stanreg)) {
    stop('"stanreg" arg must be a stanreg object')
  }

  if (!is.list(cdata)) {
    stop('"cdata" must be a list of datasets to compare')
  }

  n_datasets <- length(cdata)

  if (n_datasets < 2) {
    stop('need at least 2 contrast data sets')
  }

  cdata <- lapply(cdata, function(d) {
    if (!is.data.table(d)) {
      d <- as.data.table(d)
    }
    return(d)
  })

  cnames <- names(cdata)

  if (is.null(cnames)) {
    cnames <- paste('contrast', 1:n_datasets)
  }

  no_names <- !nzchar(cnames)

  if (any(no_names)) {
    cnames[no_names] <- paste('contrast', which(no_names))
  }

  # combine all the data to use the same draws
  cdata_bind <- rbindlist(
    cdata, use.names = TRUE, fill = TRUE,
    idcol="__cont_idx")

  posterior <- rstanarm::posterior_predict(stanreg, newdata=cdata_bind, ...)

  # get indices of which cols to pull from posterior matrix
  cont_grps <- dtbl2list(cdata_bind[, .I, by = '__cont_idx'], `__cont_idx`)

  # mean of rows at cols x
  pp <- lapply(cont_grps, function(i) {
    cidx <- i$I
    y <- posterior[, cidx]
    if (!is.null(dim(y))) {
      y <- apply(y, 1, row_stat)
    }

    matrix(y)

  })

  return(pp)
}

