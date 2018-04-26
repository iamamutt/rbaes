test_geom_post <- function() {
  x <- data.frame(
    cond = rep(letters[1:4], each = 2500),
    data = rnorm(10000), grp = rep(LETTERS[1:2], each = 5000))

  ggplot(x) + geom_posterior(aes(y = grp, x = data, color = cond, fill = grp),
    interval_width = 0.95, interval_type = "hdi",
    center_stat = "median") + theme_bw()
}

#' Like position_dodge but vertically
#'
#' @param height similar to \code{width} from position_dodge
#' @param preserve see position_dodge
#'
#' @return NULL
#' @export
#'
#' @examples
#' geom_point(aes(group=myGroup), position_spread(height = 0.5))
position_spread <- function(height = NULL, preserve = c("total", "single")) {
  ggproto(NULL, PositionSpread,
    height = height,
    preserve = match.arg(preserve))
}

PositionSpread <- ggproto(
  "PositionSpread",
  Position,
  required_aes = "y",
  height = NULL,
  preserve = "total",
  setup_params = function(self, data) {
    if (is.null(data$ymin) && is.null(data$ymax) &&
      is.null(self$height)) {
      warning("Height not defined. Set with `position_spread(height = ?)`",
        call. = FALSE)
    }

    if (identical(self$preserve, "total")) {
      n <- NULL
    } else {
      n <- max(table(data$ymin))
    }

    list(height = self$height, n = n)
  },
  compute_panel = function(data, params, scales) {
    collidev(data, params$height,
      name = "position_spread",
      strategy = pos_spread, n = params$n, check.width = FALSE)
  }
)

stat_posterior_density <- function(mapping = NULL, data = NULL,
                                   geom = "posterior", position = "spread", ...,
                                   bw = "nrd0", adjust = 1, kernel = "gaussian",
                                   center_stat = NULL, interval_width = NULL,
                                   interval_type = "hdi", na.rm = FALSE,
                                   show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = StatPosteriorDensity,
    geom = geom, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      bw = bw, adjust = adjust, kernel = kernel,
      center_stat = center_stat,
      interval_width = interval_width,
      interval_type = interval_type, na.rm = na.rm, ...))
}

StatPosteriorDensity <- ggproto(
  "StatPosteriorDensity",
  Stat,
  required_aes = c("x", "y"),
  non_missing_aes = "weight",
  compute_group = function(data, scales, height = NULL, bw = "nrd0",
                             adjust = 1, kernel = "gaussian",
                             center_stat = NULL, interval_width = NULL,
                             interval_type = "hdi", na.rm = FALSE) {
    if (nrow(data) < 3) {
      return(data.frame())
    }

    range <- range(data$x, na.rm = TRUE)
    bw <- calc_bw(data$x, bw)

    dens <- compute_density(
      data$x, data$w,
      from = range[1], to = range[2],
      bw = bw, adjust = adjust, kernel = kernel
    )

    dens$y <- mean(range(data$y))

    # Compute width if y has multiple values
    if (length(unique(data$y)) > 1) {
      height <- diff(range(data$y)) * 0.9
    }
    dens$height <- height

    interval_frame <- make_post_int_frame(
      data$x, center_stat, interval_width, interval_type
    )
    dens <- cbind(dens, interval_frame)
    dens
  },
  compute_panel = function(self, data, scales, height = NULL, bw = "nrd0",
                             adjust = 1, kernel = "gaussian",
                             center_stat = NULL, interval_width = NULL,
                             interval_type = "hdi", na.rm = FALSE) {
    data <- ggproto_parent(
      Stat, self
    )$compute_panel(
      data, scales,
      height = height,
      bw = bw, adjust = adjust,
      kernel = kernel,
      center_stat = center_stat,
      interval_width = interval_width,
      interval_type = interval_type
    )
    data$panelheight <- data$density / max(data$density)
    data
  }
)

#' Geom for plotting posterior distributions
#'
#' @inheritParams ggplot2::geom_violin
#' @param center_stat character string of function to compute center of distribution, such as \code{"median"} or \code{"mean"}
#' @param interval_width width of the confidence interval, e.g., 0.95
#' @param interval_type method of computing the interval, either \code{"hdi"} or \code{"ci"}
#' @export
#'
#' @examples
#' df <-
#' data.frame(
#' l = rep(letters[1:4], each = 2500),
#' v = rnorm(10000),
#' f = rep(LETTERS[1:2], each = 5000)
#' )
#' ggplot(df) +
#' geom_posterior(
#' aes(y = f, x = v, group = l),
#' interval_width = 0.95,
#' interval_type = 'hdi',
#' center_stat = 'median')
geom_posterior <- function(mapping = NULL, data = NULL,
                           stat = "posterior_density",
                           position = "spread", ..., na.rm = FALSE,
                           show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomPosterior,
    position = position, show.legend = show.legend,
    inherit.aes = inherit.aes, params = list(na.rm = na.rm, ...))
}

GeomPosterior <- ggproto(
  "GeomPosterior",
  Geom,
  setup_data = function(data, params) {
    data$height <- data$height %||% params$height %||%
      (resolution(data$y, FALSE) * 0.9)

    # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
    plyr::ddply(data, "group", transform,
      ymin = y - height / 2, ymax = y + height / 2)
  },
  draw_group = function(self, data, ...) {
    interp_dens_y <- function(xdens, ymin, ymax, x) {
      c(
        stats::approxfun(xdens, ymin)(x),
        stats::approxfun(xdens, ymax)(x))
    }

    post_trace <- trace_violin(data)
    post_fill <- post_trace

    post_trace$fill <- NA
    post_trace_grob <- GeomPolygon$draw_panel(post_trace, ...)

    post_fill$colour <- NA
    post_fill$fill <- "#FFFFFF"
    post_fill_grob <- GeomPolygon$draw_panel(post_fill, ...)

    if (!all(is.na(data$vlinectr))) {
      vlinectr <- data$vlinectr[1]
      post_vline <- data.frame(
        x = c(vlinectr, vlinectr),
        y = interp_dens_y(
          post_trace$x, post_trace$yminv,
          post_trace$ymaxv, vlinectr),
        group = c(vlinectr, vlinectr)
      )

      aesthetics <- post_trace[
        rep(1, nrow(post_vline)),
        setdiff(names(post_trace), c("x", "y", "group")),
        drop = FALSE
      ]
      aesthetics$alpha <- rep(1, nrow(post_vline))
      post_vline <- cbind(post_vline, aesthetics)
      post_vline$colour <- "#FFFFFF"
      post_vline$size <- post_vline$size * 1.25
      vline_grob <- GeomPath$draw_panel(post_vline, ...)
    } else {
      vline_grob <- NULL
    }

    if (!all(is.na(data$segxminc))) {
      post_segc <- trace_violin(
        data, data$segxminc[1], data$segxmaxc[1]
      )
      post_segc$fill <- gray(.84)
      post_segc$colour <- post_trace$colour[1]
      segc_grob <- GeomPolygon$draw_panel(post_segc, ...)
    } else {
      segc_grob <- NULL
    }

    if (!all(is.na(data$segxmins))) {
      post_segs <- trace_violin(
        data, data$segxmins[1], data$segxmaxs[1]
      )
      post_segs$colour <- "#FFFFFF"
      segs_grob <- GeomPolygon$draw_panel(post_segs, ...)
    } else {
      segs_grob <- NULL
    }

    ggname(
      "geom_posterior",
      grid::grobTree(
        post_fill_grob, segc_grob, segs_grob,
        vline_grob, post_trace_grob))
  },
  draw_key = draw_key_polygon,
  default_aes = aes(
    weight = 1, colour = "#B3B3B3", fill = "#999999",
    size = 0.5, alpha = NA, linetype = "solid"),
  required_aes = c("x", "y")
)



# subfunctions ------------------------------------------------------------


# Spread overlapping interval.
# Assumes that each set has the same vertical position.
pos_spread <- function(df, height, n = NULL) {
  if (is.null(n)) {
    n <- length(unique(df$group))
  }

  if (n == 1) {
    return(df)
  }

  if (!all(c("ymin", "ymax") %in% names(df))) {
    df$ymin <- df$y
    df$ymax <- df$y
  }

  d_height <- max(df$ymax - df$ymin)

  # Have a new group index from 1 to number of groups.
  # This might be needed if the group numbers in this set don't include all of 1:n
  groupidx <- match(df$group, sort(unique(df$group)))

  # Find the center for each group, then use that to calculate ymin and ymax
  df$y <- df$y + height * ((groupidx - 0.5) / n - 0.5)
  df$ymin <- df$y - d_height / n / 2
  df$ymax <- df$y + d_height / n / 2
  return(df)
}

# Detect and prevent collisions.  Powers dodging, stacking and filling.
collidev <- function(data, height = NULL, name, strategy, ...,
                     check.width = TRUE, reverse = FALSE) {
  # Determine width
  if (!is.null(height)) {
    # Width set manually
    if (!(all(c("ymin", "ymin") %in% names(data)))) {
      data$ymin <- data$y - height / 2
      data$ymax <- data$y + height / 2
    }
  } else {
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y
      data$ymax <- data$y
    }
    # Width determined from data, must be floating point constant
    heights <- unique(data$ymax - data$ymin)
    heights <- heights[!is.na(heights)]
    # # Suppress warning message since it's not reliable if (!zero_range(range(widths))) { warning(name, ' requires constant width: output may be incorrect', call. = FALSE) }
    height <- heights[1]
  }

  # Reorder by x position, then on group. The default stacking order reverses the group in order to match the legend order.
  if (reverse) {
    data <- data[order(data$ymin, data$group), ]
  } else {
    data <- data[order(data$ymin, -data$group), ]
  }

  # Check for overlap
  intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
  intervals <- intervals[!is.na(intervals)]
  if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < -1e-06)) {
    warning(name, " requires non-overlapping x intervals", call. = FALSE)
  }

  if (!is.null(data$xmax)) {
    return(plyr::ddply(data, "ymin", strategy, ..., height = height))
  } else {
    if (!is.null(data$x)) {
      data$xmax <- data$x
      data <- plyr::ddply(data, "ymin", strategy, ..., height = height)
      data$x <- data$xmax
      return(data)
    } else {
      stop("Neither x nor xmax defined")
    }
  }
}

calc_bw <- function(x, bw) {
  if (is.character(bw)) {
    if (length(x) < 2) {
      stop("need at least 2 points to select a bandwidth automatically",
        call. = FALSE)
    }
    bw <- switch(tolower(bw), nrd0 = stats::bw.nrd0(x),
    nrd = stats::bw.nrd(x), ucv = stats::bw.ucv(x),
    bcv = stats::bw.bcv(x), sj = ,
    `sj-ste` = stats::bw.SJ(x, method = "ste"),
    `sj-dpi` = stats::bw.SJ(x, method = "dpi"),
    stop("unknown bandwidth rule"))
  }
  bw
}

compute_density <- function(x, w, from, to, bw = "nrd0", adjust = 1,
                            kernel = "gaussian", n = 512) {
  nx <- length(x)
  if (is.null(w)) {
    w <- rep(1 / nx, nx)
  }

  # if less than 2 points return data frame of NAs and a warning
  if (nx < 2) {
    warning("Groups with fewer than two data points have been dropped.",
      call. = FALSE)
    return(data.frame(
      x = NA, density = NA, scaled = NA,
      count = NA, n = NA))
  }

  dens <- stats::density(x,
    weights = w, bw = bw, adjust = adjust,
    kernel = kernel, n = n, from = from, to = to)

  data.frame(
    x = dens$x, density = dens$y, scaled = dens$y /
      max(dens$y, na.rm = TRUE),
    count = dens$y * nx, n = nx)
}

make_post_int_frame <- function(x, center_stat = NULL, interval_width = NULL,
                                interval_type = "hdi") {
  interval_frame <- data.frame(
    vlinectr = NA_real_, segxmins = NA_real_, segxmaxs = NA_real_,
    segxminc = NA_real_, segxmaxc = NA_real_
  )

  no_central_line <- is.null(center_stat)
  no_interval_lines <- is.null(interval_width)

  if (no_central_line & no_interval_lines) {
    return(interval_frame)
  }

  if (no_central_line) {
    mid <- "median"
  } else {
    mid <- center_stat
  }

  if (no_interval_lines) {
    widths <- 0.8
  } else {
    widths <- interval_width[1]
  }

  markers <- post_int(x, mid = mid, int = interval_type, widths = widths)

  if (!no_central_line) {
    interval_frame$vlinectr <- markers$c
  }

  if (!no_interval_lines) {
    interval_frame$segxmins <- markers$l.sd
    interval_frame$segxmaxs <- markers$r.sd
    interval_frame$segxminc <- markers$l.wide
    interval_frame$segxmaxc <- markers$r.wide
  }

  interval_frame
}

trace_violin <- function(data, xmin = NULL, xmax = NULL) {
  if (!(is.null(xmin) | is.null(xmax))) {
    data <- data[data$x >= xmin & data$x <= xmax, ]
  }

  # Find the points for the line to go all the way around
  data <- transform(data,
    yminv = y - panelheight * (y - ymin),
    ymaxv = y + panelheight * (ymax - y))

  # Make sure it's sorted properly to draw the outline
  data <- rbind(
    plyr::arrange(transform(data, y = yminv), x),
    plyr::arrange(transform(data, y = ymaxv), -x))

  # Close the polygon: set first and last point the same Needed for coord_polar and such
  data <- rbind(data, data[1, ])
  return(data)
}

ggname <- function(prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}

`%||%` <- function(a, b) {
  if (!is.null(a)) {
    a
  } else {
    b
  }
}
