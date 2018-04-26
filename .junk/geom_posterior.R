library(ggplot2)
test <- function() {


  x <- data.frame(l=rep(letters[1:2], each=5000), v=rnorm(10000))
  ggplot(x) +
    geom_violinh(aes(y=l, x=v), draw_quantiles = c(.25, .5, .75))

}



#' PositionDodge but vertically
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
          preserve = match.arg(preserve)
  )
}

#' @export
PositionSpread <- ggproto(
  "PositionSpread", Position, required_aes = "y", height = NULL, preserve = "total",
  setup_params = function(self, data) {
    if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
      warning("Height not defined. Set with `position_spread(height = ?)`", call. = FALSE)
    }

    if (identical(self$preserve, "total")) {
      n <- NULL
    } else {
      n <- max(table(data$ymin))
    }

    list(height = self$height, n = n)
  },
  compute_panel = function(data, params, scales) {
    collidev(data, params$height, name = "position_spread",
             strategy = pos_spread, n = params$n, check.width = FALSE)
  }
)

# Spread overlapping interval.
# Assumes that each set has the same vertical position.
pos_spread <- function(df, height, n = NULL) {
  if (is.null(n)) {
    n <- length(unique(df$group))
  }

  if (n == 1)
    return(df)

  if (!all(c("ymin", "ymax") %in% names(df))) {
    df$ymin <- df$y
    df$ymax <- df$y
  }

  d_height <- max(df$ymax - df$ymin)

  # Have a new group index from 1 to number of groups.
  # This might be needed if the group numbers in this set don't include all of 1:n
  groupidx <- match(df$group, sort(unique(df$group)))

  # Find the center for each group, then use that to calculate xmin and xmax
  df$y <- df$y + height * ((groupidx - 0.5)/n - 0.5)
  df$ymin <- df$y - d_height/n/2
  df$ymax <- df$y + d_height/n/2
  return(df)
}

# Detect and prevent collisions.  Powers dodging, stacking and filling.
collidev <- function(data, height = NULL, name, strategy, ...,
                     check.width = TRUE, reverse = FALSE) {
  # Determine width
  if (!is.null(height)) {
    # Width set manually
    if (!(all(c("ymin", "ymin") %in% names(data)))) {
      data$ymin <- data$y - height/2
      data$ymax <- data$y + height/2
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
  } else if (!is.null(data$x)) {
    data$xmax <- data$x
    data <- plyr::ddply(data, "ymin", strategy, ..., height = height)
    data$x <- data$xmax
    return(data)
  } else {
    stop("Neither x nor xmax defined")
  }
}


`%||%` <- function (a, b)
{
  if (!is.null(a))
    a
  else b
}

generate <- function(fn) {
  get(fn, -1, environment(ggplot_build))
}

ggname <- function(prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}


#' Horizontal violin plot.
#'
#' Horizontal version of \code{\link[ggplot2]{geom_violin}()}.
#' @inheritParams ggplot2::geom_violin
#' @inheritParams ggplot2::geom_point
#' @export
#' @examples
#' d <- do.call(rbind, lapply(1:3, function(i) data.frame(g=letters[i], v=rnorm(1000, 0, i+1))))
#' ggplot(d, aes(y=g))+geom_violinh(aes(x=v))
geom_violinh <-
  function(mapping = NULL,
           data = NULL,
           stat = "xdensity",
           position = "spread",
           ...,
           draw_quantiles = NULL,
           trim = TRUE,
           scale = "area",
           na.rm = FALSE,
           show.legend = NA,
           inherit.aes = TRUE) {
    layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomViolinh,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        trim = trim,
        scale = scale,
        draw_quantiles = draw_quantiles,
        na.rm = na.rm,
        ...
      )
    )
  }

#' @rdname ggstance-ggproto
#' @format NULL
#' @usage NULL
#' @include legend-draw.R
#' @export
GeomViolinh <-
  ggproto(
    "GeomViolinh",
    Geom,
    setup_data = function(data, params) {
      data$width <- data$width %||% params$width %||% (resolution(data$y, FALSE) * 0.9)

      # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      plyr::ddply(data, "group", transform, ymin = y - width/2, ymax = y + width/2)
    },
    draw_group = function(self, data, ..., draw_quantiles = NULL) {

      # Find the points for the line to go all the way around
      data <- transform(data, yminv = y - violinwidth * (y - ymin), ymaxv = y + violinwidth * (ymax - y))

      # Make sure it's sorted properly to draw the outline
      newdata <- rbind(plyr::arrange(transform(data, y = yminv), x), plyr::arrange(transform(data, y = ymaxv), -x))

      # Close the polygon: set first and last point the same Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1,])

      # Draw quantiles if requested, so long as there is non-zero y range
      if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$x))) {
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))

        # Compute the quantile segments and combine with existing aesthetics
        quantiles <- create_quantile_segment_frame(data, draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- GeomPath$draw_panel(both, ...)

        browser()
        ggname("geom_violin", grid::grobTree(
          GeomPolygon$draw_panel(newdata, ...),
          quantile_grob))

      } else {
        ggname("geom_violin", GeomPolygon$draw_panel(newdata, ...))
      }
    },
    draw_key = draw_key_polygon,
    default_aes = aes(
      weight = 1,
      colour = "grey20",
      fill = "white",
      size = 0.5,
      alpha = NA,
      linetype = "solid"
    ),
    required_aes = c("x", "y")
  )


#' Density computation on x axis.
#'
#' Horizontal version of \code{\link[ggplot2]{stat_ydensity}}().
#' @inheritParams ggplot2::stat_ydensity
#' @export
stat_xdensity <- function(mapping = NULL,
                          data = NULL,
                          geom = "violinh",
                          position = "dodgev",
                          ...,
                          bw = "nrd0",
                          adjust = 1,
                          kernel = "gaussian",
                          trim = TRUE,
                          scale = "area",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = TRUE) {
  scale <- match.arg(scale, c("area", "count", "width"))

  layer(
    data = data,
    mapping = mapping,
    stat = StatXdensity,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      bw = bw,
      adjust = adjust,
      kernel = kernel,
      trim = trim,
      scale = scale,
      na.rm = na.rm,
      ...
    )
  )
}

calc_bw <- generate("calc_bw")
compute_density <- generate("compute_density")


#' @rdname ggstance-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatXdensity <- ggproto(
  "StatXdensity",
  Stat,
  required_aes = c("x", "y"),
  non_missing_aes = "weight",
  compute_group = function(data,
                           scales,
                           width = NULL,
                           bw = "nrd0",
                           adjust = 1,
                           kernel = "gaussian",
                           trim = TRUE,
                           na.rm = FALSE) {
    if (nrow(data) < 3)
      return(data.frame())
    range <- range(data$x, na.rm = TRUE)
    modifier <- if (trim) {
      0
    } else {
      3
    }

    bw <- calc_bw(data$x, bw)

    dens <-
      compute_density(
        data$x,
        data$w,
        from = range[1] - modifier * bw,
        to = range[2] + modifier * bw,
        bw = bw,
        adjust = adjust,
        kernel = kernel
      )

    # dens$y <- dens$x
    dens$y <- mean(range(data$y))

    # Compute width if y has multiple values
    if (length(unique(data$y)) > 1) {
      width <- diff(range(data$y)) * 0.9
    }
    dens$width <- width

    dens
  },

  compute_panel = function(self,
                           data,
                           scales,
                           width = NULL,
                           bw = "nrd0",
                           adjust = 1,
                           kernel = "gaussian",
                           trim = TRUE,
                           na.rm = FALSE,
                           scale = "area") {
    data <- ggproto_parent(Stat, self)$compute_panel(
      data,
      scales,
      width = width,
      bw = bw,
      adjust = adjust,
      kernel = kernel,
      trim = trim,
      na.rm = na.rm
    )

    # choose how violins are scaled relative to each other
    data$violinwidth <- switch(
      scale,
      # area : keep the original densities but scale them to a max width of 1
      #        for plotting purposes only
      area = data$density / max(data$density),
      # count: use the original densities scaled to a maximum of 1 (as above)
      #        and then scale them according to the number of observations
      count = data$density / max(data$density) * data$n / max(data$n),
      # width: constant width (density scaled to a maximum of 1)
      width = data$scaled
    )
    data
  }

)

# Interleave (or zip) multiple units into one vector
interleave <- function(...) UseMethod("interleave")
#' @export
interleave.unit <- function(...) {
  do.call("unit.c", do.call("interleave.default", plyr::llply(list(...), as.list)))
}
#' @export
interleave.default <- function(...) {
  vectors <- list(...)

  # Check lengths
  lengths <- unique(setdiff(plyr::laply(vectors, length), 1))
  if (length(lengths) == 0) lengths <- 1
  stopifnot(length(lengths) <= 1)

  # Replicate elements of length one up to correct length
  singletons <- plyr::laply(vectors, length) == 1
  vectors[singletons] <- plyr::llply(vectors[singletons], rep, lengths)

  # Interleave vectors
  n <- lengths
  p <- length(vectors)
  interleave <- rep(1:n, each = p) + seq(0, p - 1) * n
  unlist(vectors, recursive = FALSE)[interleave]
}

posterior_slices <- function(data, quantiles) {

}

# Returns a data.frame with info needed to draw quantile segments.
create_quantile_segment_frame <- function(data, draw_quantiles) {

  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$x)
  xs <- ecdf(draw_quantiles) # these are all the x-values for quantiles

  # Get the violin bounds for the requested quantiles.
  violin.yminvs <- (stats::approxfun(data$x, data$yminv))(xs)
  violin.ymaxvs <- (stats::approxfun(data$x, data$ymaxv))(xs)
  browser()
  # We have two rows per segment drawn. Each segment gets its own group.
  data.frame(
    x = rep(xs, each = 2),
    y = interleave(violin.yminvs, violin.ymaxvs),
    group = rep(xs, each = 2)
  )
}
