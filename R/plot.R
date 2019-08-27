# This software was developed at the National Institute of Standards and Technology (NIST)
# by employees of the Federal Government in the course of their official duties. Pursuant
# to title 17 Section 105 of the United States Code, this software is not subject to
# copyright protection and is in the public domain. NIST assumes no responsibility
# whatsoever for its use by other parties, and makes no guarantees, expressed or implied,
# about its quality, reliability, or any other characteristic.

#' Returns conveniently located tick marks for a log10 plot.
#'
#' @param lim log10 of the axis limits.
#' @param loc The locations of the tick marks for each factor of 10.
#' @return The locations of the tick marks.
#' @export
log.ticks <- function(lim=c(-10, 1), loc=1)
{
   c(outer(10^seq(lim[1], lim[2]), loc))
}

#' Prints a number in human-interpretable exponential format.
#'
#' @param number The numerical value to print.
#' @return The number in scientific notation as an expression object.
#' @export
scientific.notation <- function(number)
{
   # Parse number into scientific components.
   exponent    <- floor(log10((number)))
   significand <- round(number / 10^exponent)

   # Only show significand if it's not one.
   label <- ifelse(significand == 1,
                   paste("10^", exponent),
                   paste(significand, "%*% 10^", exponent))

   parse(text=label)
}

#' Returns a scale object for continuous x or y aesthetics.
#'
#' @param expand       See scale_x_continuous() in ggplot2 package.
#' @param breaks       See scale_x_continuous() in ggplot2 package.
#' @param trans        See scale_x_continuous() in ggplot2 package.
#' @param breaks       See scale_x_continuous() in ggplot2 package.
#' @param minor_breaks See scale_x_continuous() in ggplot2 package.
#' @param labels       See scale_x_continuous() in ggplot2 package.
#' @param axis "x" for horizontal axis (default), "y" for vertical axis.
#' @export
log.scale <- function(..., expand=c(0.02, 0),
                           trans='log10',
                           breaks=log.ticks(),
                           minor_breaks=log.ticks(loc=1:9),
                           labels=scientific.notation,
                           axis="x")
{
   axis.func <- switch(axis, x=scale_x_continuous, y=scale_y_continuous, axis)
   params    <- list(expand=expand, trans=trans, breaks=breaks, labels=labels,
                     minor_breaks=minor_breaks, ...)

   do.call(axis.func, params)
}

#' Sets the default theme for DET plots.
#' @export
DET.theme <- function(panel.grid.minor  = element_line(colour="grey90", size=0.2),
                      panel.grid.major  = element_line(colour="grey80", size=0.3),
                      axis.text         = element_text(colour="darkblue", size=10),
                      axis.title        = element_text(colour="darkblue"),
                      axis.line         = element_line(colour='darkblue'),
                      legend.text       = element_text(colour='darkblue'),
                      legend.key.height = unit(0.25, "cm"),
                      legend.position   = "right",
                      legend.title      = element_blank(),
                      plot.title        = element_blank(),
                      plot.margin       = unit(c(0.2, 0.5, 0.1, 0.1), "cm"),
                      panel.border      = element_rect(colour='darkblue'),
                      ...)
{
   theme_bw() +
   theme(panel.grid.minor  = panel.grid.minor,
         panel.grid.major  = panel.grid.major,
         axis.text         = axis.text,
         axis.title        = axis.title,
         axis.line         = axis.line,
         legend.text       = legend.text,
         legend.title      = legend.title,
         legend.key.height = legend.key.height,
         legend.position   = legend.position,
         plot.title        = plot.title,
         plot.margin       = plot.margin,
         panel.border      = panel.border,
         ...)
}

#' Plots one or more DET curves.
#'
#' @param DETs A data.table representing the DET curve(s) to plot.
#' @param xlim The horizontal axis limits (default: c(8e-7, 0.2)).
#' @param ylim The vertical axis limits. "auto" (the default) selects limits that
#'             fully display all curves over the range specified by xlim.
#' @param xlabel The horizontal axis label (default: FMR)
#' @param ylabel The vertical axis label (default: FNMR)
#' @param x.axis The horizontal axis scale object.
#' @param y.axis The vertical axis scale object.
#' @param print.FNMR If nonzero, prepend legend labels with the curve's FNMR at
#'                   FMR=print.FNMR. If zero (the default), nothing is prepended.
#' @param FMR.sort Sort the curves in the legend according to their FNMR at this
#'                   FMR (default: 1e-5).
#' @param plotter The ggplot2 line drawing function (default: geom_step).
#' @param theme The ggplot2 theme to apply (default: DET.theme()).
#' @param ... Additional parameters to pass to the theme object.
#' @return A ggplot plotting object.
#' @export
#' @examples
#' gen.scores <- runif(1e5, 0, 10)
#' imp.scores <- rnorm(1e5)
#' DET <- make.DET(gen.scores, imp.scores)
#' plot.DET(DET, xlim=c(1e-4, 0.2), col="darkred", legend.position='none')
plot.DET <- function(DETs, xlim=c(1e-6, 0.2), ylim="auto",
                     xlabel=intersect(c("FMR",  "FPIR"), names(DETs)),
                     ylabel=intersect(c("FNMR", "FNIR"), names(DETs)),
                     # -GW do I need to set the name.
                     x.axis=log.scale(name=xlabel, axis="x"),
                     y.axis=log.scale(name=ylabel, axis="y", labels=identity,
                                      breaks=log.ticks(loc=c(1, 2, 5))),
                     prepend.FNMR=FALSE, FMR.sort=mean(xlim), plotter=geom_step,
                     legend=scale_colour_manual(values=colors.15),
                     # -GW key, or label
                     colour=key(DETs),
                                     # -GW is this most efficient? -GW breaks=dets[,unique(mget(key(dets)))]),
                     theme=DET.theme, legend.position="auto", ...)
{
  # Avoid altering the original object.
   DETs <- copy(DETs)
 
  # Convert 1:N column names to 1:1 if necessary.
   setnames(DETs,  c("FPIR", "FNIR"), c("FMR", "FNMR"), skip_absent=TRUE)
   
   # If vertical axis range should be automatically determined.
   if (ylim[1] == "auto")
   {
      # DETs are differentiated by the DET key.
      if (is.null(key(DETs)))
         FNMRs <- get.FNMR(DETs, FMR=xlim)
      else
         FNMRs <- DETs[,get.FNMR(.SD, FMR=xlim), by=key(DETs)]$V1
       
      breaks <- y.axis$breaks
      
      ylim <- c(max(breaks[breaks < min(FNMRs)], min(breaks)),    # lower
                min(breaks[breaks > max(FNMRs)], 1))              # upper
   }

   # If a legend is to be shown.
   if (!is.null(key(DETs)) && ! "label" %in% names(DETs))
   {
      # Assign label to each curve.
      DETs[, label := paste(.BY, collapse=", "), by=key(DETs)]

      # Get FNMR at fixed FMR for each curve.
      labels <- DETs[, .(FNMR = get.FNMR(.SD, FMR=FMR.sort)), by=label]

      # Order curves by decreasing FNMR.
      setorder(labels, -FNMR)

      # Convert label column type to sorted factor.
      DETs[,label := factor(label, labels$label, ordered=TRUE)]

      # Prepend FNMR (at fixed FMR) to legend label.
      if (prepend.FNMR)
         DETs[, label := sprintf("(%.03f) %s", labels[label == ..label]$FNMR, label)]

      if (legend.position == "auto") 
         legend.position <- "right"
   }
   else if (legend.position == "auto")
   {
      # If no key provided, don't show legend.
      legend.position <- "none"
      label <- NULL
   }

   # Convert zeros to very small values (helps with plotting).
   DETs[FMR == 0, FMR := 1e-99][FNMR == 0, FNMR := 1e-99]

   # Construct plot object and return.
   g <- ggplot(data=DETs, mapping=aes(FMR, FNMR)) +
        plotter() +
        coord_cartesian(xlim=xlim, ylim=ylim) +
        x.axis + y.axis +
        theme(legend.position=legend.position, ...) +
        legend

   # Add colour mapping if requested.
   if (!is.null(key(DETs)))
        g <- g + aes(colour=label)

   g
}

#' The Geom that draws the curve annotations, called by annotate.DET().
GeomDetAnnotate <- ggproto("GeomDetAnnotate", Geom, required_aes=c("x", "y", "label"),
                           draw_key=draw_key_blank, default_aes=aes(colour="black",
                           fill="white", size=3.88, alpha=NA, family="", fontface=1,
                           lineheight=1.2),
   draw_group   = function(self, data, params, coord, parse=FALSE, na.rm=FALSE,
                           label.padding=unit(0.25, "lines"), bgfill="white", bgalpha=1)
   {
      # Convert to data.table object.
      setDT(data)[, c("FMR", "FNMR") := list(10^x, 10^y)]

      # Axes range in plot.
      xlim <- params$x.range
      ylim <- params$y.range

      # Minimum and maximum FMR locations.
      M <- floor(max(xlim))
      m <- ceiling(min(xlim))
      
      # Start in the middle of the plot.
      FMR.start <- round(mean(c(m, M))) - 1

      # Determine horizontal location of the annotation.
      data[,FMR := 10^x][,FNMR := 10^y]
      FMR <- 10^((FMR.start + data$group[1]) %% (M - m) + m + 0.5)

      # Determine vertical of label.
      FNMR <- get.FNMR(data, FMR=FMR)

      # Estimate slope of DET curve at coordinates (given log scale).
      # The FMRs are always a factor of 10 apart, making dx == 1.
      FNMR.min   <- get.FNMR(data, FMR=FMR * sqrt(10))
      FNMR.max   <- get.FNMR(data, FMR=FMR / sqrt(10))
      FNMR.slope <- log10(FNMR.max / FNMR.min)

      # Calculate the scale ratio of the plot.
      scale.ratio <- (xlim[2] - xlim[1]) / (ylim[2] - ylim[1])

      # Aspect ratio of plotting window.
      window.ratio <- dev.size()[2] / dev.size()[1]

      # Calculate the amount of rotation to apply to the label. Depends on
      # 1) the slope of the curve,
      # 2) the span of the axes,
      # 3) the aspect ratio of the plotting window.
      angle <- -atan(window.ratio * scale.ratio * FNMR.slope) * 180 / pi

      data <- data[1, -c("x", "y", "FMR", "FNMR")]
      data[, c("x", "y", "angle") := list(log10(FMR), log10(FNMR), angle)]

      data <- coord$transform(data, params)

      # Create a grob to get the dimensions.
      grob   <- textGrob(data$label, gp=gpar(fontsize=data$size * .pt))
      width  <- grobWidth(grob)  + label.padding
      height <- grobHeight(grob) + label.padding

      # Rotate and position rectangle.
      grob.rect <- rectGrob(0.5, 0.5, width, height,
                            gp = gpar(col      = alpha(data$colour, bgalpha),
                                      fill     = alpha(bgfill, bgalpha),
                                      fontsize = data$size * .pt),
                            vp = viewport(x=data$x, y=data$y, angle=data$angle))

      grob.text <- textGrob(data$label, data$x, data$y, default.units="native",
                            gp = gpar(col        = alpha(data$colour, data$alpha),
                                     fontsize   = data$size * .pt,
                                      fontfamily = data$family,
                                      fontface   = data$fontface,
                                      lineheight = data$lineheight),
                            hjust=data$hjust, vjust=data$vjust, rot=data$angle)

      gList(grob.rect, grob.text)
   })

#' This function adds a geom to a plot which annotates DET curves with their names.
#' The annotations are spaced apart so as not to occlude one another.
#'
#' See the annotation function in ggplot2 for the meaning of the parameters.
#' @export
#' @examples
#' gen.scores <- runif(1e5, 0, 10)
#' imp.scores <- rnorm(1e5)
#' DET <- make.DET(gen.scores, imp.scores)
#' DET[, curve := "Curve Name"]
#' plot.DET(DET, xlim=c(1e-4, 0.2), col="darkred", legend.position='none') +
#' annotate.DET()
annotate.DET <- function(stat="identity", data=NULL, mapping=aes(FMR, FNMR, label=label),
                         position="identity", ..., inherit.aes=TRUE, show.legend=FALSE,
                         bgfill="white", bgalpha=1, label.padding=unit(0.25, "lines"),
                         na.rm=FALSE, parse=FALSE)
{
   layer(GeomDetAnnotate, stat, data, mapping, position,
         params=list(parse=parse, label.padding=label.padding, bgfill=bgfill,
                     bgalpha=bgalpha, na.rm=na.rm, ...),
         inherit.aes, show.legend=show.legend)
}
