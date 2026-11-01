##' Plot standardized fingerprinting residuals
##'
##' Displays standardized residuals from the method of Li and Li (2026)
##' against observation index or user-supplied coordinates, with a horizontal
##' reference line at zero. Patterns can indicate departures from the model.
##' @param x A fitted object returned by \code{fingerprint()}.
##' @param index Optional numeric vector of observation coordinates, such as
##'   time, of the same length as the residuals. Defaults to observation index.
##' @param xlab,ylab Axis labels.
##' @param main Plot title.
##' @param pch Plotting symbol, default 19 (filled circles).
##' @param reference.level Central standard normal reference probability, default
##'   0.95. The band is pointwise, not simultaneous.
##' @param band.col Fill color for the normal reference band.
##' @param ylim Optional vertical axis limits; by default include the residuals
##'   and the reference band.
##' @param ... Additional graphical arguments passed to \code{graphics::plot()},
##'   such as \code{col}, \code{cex}, or \code{type}.
##' @return The fitted object \code{x}, invisibly.
##' @details Standardization uses estimated marginal variances; it does not
##'   remove dependence between observation coordinates. The residuals therefore
##'   need not be independent, even when the model is appropriate. The shaded
##'   band spans standard normal quantiles at (1 - reference.level)/2 and
##'   (1 + reference.level)/2. It is a pointwise reference interval for residuals,
##'   not a confidence interval for a fitted mean or a simultaneous envelope.
##' @seealso \code{\link{fingerprint}}
##' @examples
##' set.seed(12)
##' N <- 60
##' signal <- 3 * sin(seq_len(N) / 6)
##' ensembles <- list(sweep(matrix(rnorm(30 * N), 30), 2, signal, "+"))
##' Xtilde <- matrix(colMeans(ensembles[[1]]), ncol = 1)
##' Z <- matrix(rnorm(120 * N), nrow = 120)
##' Y <- signal + rnorm(N)
##'
##' fit <- fingerprint(Xtilde, Y, mruns = 30, ctlruns = Z,
##'                    ensembles = ensembles)
##'
##' # Observation index, with the default 95% normal reference band.
##' plot(fit)
##'
##' # Use time coordinates and customize the appearance.
##' years <- 1961 + seq_len(N) - 1
##' plot(fit, index = years, xlab = "Year", col = "steelblue",
##'      band.col = "grey90")
##'
##' # Change the pointwise reference level to 90%.
##' plot(fit, reference.level = 0.90)
##' @export
plot.fingerprint <- function(
    x,
    index = NULL,
    xlab = if (is.null(index)) "Observation index" else "Coordinate",
    ylab = "Standardized residual",
    main = "Standardized residuals",
    pch = 19,
    reference.level = 0.95,
    band.col = "grey92",
    ylim = NULL,
    ...
) {
    # Evaluate the default label before replacing a NULL index.
    xlab <- xlab
    residuals <- x$standardized.residuals

    if (is.null(index)) {
        index <- seq_along(residuals)
    }

    if (
        !is.numeric(index) || length(index) != length(residuals) ||
        any(!is.finite(index))
    ) {
        stop("index must be a finite numeric vector with one value per residual")
    }

    if (
        !is.numeric(reference.level) || length(reference.level) != 1L ||
        !is.finite(reference.level) ||
        reference.level <= 0 || reference.level >= 1
    ) {
        stop("reference.level must be a number strictly between zero and one")
    }

    band <- stats::qnorm(c(
        (1 - reference.level) / 2,
        (1 + reference.level) / 2
    ))
    if (is.null(ylim)) {
        ylim <- range(c(residuals, band), finite = TRUE)
    }

    graphics::plot(
        index,
        residuals,
        xlab = xlab,
        ylab = ylab,
        main = main,
        pch = pch,
        ylim = ylim,
        panel.first = {
            # Draw the reference band behind the residuals.
            bounds <- graphics::par("usr")
            graphics::rect(
                bounds[1], band[1], bounds[2], band[2],
                col = band.col, border = NA
            )
            graphics::abline(h = band, col = "grey65", lty = 3)
        },
        ...
    )
    graphics::abline(h = 0, col = "grey50", lty = 2)

    invisible(x)
}
