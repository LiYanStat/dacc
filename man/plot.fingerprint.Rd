% Documentation for R/plot.fingerprint.R
\name{plot.fingerprint}
\alias{plot.fingerprint}
\title{Plot standardized fingerprinting residuals}
\usage{
\method{plot}{fingerprint}(
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
)
}
\arguments{
\item{x}{A fitted object returned by \code{fingerprint()}.}
\item{index}{Optional numeric vector of observation coordinates, such as time,
of the same length as the residuals. Defaults to observation index.}
\item{xlab, ylab}{Axis labels.}
\item{main}{Plot title.}
\item{pch}{Plotting symbol, default 19 (filled circles).}
\item{reference.level}{Central standard normal reference probability, default
0.95. The band is pointwise, not simultaneous.}
\item{band.col}{Fill color for the normal reference band.}
\item{ylim}{Optional vertical axis limits; by default include the residuals
and reference band.}
\item{...}{Additional graphical arguments passed to \code{graphics::plot()},
such as \code{col}, \code{cex}, or \code{type}.}
}
\description{
Displays standardized residuals from the method of Li and Li (2026) against
observation index or user-supplied coordinates, with a reference line at zero.
Patterns can indicate departures from the model.
}
\details{
Standardization uses estimated marginal variances; it does not remove dependence
between observation coordinates. Residuals need not be independent even when
the model is appropriate. The shaded band spans the standard normal quantiles
at (1 - reference.level)/2 and (1 + reference.level)/2, approximately -1.96 and
1.96 by default. It is a pointwise reference interval for residuals, not a
confidence interval for a fitted mean or a simultaneous envelope.
}
\value{The fitted object \code{x}, invisibly.}
\seealso{\code{\link{fingerprint}}}

\examples{
set.seed(12)
N <- 60
signal <- 3 * sin(seq_len(N) / 6)
ensembles <- list(sweep(matrix(rnorm(30 * N), 30), 2, signal, "+"))
Xtilde <- matrix(colMeans(ensembles[[1]]), ncol = 1)
Z <- matrix(rnorm(120 * N), nrow = 120)
Y <- signal + rnorm(N)

fit <- fingerprint(Xtilde, Y, mruns = 30, ctlruns = Z,
                   ensembles = ensembles)

# Observation index, with the default 95% normal reference band.
plot(fit)

# Use time coordinates and customize the appearance.
years <- 1961 + seq_len(N) - 1
plot(fit, index = years, xlab = "Year", col = "steelblue",
     band.col = "grey90")

# Change the pointwise reference level to 90%.
plot(fit, reference.level = 0.90)
}
