% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/re_functions2.r
\name{plot.origin.results}
\alias{plot.origin.results}
\title{plots the output of find.origin as a heatmap using .filled.contour}
\usage{
plot.origin.results(x, n.levels = 100, color.function = heat.colors,
  color.negative = "grey", add.map = T, add.samples = T,
  add.sample.het = T, add.likely.origin = T, asp = 1, ...)
}
\arguments{
\item{x}{an object of type origin.results, as obtained by}

\item{n.levels}{the number of color levels}

\item{color.function}{a function that takes an integer argument and
returns that many colors}

\item{color.negative}{a single color to be used for negative values}

\item{add.map}{boolean whether a map should be added}

\item{add.likely.origin}{boolean, whether origin should be marked with an}

\item{asp}{aspect ratio, set to 1 to keep aspect ratio with plot
X}
}
\description{
plots the output of find.origin as a heatmap using .filled.contour
}

