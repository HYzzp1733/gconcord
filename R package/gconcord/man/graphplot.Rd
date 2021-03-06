% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gconcord.R
\name{graphplot}
\alias{graphplot}
\title{Visualization of graphical model.}
\usage{
graphplot(
  met,
  varnames,
  main,
  seed,
  mode = "circle",
  label = TRUE,
  edge.width = 1,
  color = "lightpink",
  edge.color = "gray",
  edge.size = "weights",
  ...
)
}
\arguments{
\item{met}{the estimated inverse covariance matrix.}

\item{varnames}{optional, a vector containing the labels of variables. If missing, variables will be
labeled by 1, 2, 3, ...}

\item{seed}{optional, set the seed making result reproducible.}

\item{label}{optional, a logical value indicating if labels are shown in the graph. The default is \code{TRUE}.}

\item{edge.width}{optional, a positive number controling the width of edges. The default is 1.}

\item{color}{optional, color of the nodes in the graphs. See the argument \code{node.color} in \code{[GGally]{ggnet}}. The default is "lightpink".}

\item{edge.color}{optional, color of edges in the graphs. See the argument \code{node.color} in \code{[GGally]{ggnet}}. The default is "lightpink".}

\item{edge.size}{optional, size of the edges in the graphs. See the argument \code{size} in \code{[GGally]{ggnet}}. The default is "weights".}

\item{...}{optional, other arguments. See \code{[GGally]{ggnet}}.}

\item{title}{optional, title of the plot.}

\item{model}{optional, a placement method from those provided in the \code{link[GGally]{ggnet2}}
package: see \code{\link[sna]{gplot.layout}} for details. The default is "circle".}
}
\description{
Visualize the graphical model given an estimated inverse covariance matrix.
}
\examples{
p = 10
n = 50
sigma = matrix(1,nrow = p,ncol = p)*0.5 + diag(1,p,p)*0.5
data <- round(mvrnorm(n, mu = rep(0,p), sigma),2)         # generate data
train <- sample(1:nrow(data), 40, replace = FALSE)        # select training data
omega <- gconcord(data[train,], lambda1 = 0.2, lambda2 = 0.1, method = "ista")$omega
graphplot(omega)

}
