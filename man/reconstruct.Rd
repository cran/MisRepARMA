\encoding{utf8}
\name{reconstruct}
\alias{reconstruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Reconstruct the most likely series
}
\description{
  Reconstructs the most likely series.
}
\usage{
  reconstruct(data, q, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  a numeric vector or time series giving the original data.
}
  \item{q}{
  intensity of misreporting (parameter q).
}
  \item{w}{
  frequency of misreporting (parameter w).
}
}

\value{
  the function returns a vector of the same length of \code{data} containing the reconstruction of the most likely series.
}

\references{
D. Moriña, A. Fernández-Fontelo, A. Cabaña, P. Puig (2021): New statistical model for misreported data with application to current public health challenges. arXiv preprint (https://arxiv.org/pdf/2003.09202.pdf)

D. N. Politis and J. P. Romano (1994): The Stationary Bootstrap. Journal of the American Statistical Association \bold{89}, 1303–1313.
}

\author{
David Moriña, Amanda Fernández-Fontelo, Alejandra Cabaña, Pedro Puig
}

\seealso{
\code{\link{MisRepARMA-package}}, \code{\link{fitMisRepARMA}}
}

\examples{
### Simulate underreported time series data
x <- arima.sim(model=list(ar=0.4), n=50)
ind <- rbinom(500, 1, 0.6)
y <- ifelse(ind==0, x, x*0.3)
x <- reconstruct(y, 0.6, 0.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MisRepARMA }
