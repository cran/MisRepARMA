\encoding{utf8}
\name{fitMisRepARMA}
\alias{fitMisRepARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fit ARMA model to misreported time series data
}
\description{
  Fits an ARMA model to misreported time series data.
}
\usage{
  fitMisRepARMA(y, tol, B, alpha, p_AR, q_MA, covars=NULL, misReport="U")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  a numeric vector or time series giving the original data.
}
  \item{tol}{
  tolerance limit to stop the iterative algorithm.
}
  \item{B}{
  the number of bootstrap series to compute.
}
  \item{alpha}{
  significance level.
}
\item{p_AR}{
  order of the AR part.
}
\item{q_MA}{
  order of the MA part.
}
\item{covars}{
  matrix of explanatory variables. Its default value is \code{NULL}.
}
\item{misReport}{
  direction of misreporting issue. Its default value is \code{U} for underreported data, can also take the value \code{O} for overreported data.
}
}

\details{
The stationary bootstrap scheme with mean block length \code{b} according to Politis and Romano (1994) is computed.
}

\value{
A list of class \code{resample.statistic} with the following elements is returned:

\itemize{
  \item \code{statistic}: the results of applying statistic to each of the simulated time series.
  \item \code{orig.statistic}: the results of applying statistic to the original series.
  \item \code{bias}: the bootstrap estimate of the bias of statistic.
  \item \code{se}: the bootstrap estimate of the standard error of statistic.
  \item \code{call}: the original call of tsbootstrap.
}
}

\author{
David Moriña, Amanda Fernández-Fontelo, Alejandra Cabaña, Pedro Puig
}

\references{
D. Moriña, A. Fernández-Fontelo, A. Cabaña, P. Puig (2021): New statistical model for misreported data with application to current public health challenges. arXiv preprint (https://arxiv.org/pdf/2003.09202.pdf)

D. N. Politis and J. P. Romano (1994): The Stationary Bootstrap. Journal of the American Statistical Association \bold{89}, 1303–1313.
}

\seealso{
\code{\link{MisRepARMA-package}}, \code{\link{reconstruct}}
}

\examples{
### Simulate underreported time series data
set.seed(12345)
x <- arima.sim(model=list(ar=0.4), n=50)
ind <- rbinom(50, 1, 0.6)
y <- ifelse(ind==0, x, x*0.3)
mod <- fitMisRepARMA(y, 1e-6, 3, 0.05, 1, 0, covars=NULL, misReport="U")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MisRepARMA }
