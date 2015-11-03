\name{gpsrBasic}
\alias{gpsrBasic}
\alias{gpsrBB}
\title{
Functions implementing the GPSR algorithm
}
\description{
These functions solve the convex problem :
arg min_x = 0.5*|| y - A x ||_2^2 + tau || x ||_1
subject to the additional constraint x >= 0
}
\usage{
gpsrBasic(y, A, dimx, tauVec, cont_steps = 3, initAty = TRUE,
	tolA = 0.0001, maxiter = 10000)
	
gpsrBB(y, A, dimx, tauVec, cont_steps = 4, initAty = TRUE, 
	tolA = 0.0001, maxiter = 10000)
}
\arguments{
  \item{y}{
A real valued vector y containing the observations.
}
  \item{A}{
A real valued matrix A of size \code{dimx} by len(\code{y}). 
}
  \item{dimx}{
The length of hidden vector \code{x}.
}
  \item{tauVec}{
A regularization (i.e. penalty) vector with the same size as \code{x}.
}
  \item{cont_steps}{
Number of continuation steps. Default value is 3.
}
  \item{initAty}{
Logical value to decide whether to initialize x to t(A)y. If false, x is initialized to 0. True by default.
}
  \item{tolA}{
A tolerance threshold to stop iterations. If the difference between the objective scores in successive iterations is less than this value, the algorithm terminates. Default value is 0.0001.
}
  \item{maxiter}{
Maximum number of iterations to execute before terminating. Default value is 10000.
}
}
\value{
Both functions return a list with the following elements:
\item{solFound}{A logical value indicating if a feasible solution is found}
\item{solMatrix}{A matrix of size \code{dimx} by \code{cont_steps}}
}
\references{
These algorithms are based on the GPSR-Basic and GPSR-Barzilai-borwein algorithms given in:
"Gradient Projection for Sparse Reconstruction: Application
to Compressed Sensing and Other Inverse Problems"
by Mario A. T. Figueiredo, Robert D. Nowak, Stephen J. Wright
}
\author{
Nilgun Donmez
}
\note{
Part of compas package. More information can be found at:

https://github.com/nlgndnmz/COMPAS.git
}
\seealso{
See also \code{\link{compasCore}}
}
