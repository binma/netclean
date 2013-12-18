\name{deconreg}
\alias{deconreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
deconreg(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

# network deconvolution algorithm for gene regulatory networks of DREAM5
# challenge.
# Usage: mat_nd=ND_regulatory(mat,beta,alpha,control_p)
# INPUT ARGUMENTS:
  # mat           Input matrix, it is a n_tf by n matrix where first n_tf genes are TFs.
#               Elements of the input matrix are nonnegative.
# optional parameters:
  # beta          Scaling parameter, the program maps the largest absolute eigenvalue
#               of the direct dependency matrix to beta. It should be
#               between 0 and 1. You should skip this scaling step if you know eigenvalues
#               of your matrix satisfy ND conditions.
# alpha         fraction of edges of the observed dependency matrix to be kept in
#               deconvolution process.
# control_p     If set to one, it perturbs input networks slightly to have
#               stable results in case of non-diagonalizable matrices.
#               If zero, it checks some sufficient condition and then add a small perturbation
#               (this may be slower). Default is zero.

#
# OUTPUT ARGUMENTS:
  
  # mat_nd        Output deconvolved matrix (direct dependency matrix). Its components
#               represent direct edge weights of observed interactions.
#               Choosing top direct interactions (a cut-off) depends on the application and
#               is not implemented in this code.

# In this implementation, input matrices are made symmetric.


# LICENSE: MIT-KELLIS LAB


# REFERENCES:
  #    For more details, see the following paper:
  #    Network Deconvolution as a General Method to Distinguish
#    Direct Dependencies over Networks
#    By: Soheil Feizi, Daniel Marbach,  Muriel M茅dard and Manolis Kellis
#    Nature Biotechnology
#

