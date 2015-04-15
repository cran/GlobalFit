\name{p_list}
\alias{p_list}
\docType{data}
\title{
Example parameter list for external solver.
}
\description{
Example parameter list for external solver. In this example for cplex. Note, parameters are different for every solver.
}
\usage{data(p_list)}
\format{
Example:
p_list=list(CPX_PARAM_THREADS=as.integer(1),CPX_PARAM_EPRHS=as.double(1e-9),
CPX_PARAM_NETEPRHS=as.double(1e-11),CPX_PARAM_EPINT=as.double(1e-09),
CPX_PARAM_TILIM=1e5,CPX_PARAM_PARALLELMODE=CPX_PARALLEL_OPPORTUNISTIC)
}
\seealso{
 \code{cplexAPI}
 }
 \references{
 ILOG CPLEX Parameters Reference Manual
 }
\examples{
data(p_list)
}
\keyword{datasets}
