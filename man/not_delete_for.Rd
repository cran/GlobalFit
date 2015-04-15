\name{not_delete_for}
\alias{not_delete_for}
\docType{data}
\title{
Example vector that specifies reactions that are not allowed to be removed in forward direction. 
}
\description{
Example vector that specifies reactions that are not allowed to be removed in forward direction. Usually exchange reaction, biomass objective function, ...
	
}
\usage{data(not_delete_for)}
\format{
Example:
not_delete_for=c("T[e]_import","S[e]_import","Biomass")
}
\examples{
data(not_delete_for)
## maybe str(not_delete_for) ; plot(not_delete_for) ...
}
\keyword{datasets}
