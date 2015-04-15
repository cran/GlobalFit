\name{not_delete_back}
\alias{not_delete_back}
\docType{data}
\title{
Example vector that specifies reactions that are not allowed to be removed in backward direction. 
}
\description{
Example vector that specifies reactions that are not allowed to be removed in backward direction. Usually exchange reaction, biomass objective function, ...
	
}
\usage{data(not_delete_back)}
\format{
Example:
not_delete_back=c("T[e]_import","S[e]_import","Biomass")
}
\examples{
data(not_delete_back)
## maybe str(not_delete_back) ; plot(not_delete_back) ...
}
\keyword{datasets}
