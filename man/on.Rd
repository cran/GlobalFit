\name{on}
\alias{on}
\docType{data}
\title{
Example list of growth cases.
}
\description{
Example list of growth cases. Each entry must have defined the following attributes: \cr
	on: list of exchange reactions, with according negative value.\cr
	name: character, name of non-growth case\cr
	ko_react: vector of reactions names. Upper and lower bound of these reaction will be set to zero\cr
	forced: logical value, if non-growth cases is allowed to be ignored (FALSE), or if value of biomass objective function must be lower than viability threshold (TRUE)\cr
	viability threshold: numerical >0, value that is considered as growth or non-growth\cr
	gene-copy number: integer >0, value that will be multiplied with case_ignore_penalty\cr
	
	
}
\usage{data(on)}
\format{
  Example: \cr
	off=list()\cr
	on[[1]]=list(on=influxes,name="LIVE!",ko_react=c("TtoZ"),forced=TRUE,viability_threshold=1,gene_copy_number=1)\cr
}
\examples{
data(on)
}
\keyword{datasets}
