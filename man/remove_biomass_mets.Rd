\name{remove_biomass_mets}
\alias{remove_biomass_mets}
\docType{data}
\title{
Example list containing biomass metabolites, that are allowed to be removed from the biomass objective function
}
\description{
Example list containing biomass metabolites, that are allowed to be removed from the biomass objective function.\cr
Every metabolite in this list must be in the biomass objective function.\cr
The following atrributes must be defined. \cr
	met: character, name of biomass metabolite that is allowed to be removed\cr
	pen: numerical >0, according penalty for removing
}
\usage{data(remove_biomass_mets)}
\format{
  
	Example:\cr
	remove_biomass_mets=list()\cr
	remove_biomass_mets[[1]]=c(met="S[e]",pen=0.1)\cr
	remove_biomass_mets[[2]]=c(met="T[e]",pen=0.1)\cr
	remove_biomass_mets[[3]]=c(met="Z[e]",pen=0.1)\cr
}
\examples{
data(remove_biomass_mets)
}
\keyword{datasets}
