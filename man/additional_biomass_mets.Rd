\name{additional_biomass_mets}
\alias{additional_biomass_mets}
\docType{data}
\title{
additional_biomass_mets
}
\description{
Example list containing additional biomass metabolites. \cr\cr
met: character, specifies the metabolite, which can be added to the biomass objective function (note metabolites, which are already in the biomass objective function cannot be added).\cr
pen: numeric >0, specifies the corresponding penalty for adding this metabolite.\cr
factor: numeric -1 or 1; -1 metabolite can be added as substrate; 1 metabolite can be added as product\cr
}
\usage{data(additional_biomass_mets)}
\format{
Example:\cr
 additional_biomass_mets=list()\cr
		additional_biomass_mets[[1]]=c(met="Q[e]",pen=0.1,factor=-1)\cr
		additional_biomass_mets[[2]]=c(met="B[e]",pen=0.1,factor=-1)\cr
		additional_biomass_mets[[3]]=c(met="R[e]",pen=0.1,factor=-1)
}
\keyword{datasets}
