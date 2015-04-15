\name{additional_reactions_list}
\alias{additional_reactions_list}
\docType{data}
\title{
Example list containing additional reaction
}
\description{
Example list containing additional reaction. Each entry must contain the following attributes.\cr
		id: character, id of reaction\cr
		name: character, name of reaction\cr
		eq: character, equation of reaction.\cr
		pen: numeric >0, penalty for adding this reaction to the network\cr\cr
		
		For unidirectional reactions use "=>", for bidirectional reactions use "<=>". Compartements are assigned to metabolites with "[<compartment>]" (e.g. "Glucose[c]").\cr
		For stoichiometric coefficients use "(<numerical value>)" in front of metabolite (e.g. "(1.4) Glucose[c]"). No stoichiometric coefficients are interpreted as "(1)".
}
\usage{data(additional_reactions_list)}
\format{
Example:\cr
additional_reactions_list=list()\cr
additional_reactions_list[[1]]=list(id="KtoB",name="KtoB reaction",eq="(2.1) K[e] => (1.2) B[e]",pen=7)\cr
additional_reactions_list[[2]]=list(id="TtoR",name="TtoR reaction",eq="T[e] <=> R[e] + Q[e]",pen=3)\cr
additional_reactions_list[[3]]=list(id="TtoQ",name="TtoQ reaction",eq="T[e] => Q[e]",pen=5)\cr
}
\keyword{datasets}
