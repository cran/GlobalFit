\name{reverse_reaction_list}
\alias{reverse_reaction_list}
\docType{data}
\title{
Example list of reaction that are allowed to be removed.
}
\description{
Example list of reaction that are allowed to be removed. The following attributes must be defined for each entry.\cr

	reaction: character, id of reaction that is allowed to be reversed\cr
	pen: numerical >0, according penalty for reversing\cr
}
\usage{data(reverse_reaction_list)}
\format{
Example:
reverse_reaction_list=list()\cr
reverse_reaction_list[[1]]=list(reaction="KtoT",pen=1)\cr
reverse_reaction_list[[2]]=list(reaction="TtoB",pen=1)\cr
}
\examples{
data(reverse_reaction_list)
}
\keyword{datasets}
