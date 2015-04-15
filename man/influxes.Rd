\name{influxes}
\alias{influxes}
\docType{data}
\title{
Example list of exchange reactions and corresponding values
}
\description{
List of exchange fluxes. The following attributes must be defined for each entry:\cr
exRea: character, name of exchange reaction\cr
value: numerical <0, value of influx\cr
}
\usage{data(influxes)}
\format{
  influxes=list()\cr
	influxes[[1]]=list(exRea="T[e]_import",value=-10)\cr
}
\examples{
data(influxes)
## maybe str(influxes) ; plot(influxes) ...
}
\keyword{datasets}
