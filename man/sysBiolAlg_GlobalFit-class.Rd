\name{sysBiolAlg_GlobalFit-class}
\Rdversion{1.1}
\docType{class}
\alias{sysBiolAlg_GlobalFit-class}

\title{Class \code{"sysBiolAlg_GlobalFit"}}
\description{
Class, that builds bilevel optimization object.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sysBiolAlg_GlobalFit", model,
LPvariant, useNames, cnames, rnames, pname, scaling, penalties_hin, penalties_ruck,
cancel_case_penalty, num_additional_reactions, not_delete_for, not_delete_back, off, on,
simple, reverse_hin, reverse_hin_penalty, reverse_back, reverse_back_penalty, MaxPenalty,
alternatives_list, additional_biomass_metabolites, remove_biomass_metabolites, bio_stoich,
use_indicator_constraints, writeProbToFileName, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{wu}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{wl}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{fnc}:}{Object of class \code{"integer"} ~~ }
    \item{\code{fnr}:}{Object of class \code{"integer"} ~~ }
    \item{\code{problem}:}{Object of class \code{"optObj"} ~~ }
    \item{\code{algorithm}:}{Object of class \code{"character"} ~~ }
    \item{\code{nr}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nc}:}{Object of class \code{"integer"} ~~ }
    \item{\code{fldind}:}{Object of class \code{"integer"} ~~ }
    \item{\code{alg_par}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{sysBiolAlg}"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "sysBiolAlg_GlobalFit")}: ... }
	 }
}
\author{
Daniel Hartleb
}
\examples{
showClass("sysBiolAlg_GlobalFit")
}
\keyword{classes}
