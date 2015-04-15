\name{bilevel_optimize}
\alias{bilevel_optimize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pre und post-processing of optimization
}
\description{
pre und post-processing of optimization. Adds reaction to network, creates output, creates optimzation object, interprets solution, applies modification to network.
}
\usage{
bilevel_optimize(network, on = c(), off = c(), additional_reactions = NULL, 
minimize = TRUE, simple = FALSE, verboseMode = 1, cancel_case_penalty = NULL, 
not_delete_for = c(), not_delete_back = c(), param_list = NULL, 
use_indicator_constraints = FALSE, stat_file = NULL, react_file = NULL, 
reverse_reaction_list = NULL, MaxPenalty = NULL, alternatives = 0, 
bio_stoich=1e-5, additional_biomass_metabolites = NULL, 
remove_biomass_metabolites = NULL)
}
\arguments{
  \item{network}{
metabolic network model of type \code{modelorg}
see also: \code{\link{test}}
}
  \item{on}{
list of growth cases
	each entry must contain:
		on: list of influxes\cr
		name: character, name of growth case\cr
		ko_react: vector of reaction, which are knocked out (i.e.,lower und upper bound = 0)\cr
		forced: logical, FALSE (growth case can be ignored with according penalty) or TRUE (case cannot be ignored)\cr
		viability_threshold: numerical>0 threshold, which is considered as growth\cr
		gene_copy_number= integer >0, multiplies the penalty for ignoring this growth case\cr
		
		default: NULL\cr
		see also: \code{\link{on}}
}
  \item{off}{
list of non-growth cases
	each entry must contain:
		on: list of influxes\cr
		name: character, name of non-growth case\cr
		ko_react: vector of reaction, which are knocked out (i.e.,lower und upper bound = 0)\cr
		forced: logical, FALSE (growth case can be ignored with according penalty) or TRUE (case cannot be ignored)\cr
		viability_threshold: numerical>0 threshold, which is considered as growth\cr
		gene_copy_number= integer >0, multiplies the penalty for ignoring this non-growth case\cr
		
		default: NULL\cr
		see also: \code{\link{off}}
}
  \item{additional_reactions}{
	list containing additional reaction. Each entry must contain the following attributes.\cr
		id: character, id of reaction\cr
		name: character, name of reaction \cr
		eq: character, equation of reaction.\cr
		pen: numeric >0, penalty for adding this reaction to the network\cr
		
		Example:\cr
		additional_reactions_list=list()\cr
		additional_reactions_list[[1]]=list(id="KtoB",name="KtoB reaction",eq="(2.1) K[e] => (1.2) B[e]",pen=7)\cr
		additional_reactions_list[[2]]=list(id="TtoR",name="TtoR reaction",eq="T[e] <=> R[e]",pen=3)\cr
		additional_reactions_list[[3]]=list(id="TtoQ",name="TtoQ reaction",eq="T[e] => Q[e]",pen=5)\cr
		
		default: NULL\cr
		see also: \code{\link{additional_reactions_list}}
}
  \item{minimize}{
logical, specifies if blocked reaction should be removed from network before optimizing. May decrease solving time, but it takes time to calculate blocked reactions\cr
default: FALSE
}
  \item{simple}{
logical, if run in simple mode (TRUE) only the number of contradicting cases are minimized, network changes are not penalized\cr
default: FALSE 
}
  \item{verboseMode}{
numeric, should output be printed (1 => yes; !1 => no)\cr
default: 1
}
  \item{cancel_case_penalty}{
numerical>0, penalty for ignoring (case is violating viability threshold) a single growth case\cr
	default: NULL, penalty is higher than all network changes combined
}
  \item{not_delete_for}{
vector of reaction names; forward reactions that are not allowed to be removed (e.g., biomass objective function, exchange reactions)\cr
default: NULL\cr
see also: \code{\link{not_delete_for}}
}
  \item{not_delete_back}{
	vector of reaction names; backward reactions that are not allowed to be removed (e.g., biomass objective function, exchange reactions)\cr
	default: NULL\cr
	see also: \code{\link{not_delete_back}}
}
  \item{param_list}{
	list of specific paramaters for the solver\cr
	default: NULL\cr
	see also: \code{\link{p_list}}
}
  \item{use_indicator_constraints}{
logical,indicator constraints may prevent trickle flow, only usable if cplex (cplexAPI) is used as solver\cr
default:FALSE
}
  \item{stat_file}{
path of stat file\cr
	default: NULL
}
  \item{react_file}{
path of react file, contains all network modifications (subset of stat file)\cr
	default: NULL
}
  \item{reverse_reaction_list}{
list containing reaction, that are allowed to be reversed and according penalty. The following attributes must be defined for each entry:\cr
	reaction, character: name of reaction\cr
	pen, numeric>0, penalty for reversing reaction\cr
	see also: \code{\link{reverse_reaction_list}}

Example:\cr
reverse_reaction_list=list()\cr
reverse_reaction_list[[1]]=list(reaction="KtoT",pen=1)\cr
reverse_reaction_list[[2]]=list(reaction="TtoB",pen=1)

}
  \item{MaxPenalty}{
numerical >0; maximal penalty of all network changes\cr
	default NULL, sum of penalty of all network changes will be minimized 
}
  \item{alternatives}{
integer >=0, amount of alternative solution that should be calculated\cr
default: 0
}
  \item{bio_stoich}{
	stoichiometric coefficient for additional biomass metabolites\cr
	default: 1e-5
}
  \item{additional_biomass_metabolites}{
list of additional biomass metabolites\cr
		met specifies the metabolite, which can be added to the biomass objective function (note metabolites, which are already in the biomass objective function can not be added).\cr 
		pen specifies the corresponding penalty for adding this metabolite. \cr
		factor: -1 or 1; -1 metabolite can be added as substrate; 1 metabolite can be added as product\cr
		Example:\cr
		additional_biomass_mets=list()\cr
		additional_biomass_mets[[1]]=c(met="K[e]",pen=0.1,factor=1)\cr
		see also: \code{\link{additional_biomass_mets}}
}
  \item{remove_biomass_metabolites}{
list of metabolites that can be removed from the biomass objective function\cr
		met specifies the metabolite, which can be removed from the biomass objective function\cr
		pen specifies the corresponding penalty\cr
		Example:
		remove_biomass_mets=list()\cr
		remove_biomass_mets[[1]]=c(met="Z[e]",pen=0.01)\cr
		see also: \code{\link{remove_biomass_mets}}
}
}
\value{
optimized metabolic network model of type \code{modelorg}
}

\author{
	Daniel Hartleb
}
\examples{
\dontrun{
library("sybil")
library("GlobalFit")
library("cplexAPI")

##use cplex as solver
# Note, also alternative solvers can be used (i.e.,glpk[glpkAPI],gurobi[sybilGUROBI])
#, but GlobalFit has been tested only using cplex
SYBIL_SETTINGS("SOLVER", "cplexAPI")

#test=readTSVmod(reactList=network_path,quoteChar="\"")
data(test)

# names of reactions, which are not allowed to be removed
#not_delete_for=c(react_id(findExchReact(test)),"Biomass")
#not_delete_back=c(react_id(findExchReact(test)),"Biomass")
data(not_delete_for)
data(not_delete_back)

#names of reaction, which are not allowed to be reversed
#do_not_reverse=unique(c(react_id(findExchReact(test)),"Biomass"))

#set biomass object function
#obj_coef(test)[which(react_id(test)=="Biomass")]=1

#create influxes list
#influxes=list()
#influxes[[1]]=list(exRea="T[e]_import",value=-10)
data(influxes)


#apply influxes on network
#lowbnd(test)[pos=which(react_id(test)=="T[e]_import")]=-100

#growth cases
#on=c()
#on[[1]]=list(on=influxes,name="LIVE!",ko_react=c("TtoZ")
#,forced=TRUE,viability_threshold=1,gene_copy_number=1)
data(on)

#non-growth cases
#off=c()
#off[[1]]=list(on=influxes,name="DIE!",ko_react=c("KtoT"),
#forced=FALSE,viability_threshold=1,gene_copy_number=1)
data(off)

#optional parameter list for solver, in this example for cplex
#p_list=list(CPX_PARAM_THREADS=as.integer(1),CPX_PARAM_EPRHS=as.double(1e-9),
#CPX_PARAM_NETEPRHS=as.double(1e-11),CPX_PARAM_EPINT=as.double(1e-09),CPX_PARAM_TILIM=1e5,
#CPX_PARAM_PARALLELMODE=CPX_PARALLEL_OPPORTUNISTIC)
data(p_list)

#create list of reactions, that are allowed to be reversed
#reverse_reaction_list=list()
#reverse_reaction_list[[1]]=list(reaction="KtoT",pen=1)
#reverse_reaction_list[[2]]=list(reaction="TtoB",pen=1)
data(reverse_reaction_list)

#creact list of additional reactions
#additional_reactions_list=list()
#additional_reactions_list[[1]]=list(id="KtoB",name="KtoB reaction",
#eq="(2.1) K[e] => (1.2) B[e]",pen=7)
#additional_reactions_list[[2]]=list(id="TtoR",name="TtoR reaction",eq="T[e] <=> R[e] + Q[e]",pen=3)
#additional_reactions_list[[3]]=list(id="TtoQ",name="TtoQ reaction",eq="T[e] => Q[e]",pen=5)
data(additional_reactions_list)

#create list of additional biomass metabolites
#additional_biomass_mets=list()
#additional_biomass_mets[[1]]=c(met="Q[e]",pen=0.1,factor=-1)
#additional_biomass_mets[[2]]=c(met="B[e]",pen=0.1,factor=-1)
#additional_biomass_mets[[3]]=c(met="R[e]",pen=0.1,factor=-1)
data(additional_biomass_mets)

#create list of biomass metabolites, that are allowed to be removed
#remove_biomass_mets=list()
#remove_biomass_mets[[1]]=c(met="S[e]",pen=0.1)
#remove_biomass_mets[[2]]=c(met="T[e]",pen=0.1)
#remove_biomass_mets[[3]]=c(met="Z[e]",pen=0.1)
data(remove_biomass_mets)


opt_net=bilevel_optimize(network=test,on=on,off=off,additional_reactions=additional_reactions_list
,not_delete_for=not_delete_for,not_delete_back=not_delete_back,minimize=FALSE,simple=FALSE,
verboseMode=1,param_list=p_list,cancel_case_penalty=NULL,use_indicator_constraints=FALSE,
stat_file=NULL,react_file=NULL,reverse_reaction_list=reverse_reaction_list,alternatives=0,
MaxPenalty=NULL,additional_biomass_metabolites=additional_biomass_mets,
remove_biomass_metabolites=remove_biomass_mets)
}
}

