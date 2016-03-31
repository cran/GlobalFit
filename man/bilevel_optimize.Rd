\name{bilevel_optimize}
\alias{bilevel_optimize}
\title{
pre und post-processing of optimization
}
\description{
pre und post-processing of optimization. Adds reaction to network, creates output, creates optimization object, interprets solution, applies modification to network.
}
\usage{
bilevel_optimize(network, on = c(), off = c(), algorithm = 1, additional_reactions = NULL,
 minimize = TRUE, simple = FALSE, verboseMode = 1, cancel_case_penalty = NULL, 
 not_delete_for = c(), not_delete_back = c(), param_list = NULL, 
 use_indicator_constraints = FALSE, stat_file = NULL, react_file = NULL, 
 reverse_reaction_list = NULL, MaxPenalty = NULL, alternatives = 0, bio_stoich = 1e-05, 
 additional_biomass_metabolites = NULL, remove_biomass_metabolites = NULL, 
 variable_lower_bound = NULL, forced_modifications = 0)
}

\arguments{
 \item{network}{
metabolic network model of type \code{modelorg}
}
  \item{on}{
each entry must contain:
		on: list of influxes\cr
		name: character, name of growth case\cr
		ko_react: vector of reaction, which are knocked out (i.e.,lower und upper bound = 0)\cr
		forced: logical, FALSE (growth case can be ignored with according penalty) or TRUE (case cannot be ignored)\cr
		viability_threshold: numerical>0 threshold, which is considered as growth\cr
		gene_copy_number: integer >0, multiplies the penalty for ignoring this growth case\cr
		biomass: specifies the biomass objective function of this growth case\cr
		\cr
		Example:\cr
		on=list()\cr
		on[[1]]=list(on=influxes,name="LIVE!",ko_react=c(),forced=TRUE,\cr
		viability_threshold=1,gene_copy_number=1,biomass="Biomass")
 
		default: NULL\cr
}
  \item{off}{
list of non-growth cases
	each entry must contain:
		on: list of influxes\cr
		name: character, name of non-growth case\cr
		ko_react: vector of reaction, which are knocked out (i.e.,lower und upper bound = 0)\cr
		forced: logical, FALSE (growth case can be ignored with according penalty) or TRUE (case cannot be ignored)\cr
		viability_threshold: numerical>0 threshold, which is considered as growth\cr
		gene_copy_number=:integer >0, multiplies the penalty for ignoring this non-growth case\cr
		biomass: specifies the biomass objective function of this growth case
		
		Example:\cr
		off=list()\cr
off[[1]]=list(on=influxes,name="DIE!",ko_react=c(),forced=FALSE,\cr
viability_threshold=1,gene_copy_number=1,biomass="A[e]_import")
		
 
		default: NULL\cr
}
  \item{algorithm}{
 integer, specifies which version of GlobalFit should be used:\cr
  1: fast version\cr
  2: old version, but allows to use variable_lower_bound\cr
  
  default: 1\cr
  }
  \item{additional_reactions}{
list containing additional reaction. Each entry must contain the following attributes.\cr
		id: character, id of reaction\cr
		name: character, name of reaction \cr
		eq: character, equation of reaction.\cr
		pen: numeric >0, penalty for adding this reaction to the network\cr
		
		Example:\cr
		additional_reactions_list=list()\cr
		additional_reactions_list[[1]]=list(id="KtoB",name="KtoB reaction",eq="(2.1) K[e] => B[e]",pen=7)\cr
		additional_reactions_list[[2]]=list(id="TtoR",name="TtoR reaction",eq="T[e] <=> R[e]",pen=3)\cr
		additional_reactions_list[[3]]=list(id="TtoQ",name="TtoQ reaction",eq="T[e] => Q[e]",pen=5)\cr
		
		default: NULL\cr
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

}
  \item{not_delete_back}{
	vector of reaction names; backward reactions that are not allowed to be removed (e.g., biomass objective function, exchange reactions)\cr
	default: NULL\cr
	
}
  \item{param_list}{
	list of specific paramaters for the solver\cr
	default: NULL\cr
	
}
  \item{use_indicator_constraints}{
logical,indicator constraints may prevent trickle flow, only usable if cplex (cplexAPI) is used as solver\cr
default:FALSE
}
  \item{stat_file}{
path for stat file\cr
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
	

Example:\cr
reverse_reaction_list=list()\cr
reverse_reaction_list[[1]]=list(reaction="KtoT",pen=1)\cr
reverse_reaction_list[[2]]=list(reaction="TtoB",pen=1)\cr

default: NULL\cr
}
  \item{MaxPenalty}{
integer >=0, amount of alternative solution that should be calculated\cr
default: 0
}
  \item{alternatives}{
integer >=0, amount of alternative solution that should be calculated\cr
default: 0
}
  \item{bio_stoich}{
numeric>0, stoichiometric coefficient for additional biomass metabolites\cr
	default: 1e-5
}
  \item{additional_biomass_metabolites}{
list of additional biomass metabolites\cr
		met specifies the metabolite, which can be added to the biomass objective function (note metabolites, which are already in the biomass objective function can not be added).\cr 
		pen specifies the corresponding penalty for adding this metabolite. \cr
		factor: -1 or 1; -1 metabolite can be added as substrate; 1 metabolite can be added as product\cr
		Example:\cr
		additional_biomass_mets=list()\cr
		additional_biomass_mets[[1]]=list(met="K[e]",pen=0.1,factor=1)\cr
		
		default: NULL\cr
}
  \item{remove_biomass_metabolites}{
list of metabolites that can be removed from the biomass objective function\cr
		met specifies the metabolite, which can be removed from the biomass objective function\cr
		pen specifies the corresponding penalty\cr
		
		Example:\cr
		remove_biomass_mets=list()\cr
		remove_biomass_mets[[1]]=list(met="Z[e]",pen=0.01)\cr
		
		default: NULL\cr	
}
  \item{variable_lower_bound}{
list containing reactions, which lowerbound can be optimized.\cr 
	This can be used to calculate an qualitatively optimized media formulation\cr
	Can only be used if algorithm=2\cr
	
	reaction specifies the name of the reaction\cr
	min specifies the minimal possible value\cr
	max specifies the maximal possible value\cr
	penalty specifies the penalty for changing the lowerbound\cr
	
	Example:\cr
	varying_lower_bound_list=list()\cr
	varying_lower_bound_list[[1]]=list(reaction="T[e]_import",min=-20,max=-0,penalty=0.1)\cr
	
	default: NULL\cr
}
  \item{forced_modifications}{
integer>=0,number of minimal modification; may reduce computation time if set >=1\cr
	default: 0
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




library(sybil)
library(GlobalFit)
library("cplexAPI")
SYBIL_SETTINGS("SOLVER", "cplexAPI")
#SYBIL_SETTINGS("SOLVER", "sybilGUROBI")



###################
##EXAMPLE1: RECONCILATION OF TWO FALSE PREDICTIONS


data(example_net1)


# names of reactions, which are not allowed to be removed
not_delete_for=c(react_id(findExchReact(example_net1)),"Biomass")
not_delete_back=c(react_id(findExchReact(example_net1)),"Biomass")

#set biomass object function
obj_coef(example_net1)[which(react_id(example_net1)=="Biomass")]=1

#create list of influxes
influxes=list()
influxes[[1]]=list(exRea="T[e]_import",value=-10)

#set influxes
lowbnd(example_net1)[pos=which(react_id(example_net1)=="T[e]_import")]=-10

#growth cases
on=list()
on[[1]]=list(on=influxes,name="LIVE!",ko_react=c("TtoZ"),forced=TRUE,viability_threshold=1,
 gene_copy_number=1)

#non-growth cases
off=list()
off[[1]]=list(on=influxes,name="DIE!",ko_react=c("KtoT"),forced=FALSE,viability_threshold=1,
 gene_copy_number=1)

#optional parameter list for solver, in this example for cplex
p_list=list(CPX_PARAM_THREADS=as.integer(1),CPX_PARAM_EPRHS=as.double(1e-9),
CPX_PARAM_NETEPRHS=as.double(1e-11),CPX_PARAM_EPINT=as.double(1e-09),
CPX_PARAM_TILIM=1e5,CPX_PARAM_PARALLELMODE=CPX_PARALLEL_OPPORTUNISTIC)


#create list of reactions, that are allowed to be reversed
reverse_reaction_list=list()
reverse_reaction_list[[1]]=list(reaction="KtoT",pen=1)
reverse_reaction_list[[2]]=list(reaction="TtoB",pen=1)

#create list of additional reactions
additional_reactions_list=list()
additional_reactions_list[[1]]=list(id="KtoB",name="KtoB reaction",eq="(2.1) K[e] => B[e]",pen=7)
additional_reactions_list[[2]]=list(id="TtoR",name="TtoR reaction",eq="T[e] <=> R[e] + Q[e]",pen=3)
additional_reactions_list[[3]]=list(id="TtoQ",name="TtoQ reaction",eq="T[e] => Q[e]",pen=5)

#create list of additional biomass metabolites
additional_biomass_mets=list()
additional_biomass_mets[[1]]=list(met="Q[e]",pen=0.1,factor=-1)
additional_biomass_mets[[2]]=list(met="B[e]",pen=0.1,factor=-1)
additional_biomass_mets[[3]]=list(met="R[e]",pen=0.1,factor=-1)

#create list of biomass metabolites, that are allowed to be removed
remove_biomass_mets=list()
remove_biomass_mets[[1]]=list(met="S[e]",pen=0.1)
remove_biomass_mets[[2]]=list(met="T[e]",pen=0.1)
remove_biomass_mets[[3]]=list(met="Z[e]",pen=0.1)



opt_net=bilevel_optimize(network=example_net1,on=on,off=off,algorithm=1,
additional_reactions=additional_reactions_list,not_delete_for=not_delete_for,
not_delete_back=not_delete_back,minimize=FALSE,simple=FALSE,verboseMode=1,
param_list=p_list,cancel_case_penalty=NULL,use_indicator_constraints=FALSE,
stat_file=NULL,react_file=NULL,reverse_reaction_list=reverse_reaction_list,
alternatives=0,MaxPenalty=NULL,additional_biomass_metabolites=additional_biomass_mets,
remove_biomass_metabolites=remove_biomass_mets,variable_lower_bound=NULL,forced_modifications=0)



##############
##EXAMPLE2: NETWORK CONTAINS THERMODYNAMIC CYCLE 
#(IT CAN PRODUCE UNLIMITED A[e] WITHOUT ANY INFLUX); 
#WE USE GLOBALFIT AND TWO DIFFERENT BIOMASS OBEJCTIVE FUNCTIONS 
#FOR THE GROWTH AND NON-GROWTH CASE

data(example_net2)

# names of reactions, which are not allowed to be removed
not_delete_for=c(react_id(findExchReact(example_net2)),"Biomass")
not_delete_back=c(react_id(findExchReact(example_net2)),"Biomass")

#set wild type biomass object function
obj_coef(example_net2)[which(react_id(example_net2)=="Biomass")]=1

#create 2 lists of influxes (one list is empty)
influxes=list()
influxes[[1]]=list(exRea="T[e]_import",value=-10)

influxes2=list()

#set influxes for wild type
lowbnd(example_net2)[pos=which(react_id(example_net2)=="T[e]_import")]=-10

#growth cases with wild type biomass
on=list()
on[[1]]=list(on=influxes,name="LIVE!",ko_react=c(),forced=TRUE,viability_threshold=1,
 gene_copy_number=1,biomass="Biomass")

#non-growth cases with different biomass(A[e]_imoprt)
off=list()
off[[1]]=list(on=influxes,name="DIE!",ko_react=c(),forced=FALSE,viability_threshold=1,
 gene_copy_number=1,biomass="A[e]_import")
 
# no alternative modifications allowed
reverse_reaction_list=list()
additional_reactions_list=list()
additional_biomass_mets=list()
remove_biomass_mets=list()
 
opt_net=bilevel_optimize(network=example_net2,on=on,off=off,algorithm=1,
additional_reactions=additional_reactions_list,not_delete_for=not_delete_for,
not_delete_back=not_delete_back,minimize=FALSE,simple=FALSE,verboseMode=1,
param_list=p_list,cancel_case_penalty=NULL,use_indicator_constraints=FALSE,
stat_file=NULL,react_file=NULL,reverse_reaction_list=reverse_reaction_list,
alternatives=0,MaxPenalty=NULL,additional_biomass_metabolites=additional_biomass_mets,
remove_biomass_metabolites=remove_biomass_mets,variable_lower_bound=NULL,forced_modifications=0)

##########################
##EXAMPLE3: NON-GROWTH CASE CAN ONLY BE RESOLVED BY CHANGING THE LOWER BOUND OF AN INFLUX 
#(ONLY WORKS WITH THE SLOWER IMPLEMENTATION OF GLOBALFIT; ALGORITHM=2). 
#THIS CAN BE USED TO FIND SUITABLE QUALITATIVE MEDIA COMPOSTION. 
#NOTE IN THIS SIMPLE EXAMPLE THE VIABITLITY THRESHOLD 
#OF THE NON-GROWTH CASE IS HIGHER THAN THE GROWTH CASE

data(example_net3)

# names of reactions, which are not allowed to be removed
not_delete_for=c(react_id(findExchReact(example_net3)),"Biomass")
not_delete_back=c(react_id(findExchReact(example_net3)),"Biomass")

#set wild type biomass object function
obj_coef(example_net3)[which(react_id(example_net3)=="Biomass")]=1

#create 2 lists of influxes (one list is empty)
influxes=list()
influxes[[1]]=list(exRea="T[e]_import",value=-100)

influxes2=list()

#set influxes for wild type
lowbnd(example_net3)[pos=which(react_id(example_net3)=="T[e]_import")]=-100

#growth cases with wild type biomass
on=list()
on[[1]]=list(on=influxes,name="LIVE!",ko_react=c(),forced=TRUE,viability_threshold=1,
 gene_copy_number=1)

#non-growth cases with different biomass(A[e]_imoprt)
off=list()
off[[1]]=list(on=influxes,name="DIE!",ko_react=c(),forced=FALSE,viability_threshold=2,
 gene_copy_number=1)

## set varying_lower_bound; T[e]_import is allowed to vary between 0 and -20. 
#Because the viability threshold of the non-growth case is 2 and 
#the viability threshold of the growth case is 1; 
#the optimized value should be between -2 and -1 


varying_lower_bound_list=list()
varying_lower_bound_list[[1]]=list(reaction="T[e]_import",min=-20,max=-0,penalty=0.1)

 
# no alternative modifications allowed
reverse_reaction_list=list()
additional_reactions_list=list()
additional_biomass_mets=list()
remove_biomass_mets=list()
 
opt_net=bilevel_optimize(network=example_net3,on=on,off=off,algorithm=2,
additional_reactions=additional_reactions_list,not_delete_for=not_delete_for,
not_delete_back=not_delete_back,minimize=FALSE,simple=FALSE,verboseMode=1,
param_list=p_list,cancel_case_penalty=NULL,use_indicator_constraints=FALSE,
stat_file=NULL,react_file=NULL,reverse_reaction_list=reverse_reaction_list,
alternatives=0,MaxPenalty=NULL,additional_biomass_metabolites=additional_biomass_mets,
remove_biomass_metabolites=remove_biomass_mets,
variable_lower_bound=varying_lower_bound_list,forced_modifications=0)

}
}

