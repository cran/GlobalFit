\name{GlobalFit-package}
\alias{GlobalFit-package}
\alias{GlobalFit}
\docType{package}
\title{
\packageTitle{GlobalFit}
}
\description{
Initial metabolic networks often inaccurately predict in-silico growth or non-growth if compared to in-vivo data. This package refines metabolic network models by making networks changes (i.e., removing, adding, changing reversibility of reactions; adding and removing biomass metabolites) and simultaneously matching sets of experimental growth and non-growth data (e.g., KO-mutants, mutants grown under different media conditions,...). Three versions of GlobalFit are provided; an old(algorithm=2), a newer (faster) version (algorithm=1) and a third version (algortithm=3) which can be used to remove thermodynamically infeasible cyles. 
}
\details{
\tabular{ll}{
Package: \tab GlobalFit\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2016-08-12\cr
License: \tab GPL-3\cr
}}
\author{
Daniel Hartleb

daniel.hartleb@hhu.de
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
influxes[[1]]=list(exRea="A[e]_import",value=-10)

#set influxes
lowbnd(example_net1)[pos=which(react_id(example_net1)=="A[e]_import")]=-10

#growth cases
on=list()
on[[1]]=list(on=influxes,name="LIVE!",ko_react=c("AtoB"),forced=TRUE,viability_threshold=1,
 gene_copy_number=1)
on[[2]]=list(on=influxes,name="LIVE!",ko_react=c("AtoB"),forced=TRUE,viability_threshold=1,
 gene_copy_number=1)
#non-growth cases
off=list()
off[[1]]=list(on=influxes,name="DIE!",ko_react=c("AtoR"),forced=FALSE,viability_threshold=1,
 gene_copy_number=1)
off[[2]]=list(on=influxes,name="DIE!",ko_react=c("AtoR"),forced=FALSE,viability_threshold=1,
 gene_copy_number=1)
#optional parameter list for solver, in this example for cplex
p_list=list(CPX_PARAM_THREADS=as.integer(1),CPX_PARAM_EPRHS=as.double(1e-9),
CPX_PARAM_NETEPRHS=as.double(1e-11),CPX_PARAM_EPINT=as.double(1e-09),
CPX_PARAM_TILIM=1e5,CPX_PARAM_PARALLELMODE=CPX_PARALLEL_OPPORTUNISTIC)


#create list of reactions, that are allowed to be reversed
reverse_reaction_list=list()
reverse_reaction_list[[1]]=list(reaction="AtoC",pen=1)
reverse_reaction_list[[2]]=list(reaction="TtoB",pen=1)

#create list of additional reactions
additional_reactions_list=list()
additional_reactions_list[[1]]=list(id="KtoB",name="KtoB reaction",eq="(2.1) K[e] => B[e]",pen=7)
additional_reactions_list[[2]]=list(id="TtoR",name="TtoR reaction",eq="T[e] <=> R[e] + Q[e]",pen=3)
additional_reactions_list[[3]]=list(id="CtoB",name="TtoQ reaction",eq="C[e] => B[e]",pen=5)

#create list of additional biomass metabolites
additional_biomass_mets=list()
additional_biomass_mets[[1]]=list(met="Q[e]",pen=0.1,factor=-1)
additional_biomass_mets[[2]]=list(met="R[e]",pen=0.1,factor=-1)
additional_biomass_mets[[3]]=list(met="B[e]",pen=0.1,factor=-1)

#create list of biomass metabolites, that are allowed to be removed
remove_biomass_mets=list()
remove_biomass_mets[[1]]=list(met="S[e]",pen=40.1)
remove_biomass_mets[[2]]=list(met="T[e]",pen=0.1)
remove_biomass_mets[[3]]=list(met="Z[e]",pen=0.1)

#set penalties for removing reactions (network contains nine reactions, so we have to set 9 values)
remove_penalties_hin=c(1,2.5,3,4,5,6,7,8,9)
remove_penalties_back=c(1,2.5,3,4,5,6,7,8,9)

opt_net=bilevel_optimize(network=example_net1,on=on,off=off,algorithm=1,
additional_reactions=additional_reactions_list,not_delete_for=not_delete_for,
not_delete_back=not_delete_back,minimize=FALSE,simple=FALSE,verboseMode=1,
param_list=p_list,cancel_case_penalty=NULL,use_indicator_constraints=FALSE,
stat_file=NULL,react_file=NULL,remove_penalties_hin=remove_penalties_hin,
remove_penalties_back=remove_penalties_back,reverse_reaction_list=reverse_reaction_list,
alternatives=0,MaxPenalty=NULL,additional_biomass_metabolites=additional_biomass_mets,
remove_biomass_metabolites=remove_biomass_mets,variable_lower_bound=NULL,forced_modifications=0)



##############
##EXAMPLE2: NETWORK CONTAINS THERMODYNAMIC INFEASIBLE CYCLES
# (CYC1 AND CYC2 CAN CARRY FLUX WITHOUT ANY INFLUX); 
# WE USE GLOBALFIT AND DIFFERENT BIOMASS OBEJCTIVE FUNCTIONS 
# FOR THE GROWTH AND NON-GROWTH CASE

data(example_net2)


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

#non-growth cases with different biomass ("CYC1","CYC2")
off=list()
off[[1]]=list(on=influxes,name="DIE!",ko_react=c(),forced=FALSE,viability_threshold=1,
 gene_copy_number=1,biomass=c("CYC1","CYC2"))
 
# no alternative modifications allowed
reverse_reaction_list=list()
additional_reactions_list=list()
additional_biomass_mets=list()
remove_biomass_mets=list()
 
# names of reactions, which are not allowed to be removed, including the cycle reactions
not_delete_for=c(react_id(findExchReact(example_net2)),"Biomass","CYC1","CYC2")
not_delete_back=c(react_id(findExchReact(example_net2)),"Biomass","CYC1","CYC2")




opt_net=bilevel_optimize(network=example_net2,on=on,off=off,algorithm=3,
additional_reactions=additional_reactions_list,not_delete_for=not_delete_for,
not_delete_back=not_delete_back,minimize=FALSE,simple=FALSE,verboseMode=1,
param_list=p_list,cancel_case_penalty=NULL,use_indicator_constraints=FALSE,
stat_file=NULL,react_file=NULL,reverse_reaction_list=reverse_reaction_list,
alternatives=0,MaxPenalty=NULL,additional_biomass_metabolites=additional_biomass_mets,
remove_biomass_metabolites=remove_biomass_mets,
variable_lower_bound=NULL,forced_modifications=0)





##########################
##EXAMPLE3: NON-GROWTH CASE CAN ONLY BE RESOLVED BY CHANGING THE LOWER BOUND OF AN INFLUX 
#(ONLY WORKS WITH THE SLOWER IMPLEMENTATION OF GLOBALFIT; ALGORITHM=2). 
#THIS CAN BE USED TO FIND SUITABLE QUALITATIVE MEDIA COMPOSITIONS. 
#NOTE IN THIS SIMPLE EXAMPLE THE VIABILITY THRESHOLD 
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
# Because the viability threshold of the non-growth case is 2 
# and the viability threshold of the growth case is 1; 
# the optimized value should be between -2 and -1 

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
remove_biomass_metabolites=remove_biomass_mets,variable_lower_bound=varying_lower_bound_list,
forced_modifications=0)
}
}
