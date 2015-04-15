\name{GlobalFit-package}
\alias{GlobalFit-package}
\alias{GlobalFit}
\docType{package}
\title{
GlobalFit finds minimal set of network changes for a metabolic network to improve accuracy
}
\description{
Initial metabolic networks often inaccurately predict in-silico growth or non-growth if compared to in-vivo data. This package refines metabolic network models by making networks changes (i.e., removing, adding, changing reversibility of reactions; adding and removing biomass metabolites) and simultaneously matching sets of experimental growth and non-growth data (e.g., KO-mutants, mutants grown under different media conditions,...)
}
\details{
\tabular{ll}{
Package: \tab GlobalFit\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-04-20\cr
License: \tab GPL-3\cr
}
}
\author{
Daniel Hartleb

daniel.hartleb@hhu.de
}


\examples{
\dontrun{
library("sybil")
library("GlobalFit")
library("cplexAPI")

##use cplex as solver
# Note, also alternative solvers can be used (i.e.,glpk[glpkAPI],gurobi[sybilGUROBI]), 
# but GlobalFit has been tested only using cplex
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
#on[[1]]=list(on=influxes,name="LIVE!",ko_react=c("TtoZ"),forced=TRUE,viability_threshold=1,
# gene_copy_number=1)
data(on)

#non-growth cases
#off=c()
#off[[1]]=list(on=influxes,name="DIE!",ko_react=c("KtoT"),forced=FALSE,viability_threshold=1,
# gene_copy_number=1)
data(off)

#optional parameter list for solver, in this example for cplex
#p_list=list(CPX_PARAM_THREADS=as.integer(1),CPX_PARAM_EPRHS=as.double(1e-9),
#CPX_PARAM_NETEPRHS=as.double(1e-11),CPX_PARAM_EPINT=as.double(1e-09),
#CPX_PARAM_TILIM=1e5,CPX_PARAM_PARALLELMODE=CPX_PARALLEL_OPPORTUNISTIC)
data(p_list)

#create list of reactions, that are allowed to be reversed
#reverse_reaction_list=list()
#reverse_reaction_list[[1]]=list(reaction="KtoT",pen=1)
#reverse_reaction_list[[2]]=list(reaction="TtoB",pen=1)
data(reverse_reaction_list)

#creact list of additional reactions
#additional_reactions_list=list()
#additional_reactions_list[[1]]=list(id="KtoB",
#name="KtoB reaction",eq="(2.1) K[e] => (1.2) B[e]",pen=7)
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


opt_net=bilevel_optimize(network=test,on=on,off=off,additional_reactions=additional_reactions_list,
not_delete_for=not_delete_for,not_delete_back=not_delete_back,minimize=FALSE,simple=FALSE,
verboseMode=1,param_list=p_list,cancel_case_penalty=NULL,use_indicator_constraints=FALSE,
stat_file=NULL,react_file=NULL,reverse_reaction_list=reverse_reaction_list,alternatives=0,
MaxPenalty=NULL,additional_biomass_metabolites=additional_biomass_mets,
remove_biomass_metabolites=remove_biomass_mets)
}
}
