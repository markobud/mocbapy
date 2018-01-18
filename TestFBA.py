#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 14:38:16 2017
Test ModelSeed Models
@author: mbudinich
"""
#%% Loading packages
import cobra
#%% Function Definition
def get_objective(model):
    objs=dict()
    for r in model.reactions:
        if r.objective_coefficient != 0:
            objs[r]=r.objective_coefficient
    return(objs)

def pp_obj_val(obj,opt):
    return("\tObjective:{}\n\tValue:{}".format(obj,opt))

def find_orphans_metabolites(model):
    """Metabolites without associated reactions"""
    return [mets for mets in model.metabolites if len(mets.reactions) == 0]

def find_orphans_reactions(model):
    """Reactions without associated metabolites"""
    return [reaction for reaction in model.reactions if len(reaction.metabolites) == 0]
#%% Loading models
#Define Paths

syn_WH7803_path="Synechococcus_WH7803/32051.5.sbml"
syn_WH8102_path="Synechococcus_WH8102/84588.4.sbml"
syn_CC9605_path="Synechococcus_CC9605/110662.5.sbml" #Analogous to A15-62, Roscoff Collection
pro_MIT9312_path="Prochlorococus_MIT_9312/74546.5.sbml"
pro_MIT9313_path="Prochlorococus_MIT_9313/74547.5.sbml"
pro_CCMP1986_path="Prochlorococus_CCMP1986/59919.7.sbml" #EX MED4
can_HTCC1062_path="Candidatus_Pelagibacter_ubique_HTCC1062/335992.9.sbml" #SAR11
paths = [syn_WH7803_path,syn_WH8102_path,syn_CC9605_path,pro_MIT9312_path,pro_MIT9313_path,pro_CCMP1986_path,can_HTCC1062_path]

#Load Models
models = [cobra.io.read_sbml_model(x) for x in paths]

for model, path in list(zip(models,paths)):
    model.id = path + ":" + model.id

#%%Testing FBA
for model in models:
    obj=get_objective(model)
    opt=model.optimize(solver='gurobi')
    print("{}\n{}".format(model.id,pp_obj_val(obj,opt)))
    print("\tReactions:{}, Metabolites:{}".format(len(model.reactions),len(model.metabolites)))
    print("\tOrphan Reactions:{}, Orphan Metabolites:{}".format(len(find_orphans_reactions(model)),len(find_orphans_metabolites(model))))
    print("***Prunning Model***")
    blocked_reactions = cobra.flux_analysis.find_blocked_reactions(model,solver='gurobi',zero_cutoff=1e-4,open_exchanges=True)
    model.remove_reactions(blocked_reactions,delete=True,remove_orphans=True)
    opt=model.optimize(solver='gurobi')
    print("{}\n{}".format(model.id,pp_obj_val(obj,opt)))
    print("\tReactions:{}, Metabolites:{}".format(len(model.reactions),len(model.metabolites)))
    print("\tOrphan Reactions:{}, Orphan Metabolites:{}".format(len(find_orphans_reactions(model)),len(find_orphans_metabolites(model))))
    print("\n\n")