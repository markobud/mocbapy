#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:35:50 2017
Construct Ecosystem system
@author: mbudinich
"""
#import cobra
import numpy as np
from scipy.sparse import lil_matrix, find, block_diag, eye
from benpy import vlpProblem, vlpSolution, solve as bensolve
from collections import OrderedDict
#%%
class EcosystemModel:

    from scipy.sparse import lil_matrix, block_diag, eye 
    
    def __init__(self,model_array=None,metabolic_dict=None):
        """Instantiate the EcosystemModel object
        model_array is an array of cobra models to connect
        metabolic_dict is a dictionary where its keys correspond to metabolites id's and their values to the name equivalence"""
        self.models=model_array
        self.metabolic_dict = metabolic_dict
        self._pooldict = None
        self.pool_ex_rxns = None
        self.pool_ex_mets = None
        self.Ssigma = None
        self.sysreactions = None
        self.sysmetabolites = None
        self.lb = None
        self.ub = None
        self.objectives = None

    def construct_ecosystem_pool(self):
        """Check all metabolites used in import/export exchanges and construct the pool compartment"""
        pooldict=dict()
        for model in self.models:
            for rxn_ex in model.exchanges:
                for met_ex in rxn_ex.metabolites:
                    if met_ex.id in self.metabolic_dict:
                        met_name = self.metabolic_dict[met_ex.id]
                        if  met_name not in pooldict:
                            pooldict[met_name]=[(model,rxn_ex,rxn_ex.get_coefficient(met_ex.id))]
                        else:
                            pooldict[met_name].append((model,rxn_ex,rxn_ex.get_coefficient(met_ex.id)))
        self._pooldict = pooldict        

    def populate_ecosystem_model(self):
        """Calculate the object attributes after pool construction"""
        self.pool_ex_rxns = []
        self.pool_ex_mets = []
        pool_ub = []
        pool_lb = []
        for key in self._pooldict.keys():
            self.pool_ex_rxns.append("EX_{}:pool".format(key))
            self.pool_ex_mets.append("{}:pool".format(key))
            pool_lb.append(-1000)
            pool_ub.append(1000)
        self.sysreactions=[]
        self.sysmetabolites=[]
        self.lb = []
        self.ub = []
        self.objectives = OrderedDict()
        for model in self.models:
            self.objectives[model.id] = []
        for model in self.models:
            for r in model.reactions:
                new_name = "{}:{}".format(r.id,model.id)
                self.sysreactions.append(new_name)
                self.lb.append(r.lower_bound)
                self.ub.append(r.upper_bound)
                if r.objective_coefficient != 0:
                    self.objectives[model.id].append((r.id,r.objective_coefficient))
            for m in model.metabolites:
                self.sysmetabolites.append("{}:{}".format(m.id,model.id))
        self.sysreactions.extend(self.pool_ex_rxns)
        self.sysmetabolites.extend(self.pool_ex_mets)
        self.lb.extend(pool_lb)
        self.ub.extend(pool_ub)
        array_form = [model.to_array_based_model() for model in self.models]
        self.Ssigma = block_diag([model.S for model in array_form])
        self.Ssigma = lil_matrix(block_diag([self.Ssigma,-eye(len(self.pool_ex_rxns))]))
        for met in self._pooldict.keys():
            met_name = "{}:pool".format(met)
            met_idx = self.sysmetabolites.index(met_name)
            for model,reaction,coeff in self._pooldict[met]:
                rxn_name = "{}:{}".format(reaction.id,model.id)
                rxn_idx = self.sysreactions.index(rxn_name)
                self.Ssigma[met_idx,rxn_idx]=-coeff
    
    def add_comparment(self,model):
        """Utility function to add a new agent to models. Pretty ineficient, re-runs all the steps again for each addition"""
        self.__init__(self.model_array.add(model),self.metabolic_dict)
        self.construct_ecosystem_pool()
        self.populate_ecosystem_model()
   
    @staticmethod
    def get_common_mets(model_list):
        """Naive implementation of getting common exchange metabolites, using their id's"""
        common_mets=dict()
        for model in model_list:
            for rxn_ex in model.exchanges:
                for met_ex in rxn_ex.metabolites:
                    if met_ex.id not in common_mets:
                        common_mets[met_ex.id]=dict([(model,rxn_ex)])
                    else:
                        common_mets[met_ex.id][model]=rxn_ex
        return(common_mets)
    def to_vlp(self,filename=None):
        """Returns a vlp problem from EcosystemModel"""
        # We are using bensolve-2.0.1:
        # B is coefficient matrix
        # P is objective Marix
        # a is lower bounds for B
        # b is upper bounds for B
        # l is lower bounds of variables
        # s is upper bounds of variables
        # opt_dir is direction: 1 min, -1 max
        # Y,Z and c are part of cone definition. If empty => MOLP
        vlp = vlpProblem()
        if filename != None:
            vlp.filename = filename
        m, n = self.Ssigma.shape
        q = len(self.models)
        vlp.B = self.Ssigma
        vlp.a = np.zeros((1, m))[0]
        vlp.b = np.zeros((1, m))[0]
        vlp.l = self.lb
        vlp.s = self.ub
        vlp.P = lil_matrix((q, n))
        vlp.opt_dir = -1
        for i in range(q):
            for rxn, coeff in self.objectives[self.models[i].id]:
                new_name = "{}:{}".format(rxn,self.models[i].id)
                k = self.sysreactions.index(new_name)
                print((i,k))
                vlp.P[i,k] = coeff
        vlp.Y = None
        vlp.Z = None
        vlp.c = None
        return(vlp)
#%%
