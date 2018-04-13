#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:35:50 2017
Construct Ecosystem model
@author: mbudinich
"""
from cobra.util.array import create_stoichiometric_matrix
from numpy import zeros
from scipy.sparse import lil_matrix, block_diag, eye
from benpy import vlpProblem
from collections import OrderedDict, defaultdict
from warnings import warn
from tqdm import tqdm

class EcosystemModel:

    def _construct_ecosystem_pool(self):
        """Check all metabolites used in import/export exchanges and construct the pool compartment"""
        pooldict = defaultdict(list)
        for model in self.models:
            for rxn_ex in model.exchanges:
                for met_ex in rxn_ex.metabolites:
                    if (met_ex.id, model) in self.metabolic_dict:
                        met_name = self.metabolic_dict[(met_ex.id, model)]
                        pooldict[met_name].append((model, rxn_ex, rxn_ex.get_coefficient(met_ex.id)))
        self._pooldict = dict(pooldict)

    def _populate_ecosystem_model(self):
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
        self.sysreactions = []
        self.sysmetabolites = []
        self.lb = []
        self.ub = []
        self.objectives = OrderedDict()
        for model in self.models:
            self.objectives[model] = []
        for model in self.models:
            for r in model.reactions:
                new_name = "{}:{}".format(r.id, model.id)
                self.sysreactions.append(new_name)
                self.lb.append(r.lower_bound)
                self.ub.append(r.upper_bound)
                if r.objective_coefficient != 0:
                    self.objectives[model].append((r.id, r.objective_coefficient))
            for m in model.metabolites:
                self.sysmetabolites.append("{}:{}".format(m.id, model.id))
        self.sysreactions.extend(self.pool_ex_rxns)
        self.sysmetabolites.extend(self.pool_ex_mets)
        self.lb.extend(pool_lb)
        self.ub.extend(pool_ub)

        array_form = block_diag([create_stoichiometric_matrix(model, array_type="lil")
                                 for model in self.models],
                                format="lil")

        self.Ssigma = block_diag([array_form, -eye(len(self.pool_ex_rxns))], format="lil")
        for met in self._pooldict.keys():
            met_name = "{}:pool".format(met)
            met_idx = self.sysmetabolites.index(met_name)
            for model, reaction, coeff in self._pooldict[met]:
                rxn_name = "{}:{}".format(reaction.id, model.id)
                rxn_idx = self.sysreactions.index(rxn_name)
                self.Ssigma[met_idx, rxn_idx] = -coeff

    def build_base_opt_model(self, solver=None):
        """ Builds the underlying base optimization problem Sv = 0, lb <= v <= ub """
        interfase = _choose_optlang_interfase(solver)
        model = interfase.Model(name='Base Solver Model')
        m, n = self.Ssigma.shape
        assert m == len(self.sysmetabolites)
        assert n == len(self.sysreactions)
        # Create flux variables
        flux_variables = [interfase.Variable(rxn, lb=ecosystem_model.lb[i], ub=ecosystem_model.ub[i]) for i, rxn in
                          enumerate(ecosystem_model.sysreactions)]
        model.add(flux_variables, sloppy=True)
        model.update()
#        for i in tqdm(range(m)):
#            terms_const = [flux_variables[j] * ecosystem_model.Ssigma[i, j] for j in range(n) if ecosystem_model.Ssigma[i, j] != 0]
#            mass_const = interfase.Constraint(sum(terms_const), lb=0, ub=0)
#            model.add(mass_const, sloppy=True)
#        model.update()
#        model.objective = interfase.Objective(flux_variables[0], direction="max")
        return model

    def __init__(self, model_array=None, metabolic_dict=None):

        """Instantiate the EcosystemModel object model_array is an array of cobra models to connect
        metabolic_dict is a dictionary such as:
            * Its keys correspond to tuples (metabolite_id,model)
            * Its value correspond to the id that will be used in the model
        """

        self.models = model_array
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
        if model_array is not None and metabolic_dict is not None:
            self._construct_ecosystem_pool()
            self._populate_ecosystem_model()
        elif model_array is None:
            warn("Models array is empty")
        elif metabolic_dict is None:
            warn("No metabolic dictionary is given")

    def add_comparment(self, model):
        """Utility function to add a new agent to models.
        Pretty inefficient, re-runs all the steps again for each addition"""
        self.__init__(self.models.add(model), self.metabolic_dict)

    def to_vlp(self, **kwargs):
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
        vlp = vlpProblem(**kwargs)
        m, n = self.Ssigma.shape
        q = len(self.models)
        vlp.B = self.Ssigma
        vlp.a = zeros((1, m))[0]
        vlp.b = zeros((1, m))[0]
        vlp.l = self.lb
        vlp.s = self.ub
        vlp.P = lil_matrix((q, n))
        vlp.opt_dir = -1
        for i in range(q):
            for rxn, coeff in self.objectives[self.models[i]]:
                new_name = "{}:{}".format(rxn, self.models[i].id)
                k = self.sysreactions.index(new_name)
                print((i, k))
                vlp.P[i, k] = coeff
        vlp.Y = None
        vlp.Z = None
        vlp.c = None
        return vlp
