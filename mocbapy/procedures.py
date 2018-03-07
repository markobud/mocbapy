from Ecosystem import EcosystemModel
from collections import defaultdict
from benpy import solve as bensolve, vlpSolution, vlpProblem
from cobra.util import solver as list_solvers
import optlang
from warnings import warn


def _choose_optlang_interfase(solver=None):
    avail = [sol.lower() for sol, exist in optlang.available_solvers.iteritems() if exist]
    if solver is None:
        warn("Warning: No solver selected. Available solvers: {}".format(str(avail)))
        if len(avail) == 0:
            raise RuntimeError("No solvers available. Try to install glpk or scipy packages")
        warn("Picking {} as solver".format(avail[0]))
        solver = avail[0]
    solver = solver.lower()
    if solver not in avail:
        raise RuntimeError("Solver \'{}\' not available. Available solvers: {}".format(solver, str(avail)))
    return list_solvers.solvers[solver]

def create_model(model_array=None, metabolic_dict=None):
    """Retunrs ans EcosystemModel from parameters"""
    return EcosystemModel(model_array=model_array, metabolic_dict=metabolic_dict)


def get_common_mets(model_list):
    """ Naive implementation of getting common exchange metabolites, using their id's.
        It assumes that identical metabolites have the same id.
    """
    common_mets = defaultdict(dict)
    for model in model_list:
        for rxn_ex in model.exchanges:
            for met_ex in rxn_ex.metabolites:
                common_mets[(met_ex.id, model)] = met_ex.id
    return dict(common_mets)


def mo_fba(ecosystem_model,**kwargs):
    """Solve the Ecosystem Model using bensolve procedure"""
    vlp_eco = ecosystem_model.to_vlp(**kwargs)
    return bensolve(vlp_eco)


def mo_fva(ecosystem_model, fba=None, reactions=None, alpha=0.9,solver=None):
    """Calculate the MO-FVA near the Pareto Front """
    interface = _choose_optlang_interfase(solver)
    model = interface.Model(name='Dummy')
    x1 = interface.Variable("x1", lb=0, ub=20)
    x2 = interface.Variable("x2", lb=0, ub=10)
    c1 = interface.Constraint(2 * x1 - x2, lb=0, ub=0)  # Equality constraint
    model.add([x1, x2, c1])
    model.objective = interface.Objective(x1 + x2, direction="max")
    return model




