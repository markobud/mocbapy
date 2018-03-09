# coding=utf-8
from Ecosystem import EcosystemModel
from collections import defaultdict
from benpy import solve as bensolve, vlpSolution, vlpProblem
from cobra.util import solver as list_solvers
import pandas
import optlang
from warnings import warn
from tqdm import tqdm


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
    """Returns ans EcosystemModel from parameters"""
    return EcosystemModel(model_array=model_array, metabolic_dict=metabolic_dict)


def get_common_mets(model_list):
    """
        Naive implementation of getting common exchange metabolites, using their id's.
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


def mo_fva(base_model = None, fba=None, reactions=None, alpha=0.9,solver=None):
    """
    Calculate the MO-FVA near the Pareto Front

    fba: dictionary {objective_id: obj_value} to set objective constraints

    """
    if solver is None:
        raise RuntimeError('No solver selected')
    if solver == 'gurobi':
        import gurobipy
        if base_model is None:
            raise RuntimeError("base model not given")
        if not isinstance(base_model, gurobipy.Model):
            RuntimeError("base_model is not a gurobi model")
        base_model.update()
        base_model.setParam(gurobipy.GRB.Param.LogToConsole, 0)  # non verbose output
        if fba is None:
            raise RuntimeError("No restrictions given")
        # TODO: How to verify that keys in fba dict are MO-FBA objectives?
        # TODO: Gurobi and context managemnt
        if reactions is None: #all reactions
            reactions = [r.VarName for r in base_model.getVars()]
        fva_res = {rxn:{} for rxn in reactions}
        new_cons = list()
        for obj_id,value in fba.iteritems():
            var = base_model.getVarByName(obj_id)
            new_cons.append(base_model.addConstr(var >= value*alpha))
        base_model.update()
        for dir in ("minimum", "maximum"):
            grb_sense = gurobipy.GRB.MAXIMIZE if dir == "maximum" else gurobipy.GRB.MINIMIZE
            print("Solving {} optimizations".format(dir))
            for rxn in tqdm(reactions):
                flux = base_model.getVarByName(rxn)
                base_model.setObjective(flux, sense=grb_sense)
                base_model.update()
                base_model.optimize()
                fva_res[rxn][dir] = flux.x
        base_model.remove(new_cons)
        base_model.setObjective(0)
        return pandas.DataFrame.from_dict(fva_res, orient='index')

    else:
        raise RuntimeError('solver {} interfase not implemented'.format(solver))

