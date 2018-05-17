# coding=utf-8
import optlang
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

from EcosystemModel import EcosystemModel
from collections import defaultdict
from benpy import solve as bensolve, vlpSolution, vlpProblem
from cobra.util import solver as list_solvers
import pandas
from warnings import warn
from tqdm import tqdm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from scipy.spatial import ConvexHull

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator


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


def mo_fba(ecosystem_model, **kwargs):
    """Solve the Ecosystem Model using bensolve procedure"""
    vlp_eco = ecosystem_model.to_vlp(**kwargs)
    return bensolve(vlp_eco)


def mo_fva(ecosystem_model, fba=None, reactions=None, alpha=0.9, solver=None):
    """Calculate the MO-FVA near the Pareto Front """
    # x1 = interfase.Variable("x1", lb=0, ub=20)
    # x2 = interfase.Variable("x2", lb=0, ub=10)
    # c1 = interfase.Constraint(2 * x1 - x2, lb=0, ub=0)  # Equality constraint
    # model.add([x1, x2, c1])
    # model.objective = interfase.Objective(x1 + x2, direction="max")
    interfase = _choose_optlang_interfase(solver)
    base_model = build_base_opt_model(ecosystem_model, solver=solver)
    base_model.update()
    rxn_dict = {r.name: r for r in base_model.variables}
    if fba is None:
        raise RuntimeError("No MO-FBA restriction were given")
    if reactions is None: #Go for all
        reactions = rxn_dict.keys()
    fva_res = {rxn: {} for rxn in reactions}
    for obj_id, value in fba.iteritems():
        var = rxn_dict[obj_id]
        print("Adding Constraint {} <= {}*{}".format(obj_id, value, alpha))
        base_model.add(interfase.Constraint(var, lb=value*alpha))
    base_model.update()

    for senses in ("minimum", "maximum"):
        synonyms = {"minimum": "min", "maximum": "max"}
        print("Solving {} optimizations".format(senses))
        for rxn in tqdm(reactions):
            flux = rxn_dict[rxn]
            base_model.objective = interfase.Objective(flux, direction=synonyms[senses])
            base_model.update()
            base_model.optimize()
            fva_res[rxn][senses] = flux.primal

    return pandas.DataFrame.from_dict(fva_res, orient='index')


def sum_from_list(list_expr):
    """ Dichotomous construction of expressions"""
    def sum_from_list_p(le, a, b):
        if a == b:
            return 0
        else:
            if b - a == 1:
                return le[a]
            else:
                middle = int((a+b)/2)
                return sum_from_list_p(le, a, middle) + sum_from_list_p(le, middle, b)

    return sum_from_list_p(list_expr, 0, len(list_expr))


def build_base_opt_model(ecomodel, solver=None):
    """ Builds the underlying base optimization problem Sv = 0, lb <= v <= ub """
    interfase = _choose_optlang_interfase(solver)
    model = interfase.Model(name='Base Solver Model')
    m, n = ecomodel.Ssigma.shape
    assert m == len(ecomodel.sysmetabolites)
    assert n == len(ecomodel.sysreactions)
    # Create flux variables
    flux_variables = [interfase.Variable(rxn, lb=ecomodel.lb[i], ub=ecomodel.ub[i]) for i, rxn in
                      enumerate(ecomodel.sysreactions)]
    model.add(flux_variables, sloppy=True)
    model.update()
    for i in tqdm(range(m)):
        terms_const = [flux_variables[j] * ecomodel.Ssigma[i, j] for j in range(n) if ecomodel.Ssigma[i, j] != 0]
        mass_const = interfase.Constraint(sum_from_list(terms_const), lb=0, ub=0)
        model.add(mass_const, sloppy=True)
    model.update()
    model.objective = interfase.Objective(0, direction="max")
    return model


def draw3d(polygon):
    """ Returns a 3D figure of the Pareto Front. Input: A mo-fba Polygon (eg. sol_mofba.Primal) """
    points = polygon.vertex_value[[x == 1 for x in polygon.vertex_type]]
    hull = ConvexHull(points)
    pd = Poly3DCollection([hull.points[s] for s in hull.simplices])
    fig = plt.figure(figsize=(9, 10))
    ax = fig.add_subplot(111, projection='3d')
    pd.set_facecolor('yellow')
    pd.set_alpha(0.4)
    pd.set_edgecolor('black')

    ax.add_collection3d(pd)
    return fig, ax


