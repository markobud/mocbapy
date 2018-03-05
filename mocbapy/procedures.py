from Ecosystem import EcosystemModel
from collections import defaultdict
from benpy import solve as bensolve


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


def mo_fba(ecosystem_model):
    """Solve the Ecosystem Model using bensolve procedure"""
    return bensolve(ecosystem_model.to_vlp())


def mo_fva(ecosystem_model, fba=None, reactions=None, alpha=0.9):
    """Calculate the MO-FVA near the Pareto Front """
    pass
