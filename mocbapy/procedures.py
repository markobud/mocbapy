from Ecosystem import EcosystemModel

from benpy import solve as bensolve


def create_model(model_array=None, model_dict=None):
    return (EcosystemModel(model_array=model_array, model_dict=model_dict))


def get_common_mets(model_list):
    """Naive implementation of getting common exchange metabolites, using their id's"""
    common_mets = dict()
    for model in model_list:
        for rxn_ex in model.exchanges:
            for met_ex in rxn_ex.metabolites:
                if met_ex.id not in common_mets:
                    common_mets[met_ex.id] = dict([(model, rxn_ex)])
                else:
                    common_mets[met_ex.id][model] = rxn_ex
    return (common_mets)


def mo_fba(ecosystem_model):
    """Solve the Ecosystem Model using bensolve procedure"""
    return bensolve(ecosystem_model.to_vlp())


def mo_fva(ecosystem_model, fba=None, reactions=None, alpha=0.9):
    """Calculate the MO-FVA near the Pareto Front """
    pass
