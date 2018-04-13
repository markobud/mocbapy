# mocbapy test file
import cobra.test
import mocbapy
from mocbapy.utilities import build_base_opt_model
test_arr = list()
n_test = 2
for i in range(n_test):
    model = cobra.test.create_test_model("ecoli")
    model.id = model.id + '_' + str(i+1)
    test_arr.append(model)

common_mets = mocbapy.get_common_mets(test_arr)
test_EcoSys = mocbapy.create_model(model_array=test_arr, metabolic_dict=common_mets)

index_gluc_pool = test_EcoSys.sysreactions.index("EX_glc__D_e:pool")
test_EcoSys.lb[index_gluc_pool] = -5

bensolve_opts = mocbapy.bensolve_default_options
bensolve_opts['message_level'] = 0
sol_mofba = mocbapy.mo_fba(test_EcoSys, options=bensolve_opts)
print(sol_mofba)

#Objective reactions:
#test_EcoSys.sysreactions[7]
#test_EcoSys.sysreactions[2590]

base_problem = build_base_opt_model(test_EcoSys, solver='glpk')
#Pick up one extrme point
pext_3=sol_mofba.Primal.vertex_value[3]
fva_dict = {test_EcoSys.sysreactions[7]:pext_3[0], test_EcoSys.sysreactions[2590]:pext_3[1]}
fva_res_pext3 = mocbapy.mo_fva(test_EcoSys, fba=fva_dict, solver='gurobi')

pext_4=sol_mofba.Primal.vertex_value[4]
fva_dict = {test_EcoSys.sysreactions[7]:pext_4[0], test_EcoSys.sysreactions[2590]:pext_4[1]}
fva_res_pext4 = mocbapy.mo_fva(test_EcoSys, fba=fva_dict, solver='glpk')
