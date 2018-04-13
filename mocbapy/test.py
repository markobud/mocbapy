#mocbapy test file
import cobra.test
import mocbapy

test_arr = list()
n_test = 3
for i in range(n_test):
    model = cobra.test.create_test_model("ecoli")
    model.id = model.id + str(i)
    test_arr.append(model)

common_mets = mocbapy.get_common_mets(test_arr)
test_EcoSys = mocbapy.create_model(model_array=test_arr, metabolic_dict=common_mets)
#bensolve_opts = test_EcoSys.bensolve_default_options
#bensolve_opts['message_level'] = 0
sol_mofba = mocbapy.mo_fba(test_EcoSys)
#sol_mofva = mocbapy.mo_fva(test_EcoSys, solver='gurobi')

