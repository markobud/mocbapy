# Paretogrids test file
import cobra.test
import mocbapy
from Paretogrids import *

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

mygrid = ParetoGrid(emodel = test_EcoSys, sol_mofba=sol_mofba,modelIDs=['model1','model2'], step=30j)
mygrid.randomNormal2grid(n=2)
mygrid.plotGrid(background="white",fontsize=12,figsize=(8,4),cmap='coolwarm')
