from hydro import Rational, DesignPond, Simulation
from hydro import print_assumptions, print_results, plot_everything
import json

#----------- Initiate Simulation ------------#

sim1 = Simulation()
sim1.createProject()

#------------- INPUTs -------------#


with open('catchments.json', 'r') as f:
    icatch = json.load(f)
with open('ponds.json', 'r') as f:
    ipond = json.load(f)
with open('simulations.json', 'r') as f:
    isim = json.load(f)



#-------- Initialize Rational Objects and Calculate Parameters ------#

catchments = {}

for catch in icatch:
    print (f'---- Creating {catch} Catchment ----')
    dobj = icatch[catch]
    catchments[catch] = {}
    catchments[catch] = Rational(catch, dobj['C'], dobj['i'], dobj['A'], dobj['Tc'], isim[dobj['sim']]['max_time'])
    catchments[catch].calculate()

# print(catchments)

#------- Initialize Pond Object and Calculate Parameters ---------#

ponds = {}

for pond in ipond:
    print (f'---- Creating {pond} Pond ----')
    pobj = ipond[pond]
    ponds[pond] = {}
    ponds[pond] = DesignPond(pond, catchments[pobj['us']], pobj['pond_curve'], pobj['infil'], pobj['rat_perv'], change=pobj['scale'], calc_vol=True)
    ponds[pond].calculate()





#----------- Print and Plot Results ---------#




def results(catch_e, catch_p, p):
    print_assumptions(catch_e, catch_p, p)
    print_results(catch_e, catch_p, p)
    plot_everything(catch_e, catch_p, p)

results(catchments['R-100'], catchments['R-100'], ponds['South Pond'])
results(catchments['R-200'], catchments['R-200'], ponds['West Pond'])
results(catchments['R-400'], catchments['R-400'], ponds['North Pond'])
