from hydro import Rational, DesignPond
from hydro import print_assumptions, print_results, plot_everything
import json



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
    catchments[catch] = Rational(dobj['C'], dobj['i'], dobj['A'], dobj['Tc'], isim[dobj['sim']]['max_time'])
    catchments[catch].calculate()

# print(catchments)

#------- Initialize Pond Object and Calculate Parameters ---------#

ponds = {}

for pond in ipond:
    print (f'---- Creating {pond} Pond ----')
    pobj = ipond[pond]
    ponds[pond] = {}
    ponds[pond] = DesignPond(catchments[pobj['us']], pobj['pond_curve'], pobj['infil'], pobj['rat_perv'], change=1, calc_vol=True)
    ponds[pond].calculate()



#----------- Print and Plot Results ---------#

def results(catch_e, catch_p, p):
    print_assumptions(catch_e, catch_p, p)
    print_results(catch_e, catch_p, p)
    plot_everything(catch_e, catch_p, p)

results(catchments['R-100_25'], catchments['R-100_25'], ponds['South Pond'])
results(catchments['R-200_25'], catchments['R-200_25'], ponds['West Pond'])
results(catchments['R-300_25'], catchments['R-300_25'], ponds['East Pond'])
results(catchments['R-400_25'], catchments['R-400_25'], ponds['North Pond'])
