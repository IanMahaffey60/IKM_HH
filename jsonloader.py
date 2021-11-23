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
    dobj = icatch[catch]
    catchments[catch] = {}
    catchments[catch] = Rational(dobj['C'], dobj['i'], dobj['A'], dobj['Tc'], isim[dobj['sim']]['max_time'])
    catchments[catch].calculate()

# print(catchments)

#------- Initialize Pond Object and Calculate Parameters ---------#

ponds = {}

for pond in ipond:
    pobj = ipond[pond]
    ponds[pond] = {}
    ponds[pond] = DesignPond(catchments[pobj['us']], pobj['pond_curve'], pobj['infil'], pobj['rat_perv'], change=1, calc_vol=True)
    ponds[pond].calculate()

# print (ponds)


#----------- Print and Plot Results ---------#

def results(catch_e, catch_p, p):
    print_assumptions(catchments['Q_e'], catchments['Q_p'], ponds['Pond1'])
    print_results(catchments['Q_e'], catchments['Q_p'], ponds['Pond1'])
    plot_everything(catchments['Q_e'], catchments['Q_p'], ponds['Pond1'])

results(catchments['Q_e'], catchments['Q_p'], ponds['Pond1'])
