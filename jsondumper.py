from hydro import Rational, DesignPond
from hydro import print_assumptions, print_results, plot_everything
import json


icatch = {'Q_p': {'Tc': 10, 'C': .9, 'i': 2.51, 'A':36.44, 'sim': 'pro_sim', 'event': 25},
          'Q_e': {'Tc': 10, 'C': .3, 'i': 2.51, 'A':36.44, 'sim': 'ex_sim' , 'event': 25}
         }
ipond = {'Pond1': {'infil': 1,
                   'rat_perv': 70664/270894,
                   'us': 'Q_p',
                   'pond_curve': [[4133, 145035, 0],
                                  [4134,154652, 48977],
                                  [4135, 164386, 101624],
                                  [4136, 174201, 270894],
                                  [4137, 174201, 443829]
                                  ]
                  }
        }

isim = {'pro_sim': {'max_time': 4800},
        'ex_sim' : {'max_time': 4800}}

with open('catchments.json', 'w') as f:
    json.dump(icatch, f)
with open('ponds.json', 'w') as f:
    json.dump(ipond, f)
with open('simulations.json', 'w') as f:
    json.dump(isim, f)
