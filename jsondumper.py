from hydro import Rational, DesignPond
from hydro import print_assumptions, print_results, plot_everything
from pondstagestorage import acad_stagestorage_ratingcurve
import json
import pandas as pd



wpond = acad_stagestorage_ratingcurve(r"L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\WestPond_StageStorage.txt")
epond = acad_stagestorage_ratingcurve(r"L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\EastPond_StageStorage.txt")
spond = acad_stagestorage_ratingcurve(r"L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\SouthPond_rev1_StageStorage.txt")
npond = acad_stagestorage_ratingcurve(r"L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\NorthPond_StageStorage.txt")


# epond.append([4153, 12500, 100000])

# print (epond)
icatch = {'R-100': {'Tc': 10, 'C': .69, 'i': 2.51, 'A': 28.44, 'sim': 'pro_25_sim' , 'event': 25},
          'R-200': {'Tc': 10, 'C': .87, 'i': 2.51, 'A': 1.46, 'sim': 'pro_25_sim' , 'event': 25},
          'R-400': {'Tc': 10, 'C': .87, 'i': 2.51, 'A': 1.93, 'sim': 'pro_25_sim' , 'event': 25},
         }

ipond = {'South Pond': {'infil': 2.3,
                   'rat_perv': 1,
                   'us': 'R-100',
                   'pond_curve': spond,
                   'scale': .1
                  },
         'West Pond': {'infil': 3.25,
                            'rat_perv': 1,
                            'us': 'R-200',
                            'pond_curve': wpond,
                            'scale': .2
                  },
         'North Pond': {'infil': 11,
                            'rat_perv': 1,
                            'us': 'R-400',
                            'pond_curve': npond,
                            'scale': 1
                  },
         # 'East Pond': {'infil': 21,
         #                    'rat_perv': 1,
         #                    'us': 'R-300_25',
         #                    'pond_curve': epond
         #          },
        }

isim = {'pro_5_sim': {'max_time': 10000},
        'pro_25_sim' : {'max_time': 10000},
        'pro_100_sim' : {'max_time': 10000}}

with open('catchments.json', 'w') as f:
    json.dump(icatch, f)
with open('ponds.json', 'w') as f:
    json.dump(ipond, f)
with open('simulations.json', 'w') as f:
    json.dump(isim, f)
