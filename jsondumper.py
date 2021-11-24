from hydro import Rational, DesignPond
from hydro import print_assumptions, print_results, plot_everything
from pondstagestorage import acad_stagestorage_ratingcurve
import json
import pandas as pd



wpond = acad_stagestorage_ratingcurve(r"L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\WestPond_StageStorage.txt")
epond = acad_stagestorage_ratingcurve(r"L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\EastPond_StageStorage.txt")
spond = acad_stagestorage_ratingcurve(r"L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\SouthPond_StageStorage.txt")
npond = acad_stagestorage_ratingcurve(r"L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\NorthPond_StageStorage.txt")


# epond.append([4153, 12500, 100000])

# print (epond)
icatch = {'R-100_5': {'Tc': 10, 'C': .85, 'i': 1.56, 'A': 28.44, 'sim': 'pro_5_sim', 'event': 5},
          'R-100_25': {'Tc': 10, 'C': .9, 'i': 2.51, 'A': 28.44, 'sim': 'pro_25_sim' , 'event': 25},
          'R-100_100': {'Tc': 10, 'C': .9, 'i': 3.68, 'A': 28.44, 'sim': 'pro_100_sim' , 'event': 100},
          'R-200_5': {'Tc': 10, 'C': .85, 'i': 1.56, 'A': 1.46, 'sim': 'pro_5_sim', 'event': 5},
          'R-200_25': {'Tc': 10, 'C': .9, 'i': 2.51, 'A': 1.46, 'sim': 'pro_25_sim' , 'event': 25},
          'R-200_100': {'Tc': 10, 'C': .9, 'i': 3.68, 'A': 1.46, 'sim': 'pro_100_sim' , 'event': 100},
          'R-300_5': {'Tc': 10, 'C': .85, 'i': 1.56, 'A': 12.39, 'sim': 'pro_5_sim', 'event': 5},
          'R-300_25': {'Tc': 10, 'C': .9, 'i': 2.51, 'A': 2.79, 'sim': 'pro_25_sim' , 'event': 25},
          'R-300_100': {'Tc': 10, 'C': .9, 'i': 3.68, 'A': 12.39, 'sim': 'pro_100_sim' , 'event': 100},
          'R-400_5': {'Tc': 10, 'C': .85, 'i': 1.56, 'A': 1.93, 'sim': 'pro_5_sim', 'event': 5},
          'R-400_25': {'Tc': 10, 'C': .9, 'i': 2.51, 'A': 11.53, 'sim': 'pro_25_sim' , 'event': 25},
          'R-400_100': {'Tc': 10, 'C': .9, 'i': 3.68, 'A': 1.93, 'sim': 'pro_100_sim' , 'event': 100},
         }

ipond = {'South Pond': {'infil': 1,
                   'rat_perv': 1,
                   'us': 'R-100_25',
                   'pond_curve': spond
                  },
         'West Pond': {'infil': 1,
                            'rat_perv': 1,
                            'us': 'R-200_25',
                            'pond_curve': wpond
                  },
         'North Pond': {'infil': 1,
                            'rat_perv': 1,
                            'us': 'R-400_25',
                            'pond_curve': npond
                  },
         'East Pond': {'infil': 2,
                            'rat_perv': 1,
                            'us': 'R-300_25',
                            'pond_curve': epond
                  },
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
