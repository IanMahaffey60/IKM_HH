from hydro import Rational, DesignPond
from hydro import print_assumptions, print_results, plot_everything


#------------- INPUTs -------------#

infil = 1 #inch/hour
rat_perv = 70664/270894 #Ratio of pervious to impervious area
i_25yr = 2.51 #inches/hour
Tc_e = 10 #minutes
Tc_p = 10 #minutes
C_p = .9 #unitless
C_e = .3 #unitless
A_p = 36.44 #acres
A_e = 36.44 #acres
pond_curve = [[4133, 145035, 0], [4134,154652, 48977], [4135, 164386, 101624], [4136, 174201, 270894], [4137, 174201, 443829]] #Full Size Pond
max_time = 300 #minutes

#-------- Initialize Rational Objects and Calculate Parameters ------#

Q_p = Rational(C_p, i_25yr, A_p, Tc_p, max_time)
Q_p.calculate()
Q_e = Rational(C_e, i_25yr, A_e, Tc_e, max_time)
Q_e.calculate()

#------- Initialize Pond Object and Calculate Parameters ---------#

pond = DesignPond(Q_p, pond_curve, infil, rat_perv, change=1, calc_vol=True)
pond.calculate()

#----------- Print and Plot Results ---------#

print_assumptions(Q_e, Q_p, pond)
print_results(Q_e, Q_p, pond)
plot_everything(Q_e, Q_p, pond)
