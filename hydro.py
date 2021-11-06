#!/usr/bin/env python3

''' Detention/Retention pond sizing and outlet sizing.

Current Project:
    1. Switching lists to dataframes for fast implementation

Current Limitations:
    1. Manual data entry
    2. Plenty others

Future Functionality:
    1. Access actual flow series of each individual outlet
    2. Sedimentation/water quality design parameters
    3. Automatic design tools (change pond curve to meet requirements, design outlets)
    4. Add place for external time series, flow series, and volume series input
    5. Calculate SCS method instead of just Rational
    6. Separate outlet flows to be able to analyze each indiviual performance
    7. Spatial incorporation
        a. Calculate Tc
        b. Find Rainfall Intensity based on Tc from NOAA Rasters
        c. Extract other hydrologic and hydraulic properties from GIS
    8. Utilize QT Designer to create GUI for inputting data and visualizing results
    9. Take input of pond rating curve from acad report or manual entry
    10. Import and export to hdf5 file format. Allow savable projects to create, edit, save, and load.


Author: Ian Mahaffey, P.E.
Date: 10/26/2021
Revision: 3

'''
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import warnings
import pandas as pd
from time import time




def timer_func(func):
    # This function shows the execution time of
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f'Function {func.__name__!r} executed in {(t2-t1):.4f}s')
        return result
    return wrap_func

def acad_stagestorage_ratingcurve(ifile):
    '''
    Creates a dataframe from the stage storage results file output from ACAD.
    '''
    col= ['cont elev', 'cont area', 'depth', 'incr vol1', 'cum vol1', 'incr vol2', 'incr vol2']
    df = pd.read_csv(ifc, skiprows=10, sep="\s{2}", header = None, thousands=',', engine='python')
    df.columns = col
    return df


class Simulation():
    def __init__(self, start, end, step):
        self.start = start
        self.end = end
        self.step = step

    def foo(self):
        pass


class Rational():
    '''
    Class object that creates and calculates all parameters associated with the rational hydrologic method
    Inputs:
        1. C-Value (c)
        2. Rainfall Intensity (i)
        3. Area (a)
        4. Time of Concentration (tc)
        5. Max time of analysis (max_time)

    '''
    def __init__(self, c, i, a, tc, max_time, q=None,results=None):
        self.c = c #Rational Roughness Coefficient (Unitless)
        self.i = i #Rainfall Intensity (Inches/Hour)
        self.a = a #Area (Acres)
        self.tc = tc #Time of Concentraion (Minutes)
        self.max_time = max_time #Full length of analysis (Minutes)
        self.results = results
        self.q = q #Peak Flow (cfs)
    @timer_func
    def calculate(self):
        # np.arange(0,self.max_time,1
        self.results = pd.DataFrame(columns=['time', 'flow', 'volume'])

        self.q = round(self.c*self.i*self.a,2)

        #THERES A TON OF POTENTIAL TO SPEED THIS UP - DOESN'T NEED TO BE A LOOP, MAYBE SOME SORT OF IF/THEN STATEMENT
        for i in range(0, self.tc*2+1, 1):
            rat1 = i/self.tc
            rat2 = (self.tc-(i-self.tc))/self.tc
            self.flowvol_series_append(rat1, i) if i <= self.tc else self.flowvol_series_append(rat2, i)
        for i in range(self.tc*2+1, self.max_time-2*self.tc+self.tc*2+1):
            self.results.loc[len(self.results.index)] = [i, 0, 0]
        return self.q
    @timer_func
    def total_vol(self):
        return round(self.results.volume.sum(),2)
    @timer_func
    def flowvol_series_append(self, rat, t):
        self.results.loc[len(self.results.index)] = [t, rat*self.q, rat*self.q*60]

class DesignPond():
    '''
    A Class Object that creates and calculates all the required parameters for a detention/retention pond

    Future functionality:
        1. Add orifice, spillway, weir outlet types
        2. Sedimentation/water quality design parameters
        3. Automatic design tools (change pond curve to meet requirements, design outlets)

    Inputs:
        1. Hydrologic Basin Object (Rational Object)
        2. Pond Curve (pond_curve) - Required format must be List of Lists [elevation, contour area, cumulative volume]
        3. Infiltration Rate (infil)
        4. Ratio of max pond area pervious (rat_perv)
        5. Scale pond up or down based on ratio (change)
        6. Calculate the volume of the pond (calc_vol) - Sometimes this may not be desired if pond volume was calculated elsewhere
    '''


    def __init__(self, qin, pc, infil, rat_perv, change=False, calc_vol=False,
            rating_curves=None, results=None, time_to_empty=None,
                 find_footprint=None, find_vol_fromElev=None, find_elev_fromVol=None,
                 outlet=None):
        self.qin = qin
        self.pond_curve = pc #dataframe with all hydraulic rating curves, including outlets
        self.infil = infil #Infiltration Rate (Inches / Hour)
        self.rat_perv = rat_perv #Ratio of pervious vs. impervious (decimal)
        self.change = change #Scale pond for preliminary design purposes (Default: False)
        self.calc_vol = calc_vol #Calculate pond volume (Default: True)
        self.rating_curves = rating_curves
        self.results = results
        self.time_to_empty = time_to_empty
        self.find_footprint = find_footprint
        self.find_vol_fromElev = find_vol_fromElev
        self.find_elev_fromVol = find_elev_fromVol
        self.outlet = outlet
        self.pond_curve = pd.DataFrame(self.pond_curve, columns=['elev', 'footprint', 'volume'])
        if self.change != False:
            self.scale_pond()
        if self.calc_vol == True:
            self.calc_volume()
        self.interp_pond()
        self.pond_ratingcurve(None)
    @timer_func
    def add_outlet(self, type, **input):
        if not self.outlet:
            self.outlet = {}
        if type == 'Orofice':
            self.outlet[input['name']] = Orofice(input['name'], input['min_elev'], input['max_elev'], input['shape'], input['diameter'])
        elif type == 'Weir':
            self.outlet[input['name']] = Weir(input['name'], input['min_elev'], input['max_elev'], input['length'], input['coeff'])
        self.pond_curve['elev'] = round(self.pond_curve['elev'],2)
        self.pond_ratingcurve(self.outlet[input['name']])
    @timer_func
    def interp_pond(self):
        newp = pd.DataFrame(np.arange(self.pond_curve['elev'].min(), self.pond_curve['elev'].max(), .01))
        newp.columns = ['elev']
        newp['footprint'] = interp1d(self.pond_curve['elev'], self.pond_curve['footprint'])(newp['elev'])
        newp['volume'] = interp1d(self.pond_curve['elev'], self.pond_curve['volume'])(newp['elev'])
        self.pond_curve = newp
        return self.pond_curve
    @timer_func
    def pond_ratingcurve(self, outlet):
        '''
        os = self.outlet dictionary
        '''
        if 'infiltration' not in self.pond_curve:
            self.pond_curve['infiltration'] = self.pond_curve['footprint'] * self.infil/12/60 * self.rat_perv
            self.pond_curve['overall outflow'] = self.pond_curve['infiltration']
        if outlet is not None:
            df = outlet.results
            self.pond_curve = self.pond_curve.merge(df[['elev', 'flow']], on='elev', how='left')
            self.pond_curve.rename(columns={'flow': outlet.name}, inplace=True)
            self.pond_curve[outlet.name] = self.pond_curve[outlet.name].fillna(0)
            self.pond_curve['overall outflow'] = self.pond_curve[outlet.name] + self.pond_curve['overall outflow']
        return self.pond_curve
    # @timer_func
    # def calculate(self):
    #     elev_bottom = min(self.pond_curve['elev'])
    #     elev = min(self.pond_curve['elev'])
    #
    #     for num, qi in enumerate(self.qin.flow_series):
    #         in_vol = qi*60 #turns to cubic feet
    #         infil_vol = self.infil/12/60*self.find_footprint(elev)*self.rat_perv
    #         self.infil_series.append(self.infil/12/60*self.find_footprint(elev)*self.rat_perv)
    #         out_vol = infil_vol
    #         if self.outlet:
    #             for out in self.outlet:
    #                 out_vol += self.get_flow(elev, self.outlet[out])
    #                 self.outlet_series[out] = self.get_flow(elev, self.outlet[out])
    #         if out_vol > (self.vol_series[-1]+in_vol):
    #             out_vol = self.vol_series[-1]+in_vol
    #         self.vol_in_series.append(in_vol)
    #         self.vol_out_series.append(out_vol)
    #         if in_vol-out_vol + self.vol_series[-1] < 0:
    #             self.vol_series.append(vol[-1]+0)
    #         else:
    #             self.vol_series.append(self.vol_series[-1]+in_vol-out_vol)
    #         elev = np.round(self.find_elev_fromVol(self.vol_series[-1]),2)
    #         self.time_series.append(self.qin.time_series[num])
    #         self.elevation_series.append(elev)
    #     self.pond_empty_time()
    @timer_func
    def calc_volume(self):
        last_vol = []
        last_sqft = []
        for num, i in enumerate(self.pond_curve):
            if num >0:
                last_sqft.append(i[1])
                last_vol.append(last_vol[-1]+(last_sqft[-1]+i[1])/2)
            else:
                last_vol.append(0)
        for num, i in enumerate(self.pond_curve):
            i[2] = last_vol[num]
        return self.pond_curve
    @timer_func
    # def pond_empty_time(self):
    #     val = -1
    #     if self.vol_series is not None:
    #         for num, vol in enumerate(self.vol_series):
    #             if vol < 1 and vol < val:
    #                 self.time_to_empty = self.time_series[num]
    #             val = vol
    #     else:
    #         warnings.warn('Volume series has not been calculated yet - consider executing the "calculate" function first')
    #     if self.time_to_empty == None:
    #         warnings.warn('Time to Empty returned: None, try increasing total analysis time')
    @timer_func
    def scale_pond(self, changer):
        self.pond_curve['footprint'] = self.pond_curve['footprint']*changer
        self.pond_curve['volume'] = self.pond_curve['volume']*changer
        self.pond_ratingcurve(None)
        return self.pond_curve

class Outlet():
    def __init__(self, name, low_elev, max_elev, results=None):
        self.name = name
        self.low_elev = low_elev
        self.max_elev = max_elev
        self.results = results
    @timer_func
    def calc_rating_curve(self):
        self.results = pd.DataFrame(columns=['elev', 'flow'])
        self.results['elev'] = np.round(np.arange(self.low_elev, self.max_elev, .01),2)
        self.results['flow'] = self.calc_flow(self.results['elev'])
        return self.results

class Orofice(Outlet):
    def __init__(self, name, low_elev, max_elev, shape, diameter, coeff=None, results=None):
        super().__init__(name, low_elev, max_elev, results)
        self.shape = shape
        self.diameter = diameter
        self.coeff = coeff
        self.calc_rating_curve()
    @timer_func
    def calc_coeff_of_discharge(self):
        if self.shape == 'Sharp Orifice':
            self.coeff = 0.62
            return 0.62
        elif self.shape == 'Tube':
            self.coeff = 0.80
            return 0.80
    @timer_func
    def calc_flow(self, elev):
        q = self.calc_coeff_of_discharge() * (np.pi*(self.diameter/12)**2/4) * np.sqrt(2*9.81*(elev-self.low_elev))
        return round(q,2)

class Weir(Outlet):
    def __init__(self, name, low_elev, max_elev, length, coeff, results=None):
        super().__init__(name, low_elev, max_elev)
        self.coeff = coeff
        self.length = length
        self.calc_rating_curve()
    @timer_func
    def calc_flow(self, elev):
        q = self.coeff * self.length * (elev - self.low_elev+.003)**(3/2)
        return round(q,2)


def print_assumptions(q_e, q_p, pond):
    print ('\n------ Assumptions ------\n')

    print (f'Infiltration rate = {pond.infil} inches/hour')
    print (f'Existing Rainfall Intensity (10-min Tc) = {q_e.i} inches/hour')
    print (f'Proposed Rainfall Intensity (10-min Tc) = {q_p.i} inches/hour')
    print (f'Existing Time of Concentration = {q_e.tc} minutes')
    print (f'Proposed Time of Concentration = {q_p.tc} minutes')
    print (f'Existing Area = {q_e.a} acres')
    print (f'Proposed Area = {q_p.a} acres')

def print_results(q_e, q_p, pond):
    print ('\n-------- Results ---------\n')

    print (f'Existing Peak Flow = {q_e.q} cfs')
    print (f'Proposed Peak Flow = {q_p.q} cfs')
    print (f'Existing Total Volume = {q_e.total_vol()} cu-ft')
    print (f'Proposed Total Volume = {q_p.total_vol()} cu-ft')
    print (f'Required Volume Retained = {q_p.total_vol()-q_e.total_vol()} cu-ft')
    print (f'Max Infiltration Volume Rate = {round(max(pond.vol_out_series),2)} cu-ft / minute')
    print (f'Max pond footprint = {np.round(pond.find_footprint(max(pond.elevation_series)),2)} sq-ft')
    print (f'Total time till pond empties = {pond.time_to_empty} minutes ({round(pond.time_to_empty/60,2)} hours)')
    print (f'Minimum elevation of pond = {pond.elevation_series[0]} feet')
    print (f'Max water surface elevation in pond = {np.round(max(pond.elevation_series),2)} feet')
    print (f'Max Depth in pond = {np.round(max(pond.elevation_series)-pond.elevation_series[0],2)} feet')
    print (f'Max volume in pond = {np.round(max(pond.vol_series))} cu-ft')
    print (f'Total volume retained = {round(sum(pond.vol_out_series,2))} cu-ft')

def plot_everything(q_e, q_p, pond):
    plt.figure(1)
    plt.title('Existing vs. Proposed Flow')
    plt.xlabel('Time (min)')
    plt.ylabel('Flow (cfs)')
    plt.plot(q_p.time_series, q_p.flow_series, label='Proposed')
    plt.plot(q_e.time_series, q_e.flow_series, label='Existing')
    plt.legend()
    ax = plt.gca()
    ax.set_xlim([0, max([q_e.tc*2, q_p.tc*2])])

    plt.figure(2)
    plt.title('Pond Volume Rating Curve')
    plt.xlabel('Volume (ft^3)')
    plt.ylabel('Elevation (ft)')
    plt.plot(pond.vol_curve, pond.elev_curve)

    plt.figure(3)
    plt.title('Pond Footprint Rating Curve')
    plt.xlabel('Footprint (ft^2)')
    plt.ylabel('Elevation (ft)')
    plt.plot(pond.footprint_curve, pond.elev_curve)
    plt.figure(4)
    plt.title('Pond Elevation vs Time')
    plt.xlabel('Time (min)')
    plt.ylabel('Elevation (ft)')
    plt.ticklabel_format(useOffset=False)
    plt.plot(pond.time_series,pond.elevation_series)
    ax = plt.gca()
    ax.set_xlim([0, pond.time_to_empty+5])

    plt.figure(5)
    plt.title('Inflow Volume vs Outflow Volume')
    plt.xlabel('Time (min)')
    plt.ylabel('Volume (ft^3)')
    plt.plot(pond.time_series,pond.vol_in_series, label='Inflow')
    plt.plot(pond.time_series,pond.vol_out_series, label='Outflow')
    plt.legend()
    ax = plt.gca()
    ax.set_xlim([0, pond.time_to_empty+5])

    plt.figure(6)
    plt.title('Pond Volume vs Time')
    plt.xlabel('Time (min)')
    plt.ylabel('Volume (ft^3)')
    plt.plot(pond.time_series, pond.vol_series, label='Pond Volume')
    ax = plt.gca()
    ax.set_xlim([0, pond.time_to_empty+5])
    plt.legend()

    plt.show()


'''

Testing Pieces


'''

#------------- INPUTs -------------#

infil = 1 #inch/hour
rat_perv = .5 #Ratio of pervious to impervious area
i_25yr = 2.51 #inches/hour
Tc_e = 10 #minutes
Tc_p = 10 #minutes
C_p = .9 #unitless
C_e = .3 #unitless
A_p = 36.44 #acres
A_e = 36.44 #acres
pond_curve = [[4133, 145035, 0], [4134,154652, 48977], [4135, 164386, 101624], [4136, 174201, 270894], [4137, 174201, 443829]] #Full Size Pond
max_time = 40 #minutes

#-------- Initialize Rational Objects and Calculate Parameters ------#

Q_p = Rational(C_p, i_25yr, A_p, Tc_p, max_time)
Q_p.calculate()
Q_e = Rational(C_e, i_25yr, A_e, Tc_e, max_time)
Q_e.calculate()

#------- Initialize Pond Object and Calculate Parameters ---------#

pond = DesignPond(Q_p, pond_curve, infil, rat_perv)
pond.scale_pond(.05)
# pond.calc_volume()
pond.add_outlet('Orofice', name='o1', min_elev=4133, max_elev=4137, shape='Sharp Orifice', diameter=8)
pond.add_outlet('Orofice', name='o2', min_elev=4134, max_elev=4137, shape='Sharp Orifice', diameter=6)
pond.add_outlet('Weir', name='w1', min_elev=4133, max_elev=4137, length=4, coeff=1.3)


# #----------- Print and Plot Results ---------#
#
# print_assumptions(Q_e, Q_p, pond)
# print_results(Q_e, Q_p, pond)
# plot_everything(Q_e, Q_p, pond)
