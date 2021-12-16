import numpy as np
import pandas as pd


class ModifiedRational():
    def __init__(self, name, c, i, a, tc, event=None, duration=None,
                 q=None):
        self.name = name
        self.c = c
        self.i = i
        self.a = a
        self.tc = tc
        self.event = event
        self.duration = duration
        self.q = q

        self.q = round(self.c*(self.duration)*self.a,2)

    def hydrograph(self, start, end, step):

        print (f'Developing time series data for Rational Object "{self.name}"...')
        if self.tc < step.total_seconds()/60 * 2: warnings.warn('Time step is too big - make time step smaller')
        self.results = pd.DataFrame()

        total_time = (end-start).total_seconds()/step.total_seconds()

        self.results['time'] = [start + timedelta(minutes=x) for x in range(int(total_time))]
        self.results['step'] = range(len(self.results['time']))

        slope = self.q/self.tc
        self.results.loc[self.results.index <= self.tc, f'{self.name} flow'] = slope * self.results['step']
        self.results.loc[(self.results.index > self.tc) & (self.results.index <= self.tc * 2), f'{self.name} flow'] = -slope * self.results['step'] + 2 * self.q
        self.results.loc[self.results.index > self.tc * 2, f'{self.name} flow'] = 0

        self.results = self.results.drop('step', 1)

        self.results[f'{self.name} step_volume'] = self.results[f'{self.name} flow'] * step.total_seconds() #cubic feet
        self.results[f'{self.name} cum_volume'] = self.results[f'{self.name} step_volume'].cumsum()

        print (f'Developed time series data for Rational Object "{self.name}".\n')



        # totaltime = 2 *self.tc + self.duration
        # stormduration = self.tc + self.duration
        #
        # self.hydrograph = pd.DataFrame()

        for time in np.arange(0, totaltime, 1):
            if time <= self.tc:
                # 0 - c*i(for xxx) * a


            elif time > self.tc and time <= self.tc + self.duration:
                # c*i(for xxx) * a - c*i(for xxx) * a

            elif time > self.tc + self.duration:
                # c*i(for xxx) * a - 0


    def total_vol(self):
        return self.hydrograph['flow'].sum



events = [[25, 5, 3.28], # 25 year - 5 minute storm
          [25, 10, 2.5], # 25 year - 10 minute storm
          [25, 15, 2.06], # 25 year - 15 minute storm
          [25, 30, 1.39], # 25 year - 30 minute storm
          [25, 60, .86], # 25 year - 1 hour storm
          [25, 360, .204], # 25 year - 6 hour storm
          [25, 24*60, .074] # 25 year - 24 hour storm
         ]

basins = {}

start = datetime(year=2021, day=7, month=11)
end = start + timedelta(days=1)
step = timedelta(minutes=1)

for event in events:
    name = f'{event[0]}-{event[1]}'
    basins[name] = ModifiedRational(name, .6, event[2], 10, event[0], event[1])

for basin in basins:
    basins[basin].hydrograph(start, end, step)
    print (basins[basin].total_vol)
