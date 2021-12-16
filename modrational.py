import numpy as np
import pandas as pd
from time import time
from datetime import datetime, timedelta
import matplotlib.pyplot as plt



class ModifiedRational():
    def __init__(self, name, c, i, a, tc, event, duration,
                 q=None):
        self.name = name
        self.c = c
        self.i = i
        self.a = a
        self.tc = tc
        self.event = event
        self.duration = duration
        self.q = q

        self.q = round(self.c*(i)*self.a,2)

    def hydrograph(self, start, end, step):

        print (f'Developing time series data for Modified Rational Object "{self.name}"...')
        if self.tc < step.total_seconds()/60 * 2: warnings.warn('Time step is too big - make time step smaller')
        self.results = pd.DataFrame()

        total_time = (end-start).total_seconds()/step.total_seconds()

        self.results['time'] = [start + timedelta(minutes=x) for x in range(int(total_time))]
        self.results['step'] = range(len(self.results['time']))

        slope = self.q/self.tc

        self.results.loc[self.results.index <= self.tc, f'{self.name} flow'] = slope * self.results['step']
        self.results.loc[(self.results.index > self.tc) & (self.results.index <= self.tc + self.duration), f'{self.name} flow'] = self.q
        self.results.loc[(self.results.index > self.tc + self.duration) & (self.results.index <= self.tc*2 + self.duration), f'{self.name} flow'] = -slope * (self.results['step']-(self.tc+self.duration)) + self.q
        self.results.loc[self.results.index > self.tc*2 + self.duration, f'{self.name} flow'] = 0

        # self.results = self.results.drop('step', 1)

        self.results[f'{self.name} step_volume'] = self.results[f'{self.name} flow'] * step.total_seconds() #cubic feet
        self.results[f'{self.name} cum_volume'] = self.results[f'{self.name} step_volume'].cumsum()

        print (f'Developed time series data for Modified Rational Object "{self.name}".\n')

    def total_vol(self):
        return self.results[f'{self.name} cum_volume'].max()



    def plot_hydrograph(self, start):

        plt.figure(1)
        plt.title(f'{self.name} hydrograph')
        plt.xlabel('Time (min)')
        plt.ylabel('Flow (cfs)')
        plt.plot(self.results['step'], self.results[f'{self.name} flow'])
        ax = plt.gca()
        ax.legend()
        # ax.set_xlim([start, start + timedelta(minutes=(self.tc*2+self.duration)*self.duration*.1)])
        ax.set_xlim([0, (self.tc*2+self.duration)+self.duration*.1])

        plt.figure(2)
        plt.title(f'{self.name} volume')
        plt.xlabel('Time (min)')
        plt.ylabel('Volume (cubic feet)')
        plt.plot(self.results['step'], self.results[f'{self.name} cum_volume'])
        ax = plt.gca()
        ax.legend()
        # ax.set_xlim([start, start + timedelta(minutes=(self.tc*2+self.duration)*self.duration*.01)])
        ax.set_xlim([0, (self.tc*2+self.duration)+self.duration*.1])

        plt.show()




events = [[25, 5, 3.28], # 25 year - 5 minute storm
          [25, 10, 2.5], # 25 year - 10 minute storm
          [25, 15, 2.06], # 25 year - 15 minute storm
          [25, 30, 1.39], # 25 year - 30 minute storm
          [25, 60, .86], # 25 year - 1 hour storm
          [25, 360, .204], # 25 year - 6 hour storm
          [25, 24*60, .074], # 25 year - 24 hour storm
          [25, 24*60*2, .041], # 25 year - 2 day storm
          [25, 24*60*7, .018], # 25 year - 7 day storm
          [25, 24*60*10, .014], # 25 year - 10 day storm
         ]

basins = {}

start = datetime(year=2021, day=7, month=11)
end = start + timedelta(days=12)
step = timedelta(minutes=1)

for event in events:
    name = f'{event[0]}-{event[1]}'
    basins[name] = ModifiedRational(name, .6, event[2], 5, 10, event[0], event[1])

for basin in basins:
    basins[basin].hydrograph(start, end, step)
    print (f'Total Runoff Volume for {basins[basin].name} = {basins[basin].total_vol()}')
    basins[basin].plot_hydrograph(start)
