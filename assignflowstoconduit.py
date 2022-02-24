import pandas as pd
import numpy as np
import json, math, csv
from datetime import datetime

def logger(iFC):
    with open('log.txt', 'a') as f:
        now = datetime.now()
        txt = '\n' + now.strftime("%d/%m/%Y %H:%M:%S") + ': ' + str(iFC) + '\n'
        f.write(txt)

standardsizes = [6,8,10,12,15,18,21,24,30,36,42,48,60]
class Conduit():
    def __init__(self, id, shape=None, mat=None, D=None, L=None,
               x1=None, x2=None, y1=None, y2=None, B=None, z=None,
               S=None, IE_up=None, IE_down=None, rough=None, Q=None,
               v=None, depth=None, size=False, rat_full=0.94):
        self.id = id
        self.shape = shape
        self.mat = mat
        self.D = D #inches
        self.L = L #feet
        self.x1 = x1 #feet
        self.x2 = x2  #feet
        self.y1 = y1 #feet
        self.y2 = y2 #feet
        self.B = B #feet
        self.z = z #ratio H:Z
        self.S = S #ft/ft
        self.IE_up = IE_up #feet
        self.IE_down = IE_down #feet
        self.rough = rough
        self.Q = Q #cfs
        self.v = v #fps
        self.depth = depth  #feet
        self.rat_full = rat_full #Ratio full - 94% gets highest allowable flow
        self.size = size
        self.D_orig = self.D


        self.depth = {}
        self.v = {}
    def Calc_Qcap(self):
        if self.S == None:
            self.S = abs((self.IE_up - self.IE_down)/self.L)
        self.D_ft = self.D/12
        self.r_ft = self.D/12/2
        self.theta = 2*math.acos((self.r_ft-self.D_ft*(1-self.rat_full))/self.r_ft)
        self.A = np.pi*self.r_ft**2-self.r_ft**2*(self.theta-math.sin(self.theta))/2
        self.pw = 2*np.pi*self.r_ft-self.r_ft*(self.theta)
        # logger([self.S])

        self.Q_cap = round(1.49/self.rough*self.A*(self.A/self.pw)**(2/3)*self.S**0.5,2)
        # print (f'Calculating Q_cap = {self.Q_cap}')
        if self.D == self.D_orig:
            self.Q_cap_orig = self.Q_cap

    def Calc_depth(self, storm):
        self.D_ft = self.D/12
        self.r_ft = self.D/12/2
        loops = np.linspace(0, float(self.D_ft)-0.01, int(self.D_ft/0.01))
        for i in loops:
            # print (i)
            if i == 0:
                continue
            theta = 2*math.acos((self.r_ft-self.D_ft*(1-i/self.D_ft))/self.r_ft)
            area = np.pi*self.r_ft**2-self.r_ft**2*(theta-math.sin(theta))/2
            per=2*np.pi*self.r_ft-self.r_ft*theta
            Hydror = area/per

            Qcap = 1.49/self.rough*area*Hydror**(2/3)*self.S**(0.5)
            # if Qcap < (self.Q[storm] + 0.1) and Qcap > (self.Q[storm] - 0.1):
            if Qcap > (self.Q[storm] - 0.1):
                self.depth[storm] = round(i,2)
                break
            if i == float(self.D_ft)-0.01:

                if self.size == False:
                    print (f'Pipe {self.id} is undersized, suggest upsizing pipe from a {self.D}" to a {standardsizes[standardsizes.index(self.D)+1]}" pipe.')
                    self.depth[storm] = 9999999999
                else:
                    print (f'Pipe {self.id} was undersized, changed pipe from a {self.D}" to {standardsizes[standardsizes.index(self.D)+1]}".')
                    self.D = standardsizes[standardsizes.index(self.D)+1]
                    self.Calc_Qcap()
                    self.Calc_depth(storm)

        self.v[storm] = round(self.Q[storm]/area,2)


class ConduitResults():
    def __init__(self, pdic, events):
        self.pdic = pdic
        self.events = events

    def design_Output(self):
        #pdic = dictionary of all conduit objects
        #events = storm events in list
        #ofile = outfile

        with open(r'L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\hydro\results\pipe_design_report.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            header = ['Name', 'Original Diameter', 'Diameter', 'Slope', 'Q_cap_orig', 'Q_cap']
            for event in self.events: header.append(f'Q {event}')
            for event in self.events: header.append(f'V {event}')
            for event in self.events: header.append(f'depth {event}')
            writer.writerow(header)

            for pipe in sorted(self.pdic):
                pname = pipe
                pipe = self.pdic[pipe]
                row = [pname, pipe.D_orig, pipe.D, pipe.S, pipe.Q_cap_orig, pipe.Q_cap]
                for event in self.events: row.append(pipe.Q[event])
                for event in self.events: row.append(pipe.v[event])
                for event in self.events: row.append(pipe.depth[event])
                writer.writerow(row)


    def finalTable_Output(self):
        with open(r'L:\LAProj\10067.006 - Building D-2 - Victory Logistics\Civil\Hydrology\calculations\hydro\results\pipe_report.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            header = ['Name', 'Diameter', 'Slope', 'Q_cap']
            for event in self.events: header.append(f'Q {event}')
            for event in self.events: header.append(f'V {event}')
            for event in self.events: header.append(f'depth {event}')
            writer.writerow(header)

            for pipe in sorted(self.pdic):
                pname = pipe
                pipe = self.pdic[pipe]
                row = [pname, pipe.D, pipe.S, pipe.Q_cap]
                for event in self.events: row.append(pipe.Q[event])
                for event in self.events: row.append(pipe.v[event])
                for event in self.events: row.append(pipe.depth[event])
                writer.writerow(row)


catch = pd.read_csv(r"C:\Users\imahaffey\Documents\catchments.csv")
conduit = pd.read_csv(r"C:\Users\imahaffey\Documents\conduit.csv")

catchs = catch.set_index('Catchment_ID').T.to_dict()
conduits = conduit.set_index('Name').T.to_dict()



# Add downstream flows to the conduits from catchments
for catch in catchs:
    catchname = catch
    catch = catchs[catch]
    try:
        conduits[catch['DS']]['Q_25'] = round(catch['Q_25'],2)
        conduits[catch['DS']]['Q_5'] = round(catch['Q_5'],2)
        print (f'\n\n-- Added {catch["Q_25"]} cfs from {catchname} to {catch["DS"]} --\n')
    except KeyError:
         continue


    c_pipe = conduits[catch['DS']] #current ds pipe dictionary
    c_name = catch['DS']
    try:
        ds = c_pipe['DS'] #name of ds pipe
        while ds != 'nan':
            if 'Q_25' not in conduits[ds]:
                conduits[ds]['Q_25'] = 0
                conduits[ds]['Q_5'] = 0
            print (f'- Added flows from {c_name} to {ds} -')
            print (f'Old flow = {conduits[ds]["Q_25"]} cfs')

            conduits[ds]['Q_25'] += round(catch['Q_25'],2)
            conduits[ds]['Q_5'] += round(catch['Q_5'],2)

            print (f'New flow = {conduits[ds]["Q_25"]} cfs')
            c_pipe = conduits[ds]
            c_name = ds
            ds = c_pipe['DS']

    except KeyError:
        continue


pipes = {}
for i in conduits:
    condname = i
    dic = conduits[i]
    print (f'Creating pipe object for {condname}')

    pipes[i] = Conduit(i, mat=dic['Material'], shape=dic['HydShape'], D=dic['Diameter'], L=dic['Length'],
                        x1=dic['x1'], x2=dic['x2'], y1=dic['y1'], y2=dic['y2'], S=dic['Slope'],
                        B=dic['Span'], rough=.013, Q={5: dic['Q_5'], 25: dic['Q_25']})
    pipes[i].Calc_Qcap()
    pipes[i].size = True
    pipes[i].Calc_depth(5)
    pipes[i].Calc_depth(25)


results = ConduitResults(pipes, [5,25])
results.design_Output()
results.finalTable_Output()
