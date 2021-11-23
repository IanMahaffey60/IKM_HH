class Results():
    '''
    Not really sure how I want to set this up yet. 
    '''
    def __init__(self, catchs, ponds):
        self.catchs = catchs #dictionary of catchments
        self.ponds = ponds #dictionary of ponds

    def print_assumptions(self):
        print ('\n------ Assumptions ------\n')

        for pond in self.ponds: print (f'Infiltration rate = {self.ponds[pond].infil} inches/hour')
        for catch in self.catchs: print (f'{catch} Rainfall Intensity (10-min Tc) = {self.catchs[catch].i} inches/hour')
        for catch in self.catchs: print (f'{catch} Time of Concentration = {self.catchs[catch].tc} minutes')
        for catch in self.catchs: print (f'{catch} Area = {self.catchs[catch].a} acres')

    def print_catchment_results(self, catch):
        print ('\n-------- Results ---------\n')

        print (f'Existing Peak Flow = {q_e.q} cfs')
        print (f'Existing Total Volume = {q_e.total_vol()} cu-ft')

    def print_pond_results(self):
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
