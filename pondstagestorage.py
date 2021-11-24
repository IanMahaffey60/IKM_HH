import pandas as pd

def acad_stagestorage_ratingcurve(ifile):
    '''
    Creates a dataframe from the stage storage results file output from ACAD.

    Edit: Have it output the nested lists due to current branch on - will update
    '''
    col= ['cont elev', 'cont area', 'depth', 'incr vol1', 'cum vol1', 'incr vol2', 'cum vol2']
    df = pd.read_csv(ifile, skiprows=10, sep="\s{2}", header = None, thousands=',', engine='python')
    df.columns = col
    df_sum = df.groupby(['cont elev']).agg('sum')
    df_max = df.groupby(['cont elev']).agg('max')
    df_result = df_sum.merge(df_max, on='cont elev')[['cont area_x', 'cum vol1_y']]
    df_result.rename(columns = {'cont elev':'elev', 'cont area_x':'footprint', 'cum vol1_y':'cum_vol'}, inplace = True)


    # Remove once hydro.py can utilize dataframe input stage storage instead of nested lists
    df_result.reset_index(inplace=True)
    df_result = df_result.values.tolist()

    return df_result
