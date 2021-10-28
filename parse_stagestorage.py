import pandas as pd

ifc = r"C:\Users\imahaffey\Documents\lot4_res.txt"

def create_ratingcurve(ifile):
    col= ['cont elev', 'cont area', 'depth', 'incr vol1', 'cum vol1', 'incr vol2', 'incr vol2']
    df = pd.read_csv(ifc, skiprows=10, sep="\s{2}", header = None, thousands=',', engine='python')
    df.columns = col
    return df
# print(data)

df = create_ratingcurve(ifc)
max_vol = df.loc[df['cont elev'] == 4071, 'cum vol1'].item()
print (round(max_vol/43560,2), ' acre-ft')
