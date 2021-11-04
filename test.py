import pandas as pd
import numpy as np

df = pd.DataFrame(np.arange(1,100000))
col = ['col']
df.columns = col

df['col'] = df['col'] * 5

print (df)



