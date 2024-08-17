import pandas as pd
df=pd.read_csv('l1527_prism_nrs1_s4d.tsv',sep='\t',header=None)
df.rename(columns={0:'x', 1:'y'}, inplace=True)
x=df['x']
print(x)
