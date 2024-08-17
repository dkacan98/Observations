import warnings
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import operator
import csv
import os
import pickle
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D, SpectralRegion
from specutils.fitting import fit_continuum
df=pd.read_csv('l1527_prism_nrs1_s4d.tsv',sep='\t',header=None,dtype=np.float64)
df.dropna(inplace=True)             #None values in csv file are dropped
df.rename(columns={0:'x', 1:'y'}, inplace=True)
x=df['x'].tolist()
y=df['y'].tolist()
multiply_list = [1000]*len(x)
res = list(map(operator.mul,x,multiply_list))    #to convert mm values in x to micrometer
spectrum = Spectrum1D(flux=y*u.MJy/u.sr, spectral_axis=res*u.micron)
region = [(1.2*u.micron,2.49*u.micron),(2.65*u.micron,2.69*u.micron),(2.715*u.micron,2.72*u.micron),(3.8*u.micron,4.0*u.micron),(4.06*u.micron,4.07*u.micron),(4.96*u.micron,5.3*u.micron)]
#region = [(1.2*u.micron,2.72*u.micron),(3.9*u.micron,4.0*u.micron),(5.0*u.micron,5.3*u.micron)]
with warnings.catch_warnings():  # Ignore warnings
	warnings.simplefilter('ignore')
	fitted_continuum=fit_continuum(spectrum,window=region,model=models.Polynomial1D(4))
"""
fname = "fit_function.pickle"
os.mkfifo(fname)
with open(fname,'w') as f:
	pickle.dump(fitted_continuum, f)
"""
A=[]
for x1 in res:

	if 1.2 <= x1 <= 5.3:
		A.append(x1)

y_continuum_fitted = fitted_continuum(A*u.micron)
print(len(y_continuum_fitted))
f, ax = plt.subplots()  
ax.xaxis.set_ticks(np.arange(0, 5.3, 0.1))
ax.plot(res, y)
ax.plot(A*u.micron, y_continuum_fitted,linestyle='dashed')  
ax.set_xlim(0,5.3)
ax.set_ylim(0,200)
ax.set_xlabel('Wavelength (micron)')
ax.set_ylabel('Flux (MJy/sr)')
ax.set_title("Global Continuum Fit")  
ax.grid(True)  
ax.fill_betweenx(y, 1.2, 2.49, color='lightgrey', alpha=.5) 
ax.fill_betweenx(y, 2.65, 2.69, color='lightgrey', alpha=.5) 
ax.fill_betweenx(y, 2.715, 2.72, color='lightgrey', alpha=.5) 
ax.fill_betweenx(y, 3.8, 4.0, color='lightgrey', alpha=.5)
ax.fill_betweenx(y, 4.06, 4.07, color='lightgrey', alpha=.5)
ax.fill_betweenx(y, 4.96, 5.3, color='lightgrey', alpha=.5)
#plt.savefig('best_fit.png')
plt.show()
plt.close()

#optical depth figure
xi=res.index(A[0])
xf=res.index(A[-1])
index_list=list(range(xi,xf+1,1))
flux = [y[i] for i in index_list]
flux=flux*u.MJy/u.sr
opdepth = -np.log(flux/y_continuum_fitted)
print(type(A),type(opdepth))

f, ax = plt.subplots()  
ax.xaxis.set_ticks(np.arange(0,5.3,0.1))
ax.plot(A*u.micron, opdepth)
ax.invert_yaxis()
ax.set_xlabel('Wavelength(micron)')
ax.set_ylabel('Optical Depth')
ax.set_title('Optical depth plot')
ax.grid(True)
ax.yaxis.set_ticks(np.arange(-1,5,0.2))
plt.show()
#with open('optical_depth.csv', 'w') as f:
#    writer = csv.writer(f)
#    writer.writerow(zip(A, opdepth))
f = pd.DataFrame(data={"col1": A, "col2": opdepth})
f.to_csv("./optical_depth.csv", sep=',',index=False)
