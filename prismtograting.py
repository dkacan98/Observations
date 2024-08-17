import sys
sys.path.insert(0, '/home/deno/Masaüstü/Research/Observations/PRISM/disk-integrated')
from global_continuum_fit import fitted_continuum
import warnings
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import operator
import csv
import os
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D, SpectralRegion
from specutils.fitting import fit_continuum
df=pd.read_csv('l1527_g395h_nrs2_s4d.tsv',sep='\t',header=None,dtype=np.float64)
df.dropna(inplace=True)             #None values in csv file are dropped
df.rename(columns={0:'x', 1:'y'}, inplace=True)
x=df['x'].tolist()
y=df['y'].tolist()
multiply_list = [1000]*len(x)
res = list(map(operator.mul,x,multiply_list))    #to convert mm values in x to micrometer
spectrum = Spectrum1D(flux=y*u.MJy/u.sr, spectral_axis=res*u.micron)

A=[]
for x1 in res:

	if 4.1 <= x1 <= 5.2:
		A.append(x1)

y_continuum_fitted = fitted_continuum(A*u.micron)
f, ax = plt.subplots()  
ax.xaxis.set_ticks(np.arange(1, 5.3, 0.1))
ax.plot(res, y)
ax.plot(A*u.micron, y_continuum_fitted,linestyle='dashed')  
ax.set_xlim(2,6)
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
