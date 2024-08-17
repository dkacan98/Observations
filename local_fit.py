import warnings
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import operator
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D, SpectralRegion
from specutils.fitting import fit_continuum
df=pd.read_csv('l1527_prism_nrs1_s4d.tsv',sep='\t',header=None,dtype=np.float64)
df.dropna(inplace=True)         #None values in csv file are dropped
df.rename(columns={0:'x', 1:'y'}, inplace=True)
x=df['x'].tolist()
y=df['y'].tolist()
multiply_list = [1000]*len(x)
res = list(map(operator.mul,x,multiply_list))    #to convert mm values in x to micrometer
spectrum = Spectrum1D(flux=y*u.MJy/u.sr, spectral_axis=res*u.micron)
region = [(4.33*u.micron,4.36*u.micron),(4.42*u.micron,4.48*u.micron)]
with warnings.catch_warnings():  # Ignore warnings
	warnings.simplefilter('ignore')

	fitted_continuum=fit_continuum(spectrum,window=region,model=models.Polynomial1D(4))
A=[]
for x1 in res:

	if 4.33 <= x1 <= 4.48:
		A.append(x1)

y_continuum_fitted = fitted_continuum(A*u.micron)
print(len(y_continuum_fitted))
f, ax = plt.subplots()  
ax.xaxis.set_ticks(np.arange(4.30, 4.55, 0.02))
ax.plot(res, y)
ax.plot(A*u.micron, y_continuum_fitted, linestyle='dashed', color='red') 

ax.fill_betweenx(y, 4.33, 4.36, color='lightgrey', alpha=.5) 
ax.fill_betweenx(y, 4.42, 4.48, color='lightgrey', alpha=.5) 
#ax.fill_betweenx(y, 3.27, 3.31, color='lightgrey', alpha=.5) 
#ax.fill_betweenx(y, 3.7, 4.0, color='lightgrey', alpha=.5)
ax.set_xlim(4.30,4.55)
ax.set_ylim(30,60)
ax.set_xlabel('Wavelength(micron)')
ax.set_ylabel('Flux (MJy/sr)')
ax.legend(["Observation Data","Local Continuum"],loc="upper right")
ax.set_title("4.9 local fit")  
ax.grid(True)  
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
f, ax = plt.subplots()  
ax.xaxis.set_ticks(np.arange(4.30,4.55,0.02))
ax.plot(A*u.micron, opdepth)
ax.invert_yaxis()
ax.set_xlabel('Wavelength(micron)')
ax.set_ylabel('Optical Depth')
ax.set_title('Optical depth plot')
ax.grid(True)
plt.show()

