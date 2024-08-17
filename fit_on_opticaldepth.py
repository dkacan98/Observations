import warnings
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import operator
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D, SpectralRegion
from specutils.fitting import fit_continuum
df=pd.read_csv('optical_depth.csv',sep=',',header=None,skiprows=[0],dtype=np.float64)
df.dropna(inplace=True)         #None values in csv file are dropped
df.rename(columns={0:'x', 1:'y'}, inplace=True)
x=df['x'].tolist()
y=df['y'].tolist()
spectrum = Spectrum1D(flux=y*u.MJy/u.sr, spectral_axis=x*u.micron)
region = [(4.73*u.micron,4.757*u.micron),(4.83*u.micron,4.838*u.micron),(4.855*u.micron,4.86*u.micron)]
with warnings.catch_warnings():  # Ignore warnings
	warnings.simplefilter('ignore')

	fitted_continuum=fit_continuum(spectrum,window=region,model=models.Polynomial1D(3))
A=[]
for x1 in x:

	if 4.73 <= x1 <= 4.86:
		A.append(x1)

y_continuum_fitted = fitted_continuum(A*u.micron)
y_continuum_fitted = y_continuum_fitted * u.sr / u.MJy
#for i in y_continuum_fitted:
f, ax = plt.subplots()  
#ax.xaxis.set_ticks(np.arange(4.7, 4.85, 0.02))
ax.plot(x, y)
ax.plot(A*u.micron, y_continuum_fitted, linestyle='dashed', color='red') 
ax.invert_yaxis()
ax.set_xlim(4.7,4.9)
ax.set_ylim(0.2,-0.1)
ax.set_xlabel('Wavelength(micron)')
ax.set_ylabel('Optical Depth')
ax.legend(["Observation Data","Local Continuum"],loc="upper right")
ax.set_title("4.8 local fit")  
ax.grid(True)  
#plt.savefig('best_fit.png')
plt.show()
plt.close()


#optical depth figure
xi=x.index(A[0])
xf=x.index(A[-1])
index_list=list(range(xi,xf+1,1))
flux = [y[i] for i in index_list]

new_flux = [x+0.05 for x in flux]
new_flux = np.array(new_flux)
new_y_continuum_fitted = [x+0.05 for x in y_continuum_fitted]
new_y_continuum_fitted = np.array(new_y_continuum_fitted)

#opdepth = -np.log(new_flux/new_y_continuum_fitted)
opdepth = new_y_continuum_fitted-new_flux
f, ax = plt.subplots()  
ax.xaxis.set_ticks(np.arange(4.7,4.9,0.02))
ax.plot(A*u.micron, opdepth)
#ax.invert_yaxis()
ax.set_xlabel('Wavelength(micron)')
ax.set_ylabel('Optical Depth')
ax.set_title('Optical depth plot')
ax.grid(True)
plt.show()

