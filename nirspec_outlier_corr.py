import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAperture

"""
Perform a manual outlier correction on the nirspec data
"""

# load in nirspec data
file = '../../data/nirspec/jw01798-o004_t003_nirspec_g395m-f290lp_s3d.fits'
data, header = fits.getdata(file, extname='SCI', header=True)
# convert data cube from MJy/sr to Jy
sr_per_pixel = header["PIXAR_SR"]
data = data * 1e6 * sr_per_pixel
# get primary hdu
primary_hdu = fits.open(file)[0]

shape = np.shape(data)
kernel_radius = 3 # radial of circular kernel used for median smoothing
outlier_thresh = 500 # threshold to determine if point is an outlier

data_corr = np.zeros(shape)
for chan in range(shape[0]):
    # print an update every 100 channels
    if chan % 100 == 0:
        print("Correcting Channel ", chan)
    
    img = data[chan]
    # smooth the img by using a circular median filter
    smooth_img = np.zeros(np.shape(img))
    for i in range(shape[1]):
        for j in range(shape[2]):
            # if the i, j pixel is nan, keep it a nan
            if np.isnan(img[i, j]):
                smooth_img[i, j] = np.nan
            else:
                # create a circular aperture to use as kernel
                ap = CircularAperture((i, j), r=kernel_radius)
                mask = ap.to_mask(method='exact')
                # get the data in the aperture as a 1d array
                ap_data = mask.get_values(img)
                # set the value of i, j pixel to the median of pixels in aperture
                smooth_img[i, j] = np.nanmedian(ap_data)
    
    # take the residual between the original and smoothed image
    residual = img - smooth_img
    # find the MAD (median absolute deviation) of the residual
    med = np.nanmedian(residual)
    mad = np.nanmedian(np.abs(residual - med))
    # define outliers as points that where the |residual| is more than the
    # outlier threshold * MAD
    outlier_inds = np.where(np.abs(residual) > outlier_thresh * mad)

    # create a corrected image
    img_corr = np.copy(img)
    # replace outlier points in the original image with the corresponding points
    # from the smoothed image
    img_corr[outlier_inds] = smooth_img[outlier_inds]
    data_corr[chan] = img_corr

# save outlier corrected data as a fits file
sci_hdu = fits.ImageHDU()
sci_hdu.data = data_corr
# update the header to have the correct units
header['BUNIT'] = 'Jy'
sci_hdu.header = header
hdul = fits.HDUList([primary_hdu, sci_hdu])
hdul.writeto('../p_data/nirspec/nirspec_corr.fits')
