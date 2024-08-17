from ENIIGMA.GA import optimize
filename = 'Global_Fit_Optical_Depth_new.png'
list_sp = ['H2O_40K', 'H2O_NH3_CO2_CH4_10_1_1_1_72K_b', 'd_NH3_CH3OH_50_10K_I10m_Baselined', 'CO_NH3_10K', 'H2O_CH4_10_0.6_a_V3', 'CO_CH3OH_10K', 'HNCO_NH3']
optimize.ENIIGMA(filename, 2.84, 4., list_sp, group_comb=3, skip=False, pathlib = None)
