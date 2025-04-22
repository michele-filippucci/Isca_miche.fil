import os
import subprocess
from isca.util import interpolate_output

data_directory = '/home/links/mf671/DATA/isca_data/'

for exp in ["exeter_project_mf_baseline_KS_qflux25"]:
    print(exp)
    for i in range(106,133):
        print("run"+'%04d'%i)
        infile = data_directory + exp + '/run' + '%04d'%i + '/atmos_daily.nc'
        outfile = data_directory + exp + '/run' + '%04d'%i + '/atmos_plev_daily.nc'
        
        # Perform interpolation
        interpolate_output(infile, outfile, p_levs='EVEN', var_names=['height','temp','ucomp', 'vcomp'])
        
