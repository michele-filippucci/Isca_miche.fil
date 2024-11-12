import os

import numpy as np

from isca import GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES =16

# Point to code as defined by $GFDL_BASE
cb = GreyCodeBase.from_directory(GFDL_BASE)

base_dir = os.path.dirname(os.path.realpath(__file__))


cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('exeter_project_mf_aa_KS_30_72,90_ws', codebase=cb)
exp.inputfiles = [os.path.join(base_dir,'input/ocean_qflux_anom_30_72,90.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
diag.add_file('atmos_daily', 1, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'temp', time_avg=False, files=['atmos_daily'])
diag.add_field('dynamics', 'ucomp', time_avg=False, files=['atmos_daily'])
diag.add_field('dynamics', 'vcomp', time_avg=False, files=['atmos_daily'])
diag.add_field('dynamics', 'height', time_avg=False, files=['atmos_daily'])
diag.add_field('atmosphere', 'temp_2m', time_avg=False, files=['atmos_daily'])

#Tell model which diagnostics to write
diag.add_field('dynamics', 'zsurf',time_avg=True, files=['atmos_monthly']) #need at least ps, pk, bk and zs$
diag.add_field('atmosphere', 'precipitation', time_avg=True, files=['atmos_monthly'])
diag.add_field('mixed_layer', 't_surf', time_avg=True, files=['atmos_monthly'])
diag.add_field('mixed_layer', 'flux_t', time_avg=True, files=['atmos_monthly']) #sensible heat flux
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True, files=['atmos_monthly']) #latent heat flux

exp.diag_table = diag # register diag table

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml': {
        'days'   : 30, # each run lasts one year, and then multiple runs are strung together below (loop on e.g. Line 310)
        'hours'  : 0,   # a different output file is produced for each run (in this case, each year). Data 
        'minutes': 0,   # output at the frequency specified in the diag table 
        'seconds': 0,
        'dt_atmos':480, # 900s timestep for dynamical core
        'current_date' : [1,1,1,0,0,0],
        'calendar' : 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'two_stream_gray': True, # Use RRTM radiation, not grey
        'do_rrtm_radiation': False, 
        'convection_scheme': 'SIMPLE_BETTS_MILLER', # Use the simple Betts Miller convection scheme 
        'do_damping': True, # turns on 'damping_driver', which manages the sponge at the top 
                            # -- this is not in original O'Gorman Schneider
        'turb':True, # turns on boundary layer diffusion, managed by 'vert_turb_driver' and 'gcm_vert_diff'
        'mixed_layer_bc':True, # turns on mixed layer bc 
        'do_virtual' :True, # determines whether virtual temperature is used for diffusion module 
        'roughness_mom':5.e-3, # DEFAULT: 0.05
        'roughness_heat':1.e-5, # DEFAULT: 0.05
        'roughness_moist':1.e-5, # DEFAULT: 0.05
                                 # Each of these have been set to their value in O'Gorman and Schneider 
        'land_roughness_prefactor':1.0,
        'do_simple':False,
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False, # Turn off Mellor-Yamada scheme 
        'do_diffusivity': True,  # Use non-local k-profile scheme for BL turbulence, as in O'Gorman and Schneider
    },
    
    'diffusivity_nml': {
        'do_entrain':False, 
        'free_atm_diff':False, # Turns off diffusion in free atmosphere 
        'do_simple': False, # Use virtual temperature for diffusion / when computing stability 
        'parcel_buoy': 0.0, # This change will lead to slightly reduced BL depth in unstable conditions. 
        'frac_inner': 0.1,  # Determines depth of surface layer in BL 
        'fixed_depth': False, # calculate pbl depth using stability criterion 
    },

    'surface_flux_nml': {
        'use_virtual_temp': True, # use virtual temp to compute stability of surface layer 
        'do_simple': False, # don't simplify computation of surface specific humidity at saturation 
        'old_dtaudv': False, # don't simplify derivative of surface wind stress (would set dwind_stress/du = dwind_stress/dv)  
        'gust_const':1.0, # Set constant gustiness factor at surface for computation of surface fluxes -- as in O'Gorman and Schneider 
        'land_humidity_prefactor' : 1.0, 
        'land_evap_prefactor': 1.0,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True # Use idealized_moist_phys 
    },


    'mixed_layer_nml': { 
        'depth': 2.5, # Use 2.5m mixed layer depth as in the Frierson case. 
                      # -- Different from O'Gorman and Schneider, who use 1m mixed layer 
        'albedo_value': 0.38, # As in O'Gorman and Schneider 
        'prescribe_initial_dist':True, # prescribes initial t_s with t_s = tconst - delta_T*(3*sin^2(lat)-1)/3
                                       # where tconst and delta_T are namelist options
        'tconst' : 285., # DEFAULT: 305 
                         # tconst in expression above 
        'delta_T': 40., # DEFAULT: 40. 
                        # delta_T in expression above 
        'evaporation':True, # allow surface evaporation and associated latent heat flux 
        'do_qflux': True,  # do prescribed ocean heat transport out of tropics 
                           # -- not in O'Gorman and Schneider, configured with qflux_nml
        'land_option':'none', 
        'load_qflux':True, #load triangle heating function
        'time_varying_qflux' : True, #triangle heating function will not be time-varying
        'qflux_file_name':'ocean_qflux_anom_30_72,90',
        
    },
    
    'qflux_nml':{
        'qflux_amp':50., # amplitude of q-flux heat transport out of tropics 
    }, 

    'qe_moist_convection_nml': { # Namelist for simple betts--miller convection 
        'rhbm':0.7, # relative humidity of reference profile used by convection scheme -- set following O'Gorman Schneider 
        'tau_bm':7200., # timescale for convective relaxation -- set following O'Gorman Schneider 
        'Tmin':120., # minimum temperature at LCL 
        'Tmax':360.,  # maximum temperature at LCL 
        'val_inc': 0.01, # increment in value for LCL lookup table
    },
    
    'lscale_cond_nml': {
        'do_simple':False, # if true then latent heat computation ignores latent heat of sublimation 
        'do_evap':False # turns off re-evaporation of falling precipitation -- follows O'Gorman Schneider 
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True, # This means model is using simplified expression: 
                          # esat = esat_0 * exp[ -(Lv/Rv) * (1/T - 1/T0) ] 
                          # instead of proper lookup table for esat 
                          # Follows O'Gorman Schnieder 

                          
    },    
    
    
    
    'damping_driver_nml': {
        'do_rayleigh': True, # turns on sponge (Rayleigh drag) at top of atmosphere, not included in O'Gorman and Schneider 
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  2600.,
        'do_conserve_energy': True,         
    },
    



    'two_stream_gray_rad_nml': { # each option is configured as in O'Gorman and Schneider 
        'rad_scheme': 'frierson',            #Select radiation scheme to use
        'do_seasonal': True,                #do_seasonal=false uses the p2 insolation profile
        'atm_abs': 0.22,                      # default: 0.0  
        'solar_exponent':2,
        'ir_tau_eq':7.2, 
        'ir_tau_pole':1.8, 
        'del_sol':1.2, 
        'solar_constant':1360, 
        'linear_tau':0.2, 
        'odp':1.0 # this is a multiplicative pre-factor for the longwave optical depth
    },


    'spectral_dynamics_nml': {
        'damping_order': 4, # Yields lap^8 damping 
        'water_correction_limit': 200.e2, # DEFAULT: 0. 
                                          # adds upper limit to water correction (which corrects small non-conservation
                                          # of water in the dynamical core)
        'reference_sea_level_press':1.0e5, # DEFAULT: 101325. 
                                           # used to construct hybrid coord and in implicit timestepping 
                                           # note actual mean sea level pressure is set in constants_mod 
        'num_levels':30, # Number of levels corresponding to set below 
        'valid_range_t':[100.,800.], # just set this to be a wide temperature range 
        'initial_sphum':[2.e-6], # DEFAULT: 0.0 
                                 # start the model with some water in atmosphere 
        'use_virtual_temperature':True,
        'vert_coord_option':'uneven_sigma', 
        'robert_coeff':0.03, # DEFAULT: 0.04, used in Robert filter for timestepping 
        # set to T85 resolution (default)
        'lon_max': 256, # DEFAULT: 128
        'lat_max': 128, # DEFAULT: 64
        'num_fourier': 85, # DEFAULT: 42
        'num_spherical':86, # DEFAULT: 43
        ### Slightly different from Tapio Schneider Github, but shouldn't mattter
        'surf_res': 0.05, 
        'exponent': 3., 
        'scale_heights': 5 
        
    }, 
   
    
    
    # FMS Framework configuration -- haven't modified these from defaults present in all experiment scripts
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },


    
    
})

# Lets do a run!
if __name__=="__main__":

    run_exp = exp.derive('exeter_project_mf_aa_KS_30_72,90') # derive experiment to run (this one is unchanged$
    #run_exp.run(31, use_restart='$GFDL_DATA/exeter_project_mf_aa_KS_30_72,90/restarts/res0030.tar.gz', num_cores=NCORES, overwrite_data=False)
    run_exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)
    for i in range(2,121): # run for 10 years
        run_exp.run(i, num_cores=NCORES, overwrite_data=False)

        
        


