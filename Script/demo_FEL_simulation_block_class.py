'''
SCRIPT OF USAGE OF THE FEL_SIMULATION_CLASS.

29/06/17: HMCC
'''

import ocelot
from ocelot.S2E_STFC import FEL_simulation_block

#################################################
#### Setting up the data dictionary

data = {'gen_file':'/home/qfi29231/SHORT_NOM_TD_FAST.in',
       'file_pout':'/scratch2b/qfi29231/test_fast_ASTRA_mat_aft/',
       'stat_run':1,
       'gen_launch':'genesis2.0',
       'i_edist':0,
       'i_beam':0,
        'i_dpa':0,
        'i_rad':0,
        'i_dfl':0,
       'i_astra':1,
       'astra_file':'/home/qfi29231/ASTRA_files/CLARA/clara_V10.astra',
       'i_match':1,
        'tw_match':{'alpha_x':-0.1895,
        'alpha_y' : 0,
        'beta_x':3.7666,
        'beta_y':1.4489},
       'i_rewrite':0,
        'par_rew':{'ntail':0,
                   'nslice':1200},
        'idump':1,
       'i_scan':0,
       'parameter':'quadf',
       'init':12,
       'end':25,
       'n_scan':7}

#### Making an instance of the FEL_simulation_block class
f=FEL_simulation_block.FEL_simulation_block(data)

#### Read GENESIS input file (calling the method read_GEN_input_file within the class)
A_inp=f.read_GEN_input_file()

#### Calling the method to run the simulation and making the postprocessing (GEN_simul_preproc). Input: the Input object which is the outcome of the read_GEN_input_file

A_inp_aft=f.GEN_simul_preproc(A_inp)
