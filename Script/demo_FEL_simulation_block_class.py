'''
SCRIPT OF USAGE OF THE FEL_SIMULATION_CLASS.

29/06/17: HMCC
12/07/17: HMCC Add the extra flags corresponding to the existent file paths to be read
'''

import ocelot
from ocelot.S2E_STFC import FEL_simulation_block

#################################################
#### Setting up the data dictionary

data = {'gen_file':'/scratch2b/qfi29231/original_DAVE_file/TOM/NEW_LONG_NOM_TD.in',
       'file_pout':'/scratch2b/qfi29231/original_DAVE_file/TOM/test_edist/with_latt_no_beam/',
       'stat_run':1,
       'gen_launch':'genesis2.0',
       'i_edist':0,
        'file_edist':'/home/qfi29231/e_dist_L-FEL_old.dist',
       'i_beam':0,
        'file_beam':'/scratch2b/qfi29231/original_DAVE_file/TOM/beamfile2.beam',
        'i_dpa':0,
        'file_dpa':'file',
        'i_rad':0,
        'file_rad':'file',
        'i_dfl':0,
        'file_dfl':'dfl',
       'i_astra':0,
       'astra_file':'/home/qfi29231/ASTRA_files/CLARA/clara_V10.astra',
       'i_match':0,
        'tw_match':{'alpha_x':-0.1895,
        'alpha_y' : 0,
        'beta_x':3.7666,
        'beta_y':1.4489},
       'i_rewrite':0,
        'par_rew':{'aw0':0,
                   'awd':1200},
        'idump':0,
       'i_scan':0,
       'parameter':'aw0',
       'init':0.907,
       'end':0.911,
       'n_scan':20}

#### Making an instance of the FEL_simulation_block class
f=FEL_simulation_block.FEL_simulation_block(data)

#### Read GENESIS input file (calling the method read_GEN_input_file within the class)
a_inp=f.read_GEN_input_file()

#### Calling the method to run the simulation and making the postprocessing (GEN_simul_preproc). Input: the Input object which is the outcome of the read_GEN_input_file

a_inp_aft=f.GEN_simul_preproc(a_inp)
