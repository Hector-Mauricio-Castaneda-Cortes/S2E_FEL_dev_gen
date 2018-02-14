import os,sys,shutil
import time
import numpy as np
from ocelot.S2E_STFC.s2e_ocelot import *
from ocelot.S2E_STFC.FEL_simulation_block import *
from ocelot.adaptors.genesis import *
from ocelot.gui.genesis_plot import *

def main(data0):
    print('+++++ S2E afterburner ++++++')
    f_obj = s2e_ocelot.s2e_afterburner(data0)
    print('+++++ S2E afterburner (Modulator) ++++++')
    tic0 = time.clock()
    f_mod=f_obj.modulator_amplifier('mod')
    print('Modulator simulation in {} seconds'.format(tic0))
    time.sleep(1)
    print('+++++ S2E afterburner (Amplifier) ++++++')
    f_ampl = f_obj.modulator_amplifier('amp')
    tic1=time.clock()
    print('Amplifier simulation in {} seconds'.format(tic1-tic0)) 
    time.sleep(1)
    print('+++++ S2E afterburner (after-burner) ++++++')
    f_out=f_obj.after_burner(f_ampl)
    tic2= time.clock()
    print('Amplifier simulation in {} seconds'.format(tic2-tic1)) 
    time.sleep(1)
    print('+++++ End of simulation ++++++')
    return

if __name__=="__main__":
    data = {'gen_file':'/home/qfi29231/riley_S2E/demo_mod.in',
       'file_pout':'/scratch2b/qfi29231/after_burner-class_dpa2edist2/',
       'stat_run':1,
       'gen_launch':'genesis2.0',
       'i_edist':0,
        'file_edist':'/home/qfi29231/ASTRA_files/50pC_matched_E2G_c.dist',
       'i_beam':0,
        'file_beam':'/scratch2b/qfi29231/original_DAVE_file/TOM/beamfile2.beam',
        'i_dpa':0,
        'file_dpa':'file',
        'i_rad':0,
        'file_rad':'file',
        'i_dfl':0,
        'file_dfl':'dfl',
       'i_astra':0,
       'astra_file':'/home/qfi29231/S2E_comparison_PUFFIN/test.in.128.4929.328',
       'i_match':0,
        'tw_match':{'alpha_x':-0,
        'alpha_y' : 0,
        'beta_x':0,
        'beta_y':0},
       'i_rewrite':0,
        'par_rew':{'aw0':0,
                   'awd':1200},
        'idump':0,
       'i_scan':0,
       'parameter':'quadf',
       'init':8,
       'end':15,
       'n_scan':7,
        'seedlambda':3e-7,
        'outlambda':1e-7,
        'nslice':3e3,
        'npart':2048,
        'gamma':490.2,
        'modules':10,
        'offsetradf':-16}
    main(data)

    
 
