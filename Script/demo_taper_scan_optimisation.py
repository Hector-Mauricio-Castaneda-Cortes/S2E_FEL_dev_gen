'''
SCRIPT OF USAGE OF THE FEL_taper.

08/12/17: HMCC Creation
'''


from ocelot.S2E_STFC.FEL_taper import *

data = {'gen_file': '/scratch2b/qfi29231/OPC/CLARA_2/tapered/LONG_NOM_TD_0.75m_new.in',
                'file_pout': '/scratch2b/qfi29231/test_tapering_module_flat-top/07-12/',
                'stat_run': 1,
                'gen_launch': 'genesis2.0',
                'i_edist': 0,
                'file_edist': '/home/qfi29231/e_dist_L-FEL_old.dist',
                'i_beam': 0,
                'file_beam': '/scratch2b/qfi29231/original_DAVE_file/TOM/beamfile2.beam',
                'i_dpa': 0,
                'file_dpa': 'file',
                'i_rad': 0,
                'file_rad': 'file',
                'i_dfl': 0,
                'file_dfl': 'dfl',
                'i_astra': 0,
                'astra_file': '/home/qfi29231/ASTRA_files/CLARA/clara_V10.astra',
                'i_match': 0,
                'tw_match': {'alpha_x': -0.1895,
                             'alpha_y': 0,
                             'beta_x': 3.7666,
                             'beta_y': 1.4489},
                'i_rewrite': 1,
                'par_rew': {'prad0': 12,
                            'shotnoise': 0,
                            'nsec':17,
                            'ncar':451,
                            'curpeak':400,
                            'curlen':-7.9666666666667e-5},
                'idump': 0,
                'i_scan': 0,
                'parameter': 'aw0',
                'init': 0.907,
                'end': 0.911,
                'n_scan': 20,
                'i_taper':1,
                'par_taper':{'n_starting': [7,8],
                             'num': 25,
                             'a_tmin': 0.001,
                             'a_tmax': 0.025},
        }
### Making an instance of the FEL_simulation_block class
f_obj = FEL_taper(data)

#### Taper scan method
f_obj.taper_scan()
