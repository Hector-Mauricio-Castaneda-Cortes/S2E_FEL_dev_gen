#!/usr/bin/python
'''FEL taper (daughter class of the FEL simulation block.

HMCC: 01-12-17 Reimplementation of the FEL simulation block (creation of daughter class with corresponding methods)
'''
#################################################
### import of all modules that are required.
from __future__ import print_function
from ocelot import *
from ocelot.utils.xfel_utils import *
from ocelot.gui.genesis_plot import *
from ocelot.adaptors.genesis import *
from ocelot.common import globals
from ocelot.S2E_STFC.FEL_simulation_block import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc, rcParams
from copy import deepcopy, copy
import os, sys, errno


class FEL_taper(FEL_simulation_block):
    def __init__(self, *initial_data, **kwargs):
        super(FEL_taper, self).__init__(*initial_data, **kwargs)

    @classmethod
    def taper_stepwise(cls, n, n0, aw0, a1, a2):
        for i in xrange(0, len(n0)):
            if n <= n0[i]:
                return aw0 + (n - n0[i - 1]) * a1[i - 1] + (n - n0[i - 1]) ** 2 * a2[i - 1]
        return 1.0

    @classmethod
    def taper_defined(cls, inp, tap_start, lin_tap, quad_tap):
        #inp2 = deepcopy(inp)
        inp2 = inp
        n = getattr(inp, 'nsec')
        k_x = inp2.lat.sequence[2].Kx
        n0 = [0, tap_start, n]
        a1 = [0, lin_tap * k_x]
        a2 = [0, quad_tap]
        taper_func = lambda n: cls.taper_stepwise(n, n0, k_x, a1, a2)
        setattr(inp2, 'magin', 1)
        setattr(inp2, 'lat', taper(inp2.lat, taper_func))
        return inp2

    def scan_over_noise_real(self, path_in):
        ap_out = []
        for a_files in os.listdir(path_in+'scan_0'):
            a_f = \
            [path_in +'scan_0/'+ a_files + '/' + f_out for f_out in os.listdir(path_in + 'scan_0/'+ a_files) if f_out.endswith('.gout')][0]
            ap_out.append(read_out_file(a_f))
        return ap_out

    def scan_over_taper(self, path_in):
        a_tparam = getattr(self, 'par_taper')
        return [[read_out_file(path_in + a_files + '/' + a_file2 + '/scan_0/ip_seed_-1/run.0.gout') for a_file2 in os.listdir(path_in+a_files+'/')] for a_files in os.listdir(path_in) if a_files.startswith('n_start')]

    def GEN_simul_preproc(self, a_inp, **kwargs):
        inp_arr = []
        if self.i_taper == 1:
            for key,val in kwargs.items():
                exec(key + '=val')
            tap_lin = -tap_lin
        a_inp2 = a_inp
        g_out = []
        setattr(a_inp2, 'iallharm', 1)
        setattr(a_inp2, 'lat', super(FEL_taper, self).undulator_design(a_inp2)['Magnetic Lattice'])
        if (self.i_scan == 0):
            s_scan = range(1)
            num = self.stat_run
            run_ids = xrange(0, num)
            print('++++++++ No scan ++++++++++')
        elif (self.i_scan != 0):
            s_scan = np.linspace(self.init, self.end, self.n_scan)
            num = self.stat_run
            run_ids = xrange(0, num)
            print('++++ Number of noise realisations {0} ++++++'.format(num))

        # Set up some input parameters
        if getattr(a_inp2, 'itdp') == 0:
            setattr(a_inp2, 'type', 'steady')
        else:
            setattr(a_inp2, 'type', 'tdp')

        if (getattr(self, 'idump')) == 1:
            setattr(a_inp2, 'idump', 1)
            setattr(a_inp2, 'idmpfld', 1)
        else:
            setattr(a_inp2, 'idump', 0)
            setattr(a_inp2, 'idmpfld', 0)

            # Existent dist or beam file (if exists)
        if (getattr(self, 'i_edist') == 1) and (hasattr(self, 'file_edist')):
            a_inp2 = super(FEL_taper, self).GEN_existent_beam_dist_dpa_rad(a_inp2, 'edist')
        elif (getattr(self, 'i_beam') == 1) and (hasattr(self, 'file_beam')):
            a_inp2 = super(FEL_taper, self).GEN_existent_beam_dist_dpa_rad(a_inp2, 'beam')
        elif (getattr(self, 'i_rad') == 1) and (hasattr(self, 'file_rad')):
            a_inp2 = super(FEL_taper, self).GEN_existent_beam_dist_dpa_rad(a_inp2, 'rad')
        elif (getattr(self, 'i_dpa') == 1) and (hasattr(self, 'file_dpa')):
            a_inp2 = super(FEL_taper, self).GEN_existent_beam_dist_dpa_rad(a_inp2, 'dpa')
        elif (getattr(self, 'i_dfl') == 1) and (hasattr(self, 'file_dfl')):
            a_inp2 = super(FEL_taper, self).GEN_existent_beam_dist_dpa_rad(a_inp2, 'dfl')
        else:
            print('++++ No edist or beam or dpa or rad file available ++++++')

        # Read ASTRA file.
        if hasattr(self, 'i_astra') and getattr(self, 'i_astra') == 1 and hasattr(self, 'astra_file'):
            a_inp2 = super(FEL_taper, self).convert_ASTRA_edist(a_inp2)
            setattr(a_inp2, 'beam', None)
        elif (hasattr(self, 'i_astra') and getattr(self, 'i_astra') == 1) and not (hasattr(self, 'astra_file')):
            print('Path of  the ASTRA file not provided')
            return
        else:
            print('No need to read ASTRA file')

        # Rematch beam (if the flag has been set within the data dictionary)
        if (getattr(a_inp2, 'edist') != None) and hasattr(self, 'i_match') and (getattr(self, 'i_match') == 1):
            a_inp2 = super(FEL_taper, self).rematch_edist(a_inp2)
            setattr(a_inp2, 'nslice', 0)

        # Overwrite the simulation attributes if the user has new values for them defined in the input data structure
        if (hasattr(self, 'i_rewrite')) and (hasattr(self, 'par_rew')) and (getattr(self, 'i_rewrite') == 1):
            a_inp2 = super(FEL_taper, self).GEN_rewrite_par(a_inp2)
        else:
            pass

            # Running over noise realisations and/or scan parameters
        for n_par in s_scan:
            for run_id in run_ids:
                a_inp2.runid = run_id
                if (self.stat_run == 1):
                    setattr(a_inp2, 'ipseed', -1)
                else:
                    setattr(a_inp2, 'ipseed', np.random.randint(9999))
                if self.i_taper == 1:
                    a_inp2.run_dir = getattr(self, 'file_pout') + 'n_start_' + str(n_start) + '/tap_lin_' + str(
                        tap_lin) + \
                                     '/scan_' + str(n_par) + '/ip_seed_' + str(a_inp2.ipseed) + '/'
                else:
                    a_inp2.run_dir = getattr(self, 'file_pout') + 'scan_' + str(n_par) + '/ip_seed_' + str(
                        a_inp2.ipseed) + '/'
                try:
                    os.makedirs(a_inp2.run_dir)
                except OSError as exc:
                    if (exc.errno == errno.EEXIST) and os.path.isdir(
                            self.file_pout + 'scan' + str(n_par) + '/run_' + str(run_id)):
                        pass
                    else:
                        raise
                if not (not (self.i_taper == 1) or not hasattr(self, 'par_taper')):
                    print('+++++Taper optimisation, starting module {}, linear taper {}, quadratic {}'.format(n_start,
                                                                                                              tap_lin,
                                                                                                              tap_quad))
                    a_inp2 = self.taper_defined(a_inp2, n_start, tap_lin, tap_quad)
                elif self.i_taper == 1 and not hasattr(self, 'par_taper'):
                    print('+++++ Taper flag defined but no parameters set')
                    break
                    return
                setattr(a_inp2,'magin',1)
                launcher = get_genesis_launcher(self.gen_launch)
                print('+++++ Starting simulation of noise realisation {0}'.format(run_id))
                g = run_genesis(a_inp2, launcher)
                setattr(g, 'filePath', str(a_inp2.run_dir))
                inp_arr.append(g)
                if a_inp2.itdp == 1 and self.i_taper == 0:
                    plot_gen_out_all(handle=g, savefig=True, showfig=False,
                                     choice=(1, 1, 1, 1, 3.05, 1, 0, 0, 0, 0, 0, 1, 1), vartype_dfl=complex128,
                                     debug=1)
                a_inp2.latticefile = None
                a_inp2.outputfile = None
                a_inp2.edistfile = None
                a_inp2.beamfile = None
                a_inp2.fieldfile = None
                a_inp2.radfile = None
        print('++++ postprocessing (results saved in {0}results++++'.format(self.file_pout))
        if self.i_taper == 0:
            super(FEL_taper, self).post_processing(a_inp2, s_scan)
        return inp_arr

    def taper_scan(self):
        a_inp = super(FEL_taper, self).read_GEN_input_file()
        a_inp2= a_inp
        a_tparam = getattr(self, 'par_taper')
        a_energ = np.zeros((len(a_tparam['n_starting']), a_tparam['num']))
        a_cof = np.linspace(a_tparam['a_tmin'], a_tparam['a_tmax'], a_tparam['num'])

        setattr(self, 'stat_run', 1)
        setattr(self, 'i_taper', 0)
        f_old = deepcopy(getattr(self, 'file_pout'))
        setattr(self, 'file_pout', f_old + 'untapered/')
        inp_nom = self.GEN_simul_preproc(a_inp2)

        setattr(self, 'i_taper', 1)
        setattr(self, 'file_pout', f_old + 'taper_scan/')
        setattr(self,'i_dump',0)

        for i_n, n_st in enumerate(a_tparam['n_starting']):
            for i_cof, a_tap in enumerate(a_cof):
                a_out = self.GEN_simul_preproc(a_inp2, n_start=n_st, tap_lin=a_tap, tap_quad=0)
        setattr(self, 'file_pout', f_old)
        print('++++ postprocessing tapering(results saved in {0}results++++'.format(self.file_pout+'taper_scan/'))
        self.taper_postproc()

    def taper_postproc(self):
        print('++++++++ Pulse energy curve ++++++++++++')
        self.plot_pulse_en()
        print('++++++++ Relative difference ++++++++++')
        self.plot_relative()
    

    def plot_relative(self): 
        filename = 'relative_difference_tapering'
        a_tparam = getattr(self, 'par_taper')
        a_cof = np.linspace(a_tparam['a_tmin'], a_tparam['a_tmax'], a_tparam['num'])
        
        a_untap = self.scan_over_noise_real(getattr(self, 'file_pout') + 'untapered/')
        a_tap = self.scan_over_taper(getattr(self,'file_pout')+'taper_scan/')
        
        a_en = np.zeros((len(a_tap),len(a_tap[0])))
        for i in range(len(a_tap)):
            for j in range(len(a_tap[0])):
                a_en[i,j] = np.amax(a_tap[i][j].energy,axis=0)

        a_en_unt = np.amax(np.mean([getattr(a_out, 'energy') for a_out in a_untap], 0))
        a_rel_diff = ((a_en[:, :] - a_en_unt) / a_en_unt)
        
        colors = cm.rainbow(np.linspace(0, 1, num=a_rel_diff.shape[0]))
        fig = plt.figure(filename)
        plt.clf()
        ax = fig.add_subplot(1,1,1)
        art = []
        for n_s in range(a_rel_diff.shape[0]):
            ax.plot(a_cof[:], a_rel_diff[n_s,:], color=colors[n_s], label=r'n$_{segment} = $'
                                                                   + str(a_tparam['n_starting'][n_s]))
            lg = plt.legend(loc=9, ncol=3, bbox_to_anchor=(0.5, -0.16), prop={'size': 12})
            art.append(lg)
        
        ax.set_xlabel(r'Linear coefficient', fontsize=14)
        ax.set_ylabel(r'$\Delta$E$_{tapered-untapered}$/E$_{untapered}$', fontsize=14)
        ax.set_xlim([a_tparam['a_tmin'], a_tparam['a_tmax']])
        ax.set_xticks(a_cof[::2])
        ax.set_xticklabels([str(np.round(a_tic, 4)) for a_tic in a_cof[::2]])
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(9)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(9)
        ax.grid(True)
        fig.savefig(self.file_pout + 'taper_scan/' + filename+ '.' + 'png', format='png',dpi=1200,bbox_inches="tight",additional_artists=art)
        plt.draw()
        plt.close('all')
        

    def plot_pulse_en(self):
        from matplotlib.ticker import ScalarFormatter
        frmt = ScalarFormatter(useOffset=False)
        frmt.set_scientific(True)   
        a_tparam = getattr(self, 'par_taper')
        a_cof = np.linspace(a_tparam['a_tmin'], a_tparam['a_tmax'], a_tparam['num'])
        filename = 'pulse_energy_tapering'
        fig = plt.figure(filename)
        plt.clf()
        ax_rad_pow = fig.add_subplot(1,1,1)
        art = []
        a_tparam = getattr(self, 'par_taper')
        colors = cm.rainbow(np.linspace(0, 1, len(a_tparam['n_starting'])))

        a_untap = self.scan_over_noise_real(getattr(self, 'file_pout') + 'untapered/')
        a_tap = self.scan_over_taper((getattr(self, 'file_pout') + 'taper_scan/'))
        ax_rad_pow.plot(a_untap[0].z, 1e6 * np.mean([a_pout.energy for a_pout in a_untap], axis=0), '--',
                        color='forestgreen',
                        linewidth=2.4, label='Untapered')
        for jc in range(len(a_tparam['n_starting'])):
            for i in range(a_tparam['num']):
                if jc == 0:
                    ax_rad_pow.plot(a_tap[jc][i].z, 1e6 * a_tap[jc][i].energy,
                                    color=colors[jc],
                                    label=r'n$_{segment} = $' + str(a_tparam['n_starting'][jc])+ ', linear coefficient = ' + str(np.round(a_cof[i], 4)), linewidth=1.5)
                else:
                    ax_rad_pow.plot(a_tap[jc][i].z, 1e6 * a_tap[jc][i].energy, color=colors[jc], linewidth=1.5)
            lg = plt.legend(loc=9, ncol=4, bbox_to_anchor=(0.5, -0.16), prop={'size': 9})
        art.append(lg)

        ax_rad_pow.set_ylabel(r'E [uJ]', fontsize=14)
        ax_rad_pow.set_xlabel(r'z [m]', fontsize=14)

        plt.yticks(plt.yticks()[0][0:-1])

        ax_rad_pow.grid(True)  # , which='minor')
        ax_rad_pow.tick_params(axis='y', which='both', colors='g')
        ax_rad_pow.yaxis.label.set_color('g')
        ax_rad_pow.yaxis.get_offset_text().set_color(ax_rad_pow.yaxis.label.get_color()) 

       
        fig.savefig(self.file_pout + 'taper_scan/' + filename+ '.' + 'png', format='png',dpi=1200,additional_artists=art, bbox_inches="tight") 
        plt.draw()
        plt.close('all')
