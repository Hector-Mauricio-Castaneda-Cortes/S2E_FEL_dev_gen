import numpy as np
from ocelot.S2E_STFC.FEL_simulation_block import *
from ocelot.common.globals import *
from copy import copy, deepcopy
import matplotlib.pyplot as plt
import os, shutil
'''Ming Xie class.

Has the nsec, nwig, fl,drl as attributes (can be read from the GENESIS input file. The Arfken Current is set as a 
private attribute.  The user needs to provide the undulator parameter as well, the bin width, and the total number of 
particles. The outplot attributes corresponds to the directory where the results are to be stored. If the i_plot 
attribute is set to 0, a text file is generated in the outplot/mXie directory showing the maximum value of the pulse 
energy calculated via Ming Xie. Otherwise, a plot is generated showing the saturation power as a function of z, 
pulse energy as a function of z, rho parameter a function of s and the number of particles per slice.

HMCC: 22-01-18 Implementation of the Ming-Xie class

'''

class MingXie(object):
    __lsat = []
    __beam = None
    __rho = []
    __psat = []
    __glength = []
    __npart = []
    
    def _init__(self, **initial_value):
        for ke_y, value in initial_value.iteritems():
            setattr(self, ke_y, value)
        self.__constants = {'a1': 0.45, 'a2': 0.57, 'a3': 0.55, 'a4': 1.6, 'a5': 3, 'a6': 2,
                          'a7': 0.35, 'a8': 2.9, 'a9': 2.4, 'a10': 51, 'a11': 0.95, 'a12': 3,
                          'a13': 5.4, 'a14': 0.7, 'a15': 1.9, 'a16': 1140, 'a17': 2.2, 'a18': 2.9, 'a19': 3.2}
        self.__curra = 17.045e3
        self.__nwig = 0
        self.__nsec = 0
        self.__fl =0
        self.__drl = 0
        self.__npartmax = 0
        self.__wavelength = 0
        self.__binwidth = 0
        self.__lambdu = 0
        self.__K = 0
        self.__maxpulsee = 0
        self.__outplot = None
        self.__iplot = 0
   
    @property
    def curra(self):
        return 17.045e3

    @property
    def constants(self):
        return {'a1': 0.45, 'a2': 0.57, 'a3': 0.55, 'a4': 1.6, 'a5': 3, 'a6': 2,
                'a7': 0.35, 'a8': 2.9, 'a9': 2.4, 'a10': 51, 'a11': 0.95, 'a12': 3,
                'a13': 5.4, 'a14': 0.7, 'a15': 1.9, 'a16': 1140, 'a17': 2.2, 'a18': 2.9, 'a19': 3.2}
    @property
    def lsat(self):
        return self.__lsat 
    
    @property
    def rho(self):
        return self.__rho 

    @property
    def psat(self):
        return self.__psat

    @property
    def glength(self):
        return self.__glength

    @property
    def npart(self):
        return self.__npart
    
    @property
    def nsec(self):
        return self.__nsec

    @property
    def nwig(self):
        return self.__nwig

    @property
    def fl(self):
        return self.__fl

    @property
    def drl(self):
        return self.__drl 

    @property
    def npartmax(self):
        return self.__npartmax

    @property
    def K(self):
        return self.__K

    @property
    def lambdu(self):
        return self.__lambdu

    @property
    def binwidth(self):
        return self.__binwidth

    @property
    def wavelength(self):
        return self.__wavelength

    @property
    def maxpulsee(self):
        return self.__maxpulsee

    @property
    def outplot(self):
        return self.__outplot  
    
    @property
    def iplot(self):
        return self.__iplot
    
    @property
    def beam(self):
        return self.__beam
    
    @beam.setter
    def beam(self, beamf):
        self.__beam = beamf
    
    @nsec.setter
    def nsec(self, n_sec):
        self.__nsec= n_sec

    @nwig.setter
    def nwig(self, n_wig):
        self.__nwig = n_wig

    @fl.setter
    def fl(self, f_l):
        self.__fl = f_l

    @drl.setter
    def drl(self, dr_l):
        self.__drl = dr_l 
        
    @npartmax.setter
    def npartmax(self, n_pm):
        self.__npartmax = n_pm

    @binwidth.setter
    def binwidth(self, bin_w):
        self.__binwidth = bin_w

    @lambdu.setter
    def lambdu(self, lambda_u):
        self.__lambdu = lambda_u

    @wavelength.setter
    def wavelength(self, wavelength):
        self.__wavelength = wavelength 
    
    @lsat.setter
    def lsat(self, l_sat):
        self.__lsat = deepcopy(l_sat)
    
    @psat.setter
    def psat(self, p_sat):
        self.__psat = deepcopy(p_sat)
    
    @rho.setter
    def rho(self, r_ho):
        self.__rho = deepcopy(r_ho)

    @glength.setter
    def glength(self, g_l):
        self.__glength = deepcopy(g_l)

    @npart.setter
    def npart(self,n_p):
        self.__npart = deepcopy(n_p)

    @K.setter
    def K(self, K_0):
        self.__K = K_0

    @maxpulsee.setter
    def maxpulsee(self, pulse_energy):
        self.__maxpulsee = np.amax(pulse_energy)

    @outplot.setter
    def outplot(self, out_plot):
        self.__outplot= out_plot
    
    @iplot.setter
    def iplot(self,i_plot):
        self.__iplot = i_plot

    def slice_rho(self, sl):
        beam_p = getattr(self, 'beam')
        lmd_u = getattr(self, 'lambdu')
        cur_peak = beam_p.I[sl]
        gamma = beam_p.g0[sl]
        ex = beam_p.ex[sl]
        ey = beam_p.ey[sl]
        betax = beam_p.betax[sl]
        betay = beam_p.betay[sl]
        xrms = np.sqrt(betax*ex/gamma)
        yrms = np.sqrt(betay*ey/gamma)
        sigma = float((xrms+yrms)/2.0)
        return np.cbrt((cur_peak / self.curra) * np.square(lmd_u * self.K / (2.0 * np.pi * sigma)) * (np.power((1.0 / (2.0 * gamma)), 3)))
    
    def slice_npart(self,sl):
        beam_p = getattr(self, 'beam')
        curp = beam_p.I[sl]
        t_pulse = np.amax(beam_p.z)-np.amin(beam_p.z)
        #return np.floor(getattr(self,'npartmax')*getattr(self,'binwidth')/(t_pulse))
        #return np.floor(curp*(beam_p.z[sl] /speed_of_light)/q_e)
        return np.floor(curp*self.binwidth/(speed_of_light*q_e))

    def slice_eta(self, sl):
        beam_p = getattr(self, 'beam')
        lmb = getattr(self, 'wavelength')
        lmb_u = getattr(self, 'lambdu')
        a_cof = getattr(self, 'constants')
        ex = beam_p.ex[sl]
        ey = beam_p.ey[sl]
        delgam = beam_p.dg[sl]
        gamma = beam_p.g0[sl]
        betax = beam_p.betax[sl]
        betay = beam_p.betay[sl]
        xrms = np.sqrt(betax*ex/gamma)
        yrms = np.sqrt(betay*ey/gamma)
        sigma = float((xrms+yrms)/2.0)
        emit = float(((ex + ey) / 2.0))/gamma
        n_d = (lmb_u / (4.0 * np.pi * np.sqrt(3) * self.slice_rho(sl))) / ((4.0*np.pi *np.square(sigma))/lmb)
        n_e =float((lmb_u / (4.0 * np.pi * np.sqrt(3) * self.slice_rho(sl))) /(np.square(sigma)/emit))*float(4*np.pi*emit/lmb)
        ngam = float(4.0 * np.pi) * float(1.0 / (4.0 * np.pi * np.sqrt(3) * self.slice_rho(sl))) * delgam / gamma
        eta =  (a_cof['a1'] * (n_d**a_cof['a2'])) + (a_cof['a3'] * (n_e**a_cof['a4'])) + (a_cof['a5'] * (ngam**a_cof['a6'])) +\
            (a_cof['a7'] * (n_e**a_cof['a8']) * (ngam**a_cof['a9'])) +(a_cof['a10'] * (n_d**a_cof['a11']) * (ngam**a_cof['a12'])) + \
            (a_cof['a13'] * (n_d**a_cof['a14']) * (n_e**a_cof['a15'])) + \
            (a_cof['a16'] * (n_d**a_cof['a17']) * (n_e**a_cof['a18']) *(ngam**a_cof['a19']))
        return eta

    def slice_p_sat(self, sl):
        beam_p = getattr(self, 'beam')
        curpeak = beam_p.I[sl]
        gamma = beam_p.g0[sl]
        eta =  self.slice_eta(sl)
        return 1.6e6 * self.slice_rho(sl) * (np.square(1.0 / (1.0 + eta))) * curpeak * (gamma * m_e_MeV)

    def slice_p_noise(self, sl):
        beam_p = self.beam
        gamma = beam_p.g0[sl]
        lbd = self.wavelength
        return float(1.0e6*speed_of_light * m_e_MeV * q_e * gamma * np.square(self.slice_rho(sl)) / lbd)

    def slice_gain_length(self, sl):
        lbd_u = float(getattr(self, 'lambdu'))
        return float((1.0 + self.slice_eta(sl)) * lbd_u / (4.0 * np.pi * np.sqrt(3.0) * self.slice_rho(sl)))

    def slice_sat_length(self, sl):
        return self.slice_gain_length(sl) * np.log(9.0 * self.slice_p_sat(sl) / self.slice_p_noise(sl))

    def ming_xie_plotting(self, z, p_energ, p_noise, filename):
        import matplotlib
        from matplotlib.ticker import ScalarFormatter
        
        matplotlib.rcParams.update({'font.size':10})
        frmt = ScalarFormatter(useOffset = False)
        frmt.set_scientific(True)

        beam_z = getattr(self,'beam')
        setattr(self,'outplot',filename)
        fig , ax_mxie= plt.subplots(2,2,squeeze=False)
        ln1 = ax_mxie[0,0].plot(z,1e6*p_energ,color='navy')
        ax_mxie[0,0].scatter(z, 1e6 * p_energ, color='navy')
        ax_mxie[0,0].set_ylabel(r'E$_{pulse}$[$\mu$J]')
        ax_mxie[0,0].set_xlabel(r'z[m]')
        ax_mxie[0,0].grid(True)
        for ticks in ax_mxie[0,0].xaxis.get_major_ticks():
            ticks.label.set_fontsize(8)
        for ticks in ax_mxie[0,0].yaxis.get_major_ticks():
            ticks.label.set_fontsize(8)
        ax_mxie[0,0].yaxis.offsetText.set_fontsize(8)

        #ax_nois = ax_mxie[0, 0].twinx()
        #ax_nois.plot(z,1e6*np.sum(p_noise,axis=0)*getattr(self,'binwidth')/speed_of_light,'-.', color='lightcoral')
        #ln2 = ax_nois.scatter(z, 1e6 * np.sum(p_noise, axis=0) * getattr(self, 'binwidth') / speed_of_light,'-.', color='lightcoral', label=r'Noise Pulse Energy')
        #ax_nois.set_ylabel(r'E[$\mu$J]', fontsize=12)
        #ax_mxie[0, 0].legend(ln1+ln2,[l.get_label for l in (ln1+ln2)],loc = 0,prop={'size':10})
       # plt.suptitle(r'Properties of the Ming-Xie 1D model', color= 'red',fontsize=13)
        #ax_mxie[0, 0].legend(loc=2,prop={'size':10})
        ln01 = ax_mxie[0,1].plot(1e6*np.array(beam_z.z),getattr(self,'rho'),'--', color = 'forestgreen')
        ax_mxie[0,1].tick_params(r'$\rho$',color = 'forestgreen')
        ax_mxie[0,1].set_xlabel(r's[$\mu$m]')
        ax_mxie[0,1].set_ylabel(r'$\rho$',color='forestgreen')
        ax_gl = ax_mxie[0,1].twinx()
        ln02 = ax_gl.plot(1e6*np.array(beam_z.z),getattr(self,'glength'),'-.',color='orangered')
        ax_gl.set_ylabel(r'L$_g$[m]', color='orangered')
        ax_gl.tick_params(r'L$_g$[m]', color='orangered')
        #ax_gl.legend([ln01,ln02],[r'$\rho$',r'Gain length'], loc=1, prop = {'size':8})
        ax_mxie[0,1].tick_params(axis='y',which='both',colors= 'forestgreen')
        ax_gl.tick_params(axis='y',which='both',colors= 'orangered') 
        ax_mxie[0,1].yaxis.label.set_color('forestgreen')
        ax_gl.yaxis.label.set_color('orangered')
        for ticks in ax_mxie[0,1].xaxis.get_majorticklabels():
            ticks.set_fontsize(8) 
        for ticks in ax_mxie[0,1].yaxis.get_majorticklabels():
            ticks.set_fontsize(8)
        for ticks in ax_gl.yaxis.get_majorticklabels():
            ticks.set_fontsize(8) 
        ax_mxie[0,1].yaxis.offsetText.set_fontsize(8)
        ax_gl.yaxis.offsetText.set_fontsize(8)
        
        ax_mxie[1,0].plot(1e6*np.array(beam_z.z),p_noise,'--',color='brown')
        ax_mxie[1,0].scatter(1e6*np.array(beam_z.z),p_noise,color='brown')
        ax_mxie[1,0].tick_params(r'P$_{saturation}$[W]', color='brown')
        ax_mxie[1,0].set_xlabel(r's[$\mu$m]')
        ax_mxie[1,0].set_ylabel(r'P$_{saturation}$[W]', color='brown')
        for ticks in ax_mxie[1,0].xaxis.get_major_ticks():
            ticks.label.set_fontsize(8) 
        for ticks in ax_mxie[1,0].yaxis.get_major_ticks():
            ticks.label.set_fontsize(8) 
        ax_mxie[1,0].yaxis.offsetText.set_fontsize(8)
        ax_mxie[1,1].plot(1e6*np.array(beam_z.z),self.npart,'--',color='navy')
        ax_mxie[1,1].scatter(1e6*np.array(beam_z.z),self.npart,color='navy')
        ax_mxie[1,1].tick_params('N', color='navy')
        ax_mxie[1,1].set_xlabel(r's[$\mu$m]')
        ax_mxie[1,1].set_ylabel(r'N', color='navy')
        for ticks in ax_mxie[1,1].xaxis.get_major_ticks():
            ticks.label.set_fontsize(8) 
        for ticks in ax_mxie[1,1].yaxis.get_major_ticks():
            ticks.label.set_fontsize(8)
        plt.tight_layout()
        fig.savefig(filename,format='png',dpi=1200)

    def saturation_power_length(self):
        l_sat = []
        r_ho = []
        n_part = []
        p_sat = []
        p_noise = []
        n_slice = []
        lg = []
        ref = []
        p_min_lsat = []
        beam_z = getattr(self, 'beam')
        curr_mx = np.amax(beam_z.I)
        for i_sl in range(len(beam_z.z)):
            n_part.append(self.slice_npart(i_sl))
            if beam_z.I[i_sl] > 0.15 * curr_mx:
                p_sat.append(self. slice_p_sat(i_sl))
                p_noise.append(self.slice_p_noise(i_sl))
                n_slice.append(beam_z.I[i_sl] * q_e * self.binwidth / speed_of_light)
                lg.append(self.slice_gain_length(i_sl))
                l_sat.append(self.slice_sat_length(i_sl))
                ref.append(1.0) 
                r_ho.append(self.slice_rho(i_sl))
            else:
                p_sat.append(0)
                p_noise.append(0)
                n_slice.append(0)
                lg.append(0)
                l_sat.append(0)
                ref.append(0)
                r_ho.append(0)
        minlsat = np.amax(l_sat)
        for i_sl in range(len(beam_z.z)):
            if ref[i_sl] == 1 and l_sat[i_sl] <minlsat:
                minlsat = l_sat[i_sl]
        for i_sl in range(len(beam_z.z)):
            if ref[i_sl]==1:
                p_min_lsat.append(self.slice_p_noise(i_sl) * np.exp(minlsat / self.slice_gain_length(i_sl)) / 9.0)
            else:
                p_min_lsat.append(0.0)
        setattr(self,'lsat',np.array(l_sat))
        setattr(self, 'psat', np.array(p_sat))
        setattr(self, 'glength', np.array(lg))
        setattr(self,'rho',np.array(r_ho))
        setattr(self,'npart',np.array(n_part))
        p_min_lsat = np.array(p_min_lsat)
        p_noise = np.array(p_noise)
        ref = np.array(ref)
        return p_min_lsat, p_noise, ref

    def ming_xie_calculation(self):
        und_l = []
        l_und_l =[]
        p_und_l = []
        p_energy = []
        try:
            os.makedirs(getattr(self,'outplot'))
        except OSError as exc:
            if (exc.errno == errno.EEXIST):
                pass
                #shutil.rmtree(getattr(self,'outplot')) 
                #os.makedirs(getattr(self,'outplot'))
            else:
                raise
        
        p_min_s, p_noise, ref = self.saturation_power_length()

        beam_z = getattr(self,'beam')
        l_sat = getattr(self,'lsat')
        g_length = getattr(self,'glength')
        n_num =  1+int(getattr(self, 'nsec')/2)
        p_und_l = np.zeros((len(beam_z.I)))
        l_und_l = np.zeros_like(p_und_l)
        p_energy=np.zeros((n_num))
        und_l = np.zeros_like(p_energy)
        for i_ndx, i_el in enumerate(np.linspace(0, float(getattr(self, 'nsec')*(getattr(self, 'nwig') / (getattr(self, 'fl') + getattr(self, 'drl')))), num=n_num)):
            for i_sl in range(len(beam_z.I)):
                if ref[i_sl]==1 and l_sat[i_sl] > i_el:
                    l_und_l[i_sl]=i_el
                else:
                    l_und_l[i_sl]=l_sat[i_sl]
                if ref[i_sl]==1:
                    p_und_l[i_sl]=float((1.0/9.0)*p_noise[i_sl]*np.exp(l_und_l[i_sl]/g_length[i_sl]))
                else:
                    p_und_l[i_sl]=0.0
            p_energy[i_ndx]=float(0.5*np.sum(p_und_l)*getattr(self,'binwidth')/speed_of_light)
            und_l[i_ndx]=i_el*(getattr(self, 'fl') + getattr(self, 'drl'))/getattr(self, 'nwig')
        pulse_energy =float(0.5 * np.sum(getattr(self,'psat'), axis=0)*float(getattr(self,'binwidth'))/speed_of_light)
        setattr(self,'maxpulsee',np.amax(p_energy))
        if getattr(self,'iplot') ==1:
            self.ming_xie_plotting(und_l, p_energy, getattr(self,'psat'), getattr(self,'outplot')+'ming_xie_calculation.png')
        else:
            print('+++++++++++++ No plotting of MingXie ++++++++++')
            with open(getattr(self,'outplot')+'/ming_xie_calculation.txt','a') as f_text:
                
                f_text.write('%.15f\t' % getattr(self,'maxpulsee'))
                f_text.write('%.8f\n' % np.average(getattr(self,'rho'),weights=beam_z.I))
        return p_energy
