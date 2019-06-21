import numpy as np
from ocelot.S2E_STFC.FEL_simulation_block import *
from ocelot.common.globals import *
from copy import copy, deepcopy
import matplotlib.pyplot as plt
import os, shutil
import scipy.special as sf
'''Ming Xie class.

Has the nsec, nwig, fl,drl as attributes (can be read from the GENESIS input file. The Arfken Current is set as a
private attribute.  The user needs to provide the undulator parameter as well, the bin width, and the total number of
particles. The outplot attributes corresponds to the directory where the results are to be stored. If the i_plot
attribute is set to 0, a text file is generated in the outplot/mXie directory showing the maximum value of the pulse
energy calculated via Ming Xie. Otherwise, a plot is generated showing the saturation power as a function of z,
pulse energy as a function of z, rho parameter a function of s and the number of particles per slice.

HMCC: 22-01-18 Implementation of the Ming-Xie class
HMCC: 05-03-18: Bug fixing of particles per slice plot
HMCC: 16-08-18: Implementation of the correction of the Ming Xie scaling due to bunch length (S. Bajlekov, Ph.D. thesis,
                University of Oxford (2011) and Seggebrock, T. et al. Phys. Rev. ST Accel.Beams 16,070703(2013)
HMCC: 03-08-18: Add cooperation length calculation per slice
HMCC: 19-03-19  Add brightness calculation per slice
HMCC: 29-03-19  Redefinition of the  brightness using Saldin et al. (NJP 12(2010) 035010)
                formalism (coherence properties). Bugs in the calculation of the cooperation
                length and saturation length fixed (still bug of the first slice present).
'''

class MingXie(object):
    __lsat = []
    __beam = None
    __rho = []
    __psat = []
    __glength = []
    __npart = []
    __clength=[]
    __brightness={}

    def _init__(self, **initial_value):
        for ke_y, value in initial_value.iteritems():
            setattr(self, ke_y, value)
        self.__constants = {'a1': 0.45, 'a2': 0.57, 'a3': 0.55, 'a4': 1.6, 'a5': 3, 'a6': 2,
                          'a7': 0.35, 'a8': 2.9, 'a9': 2.4, 'a10': 51, 'a11': 0.95, 'a12': 3,
                          'a13': 5.4, 'a14': 0.7, 'a15': 1.9, 'a16': 1140, 'a17': 2.2, 'a18': 2.9, 'a19': 3.2}
        self.__consteta = {'b1':16.7512,'b2':-3.0430,'b3':0.3267}
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
        self.__igen =0


    @property
    def curra(self):
        return 17.045e3

    @property
    def constants(self):
        return {'a1': 0.45, 'a2': 0.57, 'a3': 0.55, 'a4': 1.6, 'a5': 3, 'a6': 2,
                'a7': 0.35, 'a8': 2.9, 'a9': 2.4, 'a10': 51, 'a11': 0.95, 'a12': 3,
                'a13': 5.4, 'a14': 0.7, 'a15': 1.9, 'a16': 1140, 'a17': 2.2, 'a18': 2.9, 'a19': 3.2}

    @property
    def consteta(self):
        return {'b1':16.7512,'b2':-3.0430,'b3':0.3267}

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
    def clength(self):
        return self.__clength

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

    @property
    def igen(self):
        return self.__igen

    @property
    def brightness(self):
        return self.__brightness

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

    @clength.setter
    def clength(self, c_l):
        self.__clength = deepcopy(c_l)

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

    @igen.setter
    def igen(self,i_gen):
        self.__igen = i_gen

    @brightness.setter
    def brightness(self,brw):
        self.__brightness = brw

    def slice_rho(self, sl):
        beam_p = getattr(self, 'beam')
        lmd_u = getattr(self, 'lambdu')
        cur_peak = beam_p.I[sl]
        gamma = beam_p.g0[sl]
        ex = beam_p.ex[sl]
        ey = beam_p.ey[sl]
        betax = beam_p.betax[sl]
        betay = beam_p.betay[sl]
        xrms = np.sqrt(betax*ex)
        yrms = np.sqrt(betay*ey)
        sigma = float((xrms+yrms)/2.0)

        if self.i_gen==0:
            A_Wf = self.K;
        else:
            zeta = np.square(self.K)/(1.0+np.square(self.K))/2.0;
            A_Wf = (sf.j0(zeta)-sf.j1(zeta))*self.K*np.sqrt(2.0);

#        if self.i_gen==0:
#            fc = 1
#        else:
#            a= np.square(self.K)/(2.0*(2.0+np.square(self.K)))
#            fc = sf.j0(a)-sf.j1(a)
        return np.cbrt((cur_peak/self.curra)*np.square(lmd_u*A_Wf/(8.0*pi*sigma)))/gamma

        #return np.cbrt(float(cur_peak/self.curra)*np.square(lmd_u*self.K*fc/(2*np.pi*sigma))*np.power(1.0/(2.0*gamma),3))

    def slice_npart(self,sl):
        beam_p = getattr(self, 'beam')
        curp = beam_p.I[sl]
        t_pulse = np.amax(beam_p.z)-np.amin(beam_p.z)
        nl_wv = getattr(self,'npartmax')*getattr(self,'wavelength')
        q_t = np.trapz(beam_p.I,x= np.array(beam_p.z)/speed_of_light)
        #fact_den = q_t/(nl_wv)
        return int(np.floor(curp*self.binwidth/(q_e*speed_of_light)))

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
        xrms = np.sqrt(betax*ex)
        yrms = np.sqrt(betay*ey)
        sigma = float((xrms+yrms)/2.0)
        emit = float(((ex + ey) / 2.0))
        n_d = (lmb_u / (4.0 * np.pi * np.sqrt(3) * self.slice_rho(sl))) / ((4.0*np.pi *np.square(sigma))/lmb)
        n_e =float((lmb_u / (4.0 * np.pi * np.sqrt(3) * self.slice_rho(sl))) /(np.square(sigma)/emit))*float(4*np.pi*emit/lmb)
        ngam = float(4.0 * np.pi) * float(1.0 / (4.0 * np.pi * np.sqrt(3) * self.slice_rho(sl))) * delgam / gamma
        eta =  (a_cof['a1'] * (n_d**a_cof['a2'])) + (a_cof['a3'] * (n_e**a_cof['a4'])) + (a_cof['a5'] * (ngam**a_cof['a6'])) +\
            (a_cof['a7'] * (n_e**a_cof['a8']) * (ngam**a_cof['a9'])) +(a_cof['a10'] * (n_d**a_cof['a11']) * (ngam**a_cof['a12'])) + \
            (a_cof['a13'] * (n_d**a_cof['a14']) * (n_e**a_cof['a15'])) + \
            (a_cof['a16'] * (n_d**a_cof['a17']) * (n_e**a_cof['a18']) *(ngam**a_cof['a19']))
        return eta

    def slice_eta_z(self,sl):
        beam_p=getattr(self,'beam')
        b_cof = getattr(self,'consteta')
        l_c = getattr(self,'wavelength')/(4.0*np.pi*np.sqrt(3.0)*self.slice_rho(sl))
        #t_pulse = np.mean(np.square(beam_p.z[:sl]))
        t_pulse = beam_p.z[sl]-np.amin(beam_p.z)
        eta_z= b_cof['b1']*np.exp(b_cof['b2']*(np.power(t_pulse,b_cof['b3'])/np.power(l_c,b_cof['b3'])))
        return eta_z

    def slice_p_sat(self, sl):
        beam_p = getattr(self, 'beam')
        curpeak = beam_p.I[sl]
        gamma = beam_p.g0[sl]
        eta =  self.slice_eta(sl)
        eta_z = self.slice_eta_z(sl)
        return 1.6e6 * self.slice_rho(sl) * (np.square(1.0 / ((1.0)*(1.0 + eta)))) * curpeak * (gamma * m_e_MeV)

    def slice_p_noise(self, sl):
        beam_p = self.beam
        gamma = beam_p.g0[sl]
        lbd = self.wavelength
        return float(1e6*speed_of_light * m_e_MeV * q_e * gamma * np.square(self.slice_rho(sl)) / lbd)

    def slice_gain_length(self, sl):
        lbd_u = float(getattr(self, 'lambdu'))
        eta_z = self.slice_eta_z(sl)
        return float((lbd_u*(1.0) * (1.0 + self.slice_eta(sl)))/ (4.0 * np.pi * np.sqrt(3.0) * self.slice_rho(sl)))

    def slice_coop_length(self, sl):
        lbd_u = float(getattr(self, 'lambdu'))
        lg = self.slice_gain_length(sl)
        return float(getattr(self,'wavelength')*lg/lbd_u)

    def slice_sat_length(self, sl):
        return self.slice_gain_length(sl) * np.log(9.0 * self.slice_p_sat(sl) / self.slice_p_noise(sl))

    def slice_brightness_over_lsat(self,sl):
        beam_p = getattr(self, 'beam')
        lmd_u = getattr(self, 'lambdu')
        cur_peak = beam_p.I[sl]
        gamma = beam_p.g0[sl]
        ex = beam_p.ex[sl]
        ey = beam_p.ey[sl]
        ex_scale = 2*np.pi*ex/(self.wavelength)
        ey_scale = 2*np.pi*ey/(self.wavelength)
        betax = beam_p.betax[sl]
        betay = beam_p.betay[sl]
        xrms = np.sqrt(betax*ex)
        yrms = np.sqrt(betay*ey)
        sigma = float((xrms+yrms)/2.0)
        p_sat= self.slice_p_sat(sl)
        Nc = cur_peak/(q_e*self.slice_rho(sl)*(2.0*np.pi*speed_of_light/self.wavelength))

        if self.i_gen==0:
            A_Wf = self.K
        else:
            zeta = np.square(self.K)/(1.0+np.square(self.K))/2.0
            A_Wf = (sf.j0(zeta)-sf.j1(zeta))*self.K*np.sqrt(2.0)

        #if self.i_gen==0:
        #    fc = 1
        #else:
        #    a= np.square(self.K)/(2.0*(2.0+np.square(self.K)))
        #    fc = sf.j0(a)-sf.j1(a)
        gain_par =  np.sqrt(8*cur_peak*np.square(np.pi*A_Wf)/(self.curra*self.wavelength*lmd_u*np.power(gamma,3)))
        D_par = 4*np.pi*gain_par*np.square(sigma)/self.wavelength
        rho_scale =np.cbrt(D_par)*self.slice_rho(sl)
        transv_coh_x = 1.1*np.power(ex_scale,0.25)/(1.0+(0.15*np.power(ex_scale,float(9.0/4.0))))
        transv_coh_y = 1.1*np.power(ey_scale,0.25)/(1.0+(0.15*np.power(ex_scale,float(9.0/4.0))))
        #coh_time = (getattr(self,'wavelength')/(2.0*np.pi*speed_of_light*rho_scale))*np.sqrt(np.pi*np.log(self.wavelength*cur_peak/(2.0*np.pi*rho_scale*q_e*speed_of_light)))
        coh_time = 1./(self.slice_rho(sl)*(2.0*np.pi*speed_of_light/self.wavelength))*np.sqrt(np.pi*np.log(Nc)/18)
        photon_flux = p_sat/(h_J_s*speed_of_light/self.wavelength)

        #photon_flux = 2.0*p_sat*coh_time*np.sqrt(np.pi)/h_J_s
        delta_x_par = photon_flux*coh_time*transv_coh_x
        delta_y_par = photon_flux*coh_time*transv_coh_y
        delta_x_hat = delta_x_par*(0.34*np.pi*speed_of_light/self.wavelength)*(rho_scale*ex_scale)/(photon_flux)
        delta_y_hat = delta_y_par*(0.34*np.pi*speed_of_light/self.wavelength)*(rho_scale*ey_scale)/(photon_flux)
        return 4.0e-15*np.sqrt(2)*speed_of_light*delta_x_par/(np.power(getattr(self,'wavelength'),3)),\
            4.0e-15*np.sqrt(2)*speed_of_light*delta_y_par/(np.power(getattr(self,'wavelength'),3))

    def calculate_brightness_over_lsat(self):
        br_lsat_x= []
        br_lsat_y = []
        beam_z = getattr(self,'beam')
        l_sat = getattr(self,'lsat')
        g_length = getattr(self,'glength')
        n_num =  1+int(getattr(self, 'nsec')/2)
        p_und_l = np.zeros((len(beam_z.I)))
        l_und_l = np.zeros_like(p_und_l)
        p_energy=np.zeros((n_num))
        und_l = np.zeros_like(p_energy)
        for i_sl in range(len(beam_z.z)):
            br_x,br_y = self.slice_brightness_over_lsat(i_sl)
            br_lsat_x.append(br_x/self.slice_sat_length(i_sl))
            br_lsat_y.append(br_y/self.slice_sat_length(i_sl))
        br_lsat_x = np.array(br_lsat_x)
        br_lsat_y = np.array(br_lsat_y)
        setattr(self,'brightness',{'br_x':br_lsat_x, 'br_y':br_lsat_y})
        return

    def ming_xie_plotting(self, z, p_energ, p_noise, filename):
        import matplotlib
        from matplotlib.ticker import ScalarFormatter

        matplotlib.rcParams.update({'font.size':8})
        frmt = ScalarFormatter(useOffset = False)
        frmt.set_scientific(True)

        beam_z = getattr(self,'beam')
        setattr(self,'outplot',filename)
        fig , ax_mxie= plt.subplots(2,2,squeeze=False)
        fig.set_figheight(7)
        fig.set_figwidth(7)
        ln1 = ax_mxie[0,0].plot(z,1e6*p_energ,color='navy')
        ax_mxie[0,0].scatter(z, 1e6 * p_energ, color='navy')
        ax_mxie[0,0].set_ylabel(r'E$_{pulse}$[$\mu$J]',color='navy')
        ax_mxie[0,0].set_xlabel(r'z[m]')
        ax_mxie[0,0].grid(True)
        for ticks in ax_mxie[0,0].xaxis.get_major_ticks():
            ticks.label.set_fontsize('large')
        for ticks in ax_mxie[0,0].yaxis.get_major_ticks():
            ticks.label.set_fontsize('large')
        ax_mxie[0,0].yaxis.get_offset_text().set_fontsize('medium')
        ax_mxie[0,0].set_ylim(ymin=0.95e6*np.amin(p_energ),ymax=1.01e6*np.amax(p_energ))

        ln01 = ax_mxie[0,1].plot(1e6*np.array(beam_z.z),getattr(self,'rho'),'--', color = 'forestgreen')
        ax_mxie[0,1].set_xlabel(r's[$\mu$m]')
        ax_mxie[0,1].set_ylabel(r'$\rho$',color='forestgreen')
        ax_mxie[0,1].tick_params(axis='y',which='both',color='forestgreen',labelcolor = 'forestgreen')
        ax_mxie[0,1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        ax_mxie[0,1].set_ylim(ymin=0.95*np.amin(self.rho),ymax=1.01*np.amax(self.rho))
        ax_mxie[0,1].yaxis.get_offset_text().set_fontsize('medium')
        ax_mxie[0,1].yaxis.get_offset_text().set_color('forestgreen')
        ax_gl = ax_mxie[0,1].twinx()
        ln02 = ax_gl.plot(1e6*np.array(beam_z.z),getattr(self,'glength'),'-.',color='orangered')
        ax_gl.set_ylabel(r'L$_g$[m]',color='orangered')
        ax_gl.tick_params(axis='y', which='both',color='orangered',labelcolor='orangered')
        for ticks in ax_mxie[0,1].xaxis.get_majorticklabels():
            ticks.set_fontsize('large')
        for ticks in ax_mxie[0,1].yaxis.get_majorticklabels():
            ticks.set_fontsize('large')
        for ticks in ax_gl.yaxis.get_majorticklabels():
            ticks.set_fontsize('large')
        ax_gl.yaxis.offsetText.set_fontsize('medium')
        ax_gl.set_ylim(ymin=0.7*np.amin(self.glength),ymax=2*np.amin(self.glength))

        ax_mxie[1,0].plot(1e6*np.array(beam_z.z),1e-9*p_noise,'--',color='brown')
        ax_mxie[1,0].scatter(1e6*np.array(beam_z.z),1e-9*p_noise,color='brown')
        ax_mxie[1,0].set_xlabel(r's[$\mu$m]')
        ax_mxie[1,0].set_ylabel(r'P$_{saturation}$[GW]', color='brown')
        ax_mxie[1,0].tick_params(axis='both',which='both', color='brown',labelcolor='brown')
        for ticks in ax_mxie[1,0].xaxis.get_major_ticks():
            ticks.label.set_fontsize(8)
        for ticks in ax_mxie[1,0].yaxis.get_major_ticks():
            ticks.label.set_fontsize(8)
        ax_mxie[1,0].yaxis.offsetText.set_fontsize(8)
        ax_mxie[1,0].set_ylim(ymin=0.95e-9*np.amin(p_noise),ymax=1.01e-9*np.amax(p_noise))
        ax_mxie[1,0].grid(True,color='maroon')
        ax_mxie[1,1].plot(1e6*np.array(beam_z.z),1e6*self.clength,'--',color='navy')
        ax_mxie[1,1].set_ylabel(r'L$_{cooperation}[\mu$m$]$',color='navy')
        ax_mxie[1,1].tick_params(axis='y', which='both',labelcolor='navy',color='navy')
        ax_mxie[1,1].set_xlabel(r's[$\mu$m]')
        for ticks in ax_mxie[1,1].xaxis.get_major_ticks():
            ticks.label.set_fontsize('large')
        for ticks in ax_mxie[1,1].yaxis.get_major_ticks():
            ticks.label.set_fontsize('large')
        #ax_mxie[1,1].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
        ax_mxie[1,1].xaxis.get_offset_text().set_fontsize('medium')
        ax_mxie[1,1].yaxis.get_offset_text().set_fontsize('medium')
        ax_mxie[1,1].set_xlim(0,1.01e6*np.amax(np.array(beam_z.z)))
        ax_lsat = ax_mxie[1,1].twinx()
        ln_n = ax_lsat.plot(1e6*np.array(beam_z.z),getattr(self,'lsat'),'-.',color='orchid')
        ax_lsat.scatter(1e6*np.array(beam_z.z),self.lsat,color='orchid')
        for ticks in ax_lsat.yaxis.get_major_ticks():
            ticks.label.set_fontsize('large')
        for ticks in ax_lsat.xaxis.get_major_ticks():
            ticks.label.set_fontsize('large')
        ax_lsat.set_ylabel(r'L$_{sat}$[m]',color='orchid')
        ax_lsat.tick_params(axis='y',which='both',color='orchid', labelcolor='orchid',labelsize='large')
        ax_lsat.set_ylim(ymin=0.7*np.amin(self.lsat),ymax=2*np.amin(self.lsat))
        plt.tight_layout()
        fig.savefig(filename,format='png',dpi=120,bbox_inches='tight')
        return fig

    def saturation_power_length(self):
        l_sat = []
        r_ho = []
        n_part = []
        p_sat = []
        p_noise = []
        n_slice = []
        lg = []
        clg= []
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
                clg.append(self.slice_coop_length(i_sl))
                l_sat.append(self.slice_sat_length(i_sl))
                ref.append(1.0)
                r_ho.append(self.slice_rho(i_sl))
            else:
                p_sat.append(0)
                p_noise.append(0)
                n_slice.append(0)
                lg.append(0)
                clg.append(0)
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
        setattr(self, 'clength', np.array(clg))
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
