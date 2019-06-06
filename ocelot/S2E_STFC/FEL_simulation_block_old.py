#!/usr/bin/python
'''FEL simulation block.

HMCC: 25-05-17 Reimplementation of the FEL simulation block (creation of class with corresponding methods)
HMCC: 04-07-17 Fixing bugs of time window (reading nslice from input file)
HMCC: 12-07-17: Fixing bugs reading from input file, adding support to read beam file and fixing bugs 
               of reading existent beam file.
HMCC: 21-07-17: Fixing bug related to centering when an ASTRA file is read. Fixing quad scan bug. Add slice plot when
                a distribution is given (via ASTRA or existent dist file)
HMCC: 29-01-17: Fixing bugs in the reading routine in order to run the S2E for the afterburner scheme
HMCC: 23-03-17: Fixing a bug when there are several noise realisations in order to avoid having the same seed
HMCC: 26-10-18: Workaround for beta scan (fixed for aw0)
HMCC: 26-06-18: Add tapering (wcoefz) support
HMCC: 26-10-18: Add support for fresh slice. 
'''
#################################################
### import of all modules that are required.
from __future__ import print_function
from copy import deepcopy,copy
from ocelot import *
from ocelot.utils.xfel_utils import *
from ocelot.gui.accelerator import *
from ocelot.gui.genesis_plot import *
from ocelot.optics.elements import Filter_freq
from ocelot.adaptors.genesis import *
from ocelot.common import globals
from ocelot.adaptors.genesis import generate_input,generate_lattice,get_genesis_launcher,run_genesis,filename_from_path, read_edist_file
from ocelot.gui.genesis_plot import plot_gen_out_all,plot_gen_stat
import ocelot.cpbd.elements
import numpy as np
import matplotlib.pyplot as plt
import os,sys, errno, shutil
import time
    

class FEL_simulation_block(object):
    '''Class which corresponds to the FEL simulation block 
     The instance receives a data structure as 
     The FEL simulation block class gets a "data" dictionary, with the following entries:
         gen_file: path where the genesis input file is
         file_pout: Path where the results are to be saved
         stat_run: Number of noise realisations
         gen_launch: GENESIS launcher (Could be version 2 or version 3)
         i_edist: Flag to set if there is an existent dist file.
         file_edist: Path where the existent dist file is located.
         i_astra: Flag in order to read an ASTRA file (bunch)
         astra_file: Path corresponding to the astra file
         i_beam: Flag to set if there is an existent beam file
         file_beam: Path where the existent beam file is located.
         i_rad: Flag to set if there is an existent rad file.
         file_rad: Path where the existent rad file is located.
         i_dpa: Flag to set if there is an existent dpa file
         file_dpa: Path where the existent dpa file is located
         i_dfl: Flag to set if there is an existent dfl file
         file_dfl: Path where the existent dfl file is located
         i_match: Flag in order to match the distribution
         tw_match: Dictionary which contains the twiss parameters to match upon:
             - alpha_x
             - alpha_y
             - beta_x
             - beta_y
         i_rewrite: Flag to rewrite any input parameter that sets up the simulation 
                    (GENESIS simulation parameters, lattice structure)
         par_rew: Dictionary which contain the    parameters of the input file to be overwritten.
         idump: Flag to dump dpa and fld files
         i_scan : flag to determine whether a scan is going to take place
         parameter: Parameter scan        
         init: start of the scan parameter
         end: End value of the scan parameter
         n_scan : Number of scans
         
    The instance is made by setting up a dictionary with the necessary attributes
        FEL_simulation_block(data)
    '''
    def __init__(self, *initial_data, **kwargs):
        for dictionary in initial_data:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])
    
    def read_GEN_input_file(self):
        from ocelot.adaptors.genesis import GenesisInput
        
        A_input = GenesisInput()
        with open(self.gen_file, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                else:
                    splitLine = line.rsplit('=')
                    if splitLine[0][0]!=' ':
                        splitLine[0] = ' '+splitLine[0]
                    if splitLine[0].startswith('alphax'):
                        splitLine[0]=splitLine[0].replace('=','').rstrip()
                    else:
                        splitLine[0]=splitLine[0][1:].replace('=','').rstrip()         
                    splitLine[-1] = splitLine[-1].replace('D+','e').replace('D-','e-').replace('\n','').replace('\r','').rstrip()
                    splitLine[-1] = splitLine[-1].replace('E+','e').replace('E-','e-').replace('\n','').replace('\r','').rstrip()
                    if(splitLine[0][:].endswith('file') and splitLine[0][:].startswith('mag')):
                        val_attr = str(splitLine[-1])
                    elif (splitLine[0].startswith('$n')) or (splitLine[0].find('filetype')!=-1) or (splitLine[0].find('distfile')!=-1) or (splitLine[0].find('beamfile')!=-1):
                        continue
                    elif(splitLine[0].startswith('$end')):
                        break
                    elif ((str(splitLine[0]).startswith('i') and not str(splitLine[0])=='ibfield' and not str(splitLine[0])=='imagl') or str(splitLine[0])=='itdp' or str(splitLine[0]).startswith('n') 
                          or str(splitLine[0]).startswith('lbc') or str(splitLine[0]).startswith('magin') or str(splitLine[0]).startswith('magout')   
                          or str(splitLine[0]).startswith('ffspec') or str(splitLine[0]).startswith('convharm') or str(splitLine[0]).startswith('ippart')
                          or str(splitLine[0]).startswith('ispart') or str(splitLine[0]).startswith('ipradi') or str(splitLine[0]).startswith('isradi')
                          or str(splitLine[0])=='alignradf'
                          or str(splitLine[0]).startswith('multconv') or str(splitLine[0])=='offsetradf'):
                        val_attr = int(float(splitLine[-1].replace('=',"")))
                    elif (str(splitLine[0]).startswith('fbess')):
                        val_attr = int(float(splitLine[-1].replace('=',"")))
                    elif (str(splitLine[0]).startswith('ibfield')) or (str(splitLine[0]).startswith('imagl')) :
                        val_attr = float(splitLine[-1].replace('=',""))
                    elif (str(splitLine[0])=='xkx') or (str(splitLine[0])=='xky'):
                        val_attr = float(splitLine[-1].replace('=',""))
                    elif (str(splitLine[0])=='f1st'):
                        val_attr = int(float(splitLine[-1].replace('=',"")))
                    elif (splitLine[0].startswith('outputfile')):
                        continue
                    elif (splitLine[0].startswith('wcoefz')):
                        val_attr = [float(splitLine[-1][2:].rsplit(' ')[i]) for i in range(0,7,3)] 
                    elif (splitLine[0].startswith('lout')):  
                        val_attr = [int(float(l_out)) for l_out in splitLine[-1][1:].rsplit(' ')]
                    elif(splitLine[0].startswith('alpha')):
                        val_attr = float(splitLine[-1].replace('=',""))
                    elif((splitLine[0].startswith('magin')) and (splitLine[0].endswith('file'))):
                        val_attr = str(splitLine[-1].replace('=',""))
                    elif (splitLine[0].startswith('partfile')) or ((splitLine[0].startswith('fieldfile'))):
                        continue
                    else:
                        val_attr = float(splitLine[-1].replace('=',""))        
                    setattr(A_input,str(splitLine[0]),val_attr)
                    if hasattr(A_input,'maginfile') and A_input.maginfile != 0:
                        vfile = str(A_input.maginfile).replace(' ',"").replace("'","")
                        setattr(A_input,'latticefile',vfile)
                    elif not hasattr(A_input,'maginfile'):
                        setattr(A_input,'latticefile',None)
        f.close()
        return A_input
        
    def undulator_design(self,A_contents,i_fs=0):
        from ocelot.cpbd.elements import Drift, Quadrupole, Undulator
        from ocelot.cpbd.magnetic_lattice import MagneticLattice
        from ocelot.common.globals import m_e_GeV, speed_of_light, h_eV_s
        from ocelot.rad.undulator_params import Ephoton2K
        from ocelot.rad.undulator_params import UndulatorParameters
        
        ## Taking it from the Notebook
        xlamd = getattr(A_contents,'xlamd')
        nwig = getattr(A_contents,'nwig')
        gamma0 = getattr(A_contents,'gamma0')
        xlamds = getattr(A_contents,'xlamds')
        fl = getattr(A_contents,'fl')
        dl = getattr(A_contents,'dl')
        drl = getattr(A_contents,'drl')
        quadf = getattr(A_contents,'quadf')
        quadd = -getattr(A_contents,'quadd')
        nsec = getattr(A_contents,'nsec') 
        f1st = getattr(A_contents,'f1st')
        E_beam = gamma0*m_e_GeV
        E_photon=h_eV_s*speed_of_light/xlamds
        p_beam = np.sqrt(E_beam**2-m_e_GeV**2)
        drl1 = int(np.rint((drl-nwig)*0.5))
        drl2 = drl-nwig-drl1
        kw0 = float(getattr(A_contents,'aw0'))/np.sqrt(0.5)
       
       # Instance of the Undulator object
        und= Undulator(lperiod=xlamd, nperiods=nwig, Kx=kw0)
       
       # Calculation of the Undulator parameter from the Photon and Beam Energies)
        und.Kx = kw0    
       
       # Drift sections (if they are the same)
        d_rift = Drift(l=drl1*und.lperiod)
        d_rift2 = Drift(l=drl2*und.lperiod)
       
       # Definition of Quads
     
        qf= Quadrupole(l=fl*und.lperiod,k1=0.3*quadf/p_beam)
        qd = Quadrupole(l=dl*und.lperiod,k1=0.3*quadd/p_beam)
        if i_fs ==0:
            qdh=deepcopy(qf)
        else:
            qdh=deepcopy(qd)
        qdh.l/=2
        
       # Creating the cell
     
        extra_fodo = (und,d_rift,qdh)
        if i_fs ==0:
            cell_ps = (und,d_rift, qd, d_rift2,und, d_rift,qf,d_rift2)
        else:
            cell_ps = (und,d_rift, qf, d_rift2,und, d_rift,qd,d_rift2) 
        l_fodo = MagneticLattice(cell_ps).totalLen/2 ##Length of fodo cell
        if f1st ==1:
            sase3= MagneticLattice((qdh,d_rift,)+(np.ceil(nsec/2)*cell_ps)) # Lattice with nsec modules
        elif f1st ==0 and nsec >1:
            sase3=  MagneticLattice(np.rint(nsec/2)*cell_ps) # Lattice with nsec modules
        elif f1st ==0 and nsec ==1:
            sase3=  MagneticLattice(cell_ps) # Lattice with nsec modules
        up = UndulatorParameters(und,E_beam) # Instance of the Class UndulatorParameters
        print('++++ Undulator Parameters +++')
        up.printParameters()
        return {'Undulator Parameters':up,'Magnetic Lattice':sase3}
        
    def BeamDefinition(self,A_contents):
        from ocelot.common.globals import m_e_GeV, speed_of_light
        from ocelot.cpbd.beam import Beam
    
        beamf = Beam()
        A_dict= {'E':getattr(A_contents,'gamma0')*m_e_GeV,'sigma_E':getattr(A_contents,'delgam')*m_e_GeV,'beta_x':getattr(A_contents,'gamma0')*(getattr(A_contents,'rxbeam')**2)/(getattr(A_contents,'emitx')), 
            'beta_y':getattr(A_contents,'gamma0')*(getattr(A_contents,'rybeam')**2)/getattr(A_contents,'emity'), 'alpha_x':getattr(A_contents,'alphax'),'alpha_y':getattr(A_contents,'alphay'),
            'emit_x':getattr(A_contents,'emitx')/getattr(A_contents,'gamma0'),'emit_y' : getattr(A_contents,'emity')/getattr(A_contents,'gamma0'),'emit_xn':getattr(A_contents,'emitx'),'emit_yn':getattr(A_contents,'emity'),
             'x' :  0.000000e+00,'y' : 0.000000e+00,'px':0,'py':0,'I':getattr(A_contents,'curpeak'),'tpulse':1e15*(getattr(A_contents,'curlen')/speed_of_light)}
                
        for item in A_dict:
            setattr(beamf, item,A_dict[item])
        return beamf
    
    def beta_matching(self,inp,f_path):
        import os
        import shutil
        from ocelot.adaptors.genesis import filename_from_path
        from ocelot.common.globals import m_e_GeV
        from copy import deepcopy
        file_or = getattr(self,'gen_file')
        A_params = ['rxbeam', 'rybeam','alphax','alphay','emitx','emity']
        inp0 = inp
        os.system('cp /scratch2b/qfi29231/betamatch_dir/betamatch %s' %(f_path))
        os.chdir(f_path) 
        with open(f_path+'/mod_file.in','w') as f:
            f.write(inp0.input())
        f.close()
        with open(f_path+'/beta_input_file.in','w') as f:
            f.write('mod_file.in')
        f.close()
       
        if inp.latticefile != None or inp0.lat != None:
            with open(f_path+'/'+inp0.latticefile,'w') as f:
                f.write(generate_lattice(inp0.lat,unit=inp0.xlamd*inp0.delz,energy=inp0.gamma0*m_e_GeV))   
                f.close()
        try:
            os.system('./betamatch < beta_input_file.in')
        except:
            print('Betamatch did not work')
            
        setattr(self,'gen_file',f_path+'TEMPLATE.IN')
            
        inp2 = self.read_GEN_input_file()
        for params in A_params:
            if getattr(inp0,params)==getattr(inp2,params) and  (getattr(inp0,params)==getattr(inp0,params)):
                print('Betamatch not change parameters')
            elif (getattr(inp0,params)==getattr(inp0,params)) and getattr(inp0,params)!=getattr(inp2,params) and getattr(inp0,params)!='nan':
                print('Betamatch changed parameters')
                setattr(inp0,params,getattr(inp2,params))
            elif (getattr(inp0,params)!=getattr(inp0,params)):
                print('Nan as a solution. Parameters of original input file')
                for params2 in A_params:
                    setattr(inp0, params, getattr(inp2, params))
                break
        os.remove(f_path+'/beta_input_file.in')
        os.remove(f_path+'/TEMPLATE.IN')
        os.remove(f_path+'/mod_file.in')
        setattr(self,'gen_file',file_or)
        return inp0
        
    def rematch_edist(self,inp):
        from ocelot.adaptors.genesis import rematch_edist
        from ocelot.cpbd.beam import Twiss
        
        tw = Twiss()        
        inp_out = inp
        for key,value in getattr(self,'tw_match').items():
            setattr(tw,key,value)
        e_distf = rematch_edist(getattr(inp_out,'edist'),tw)
        setattr(e_distf,'filePath',getattr(self,'file_pout')+'rematched_edist')    
        setattr(inp_out,'edist',e_distf)
        return inp_out
    
    def convert_ASTRA_edist(self,inp):       
        from ocelot.adaptors.genesis import read_astra_dist, astra2edist
        
        inp_out = inp
        edist = read_astra_dist(getattr(self,'astra_file'))
        edist = astra2edist(edist, center=1)
        setattr(edist,'filePath',getattr(self,'file_pout')+'read_edist_astra')
        setattr(inp_out,'edist',edist)
        return inp_out 

    def GEN_rewrite_par(self,in_p):
        inp = in_p
        for key in getattr(self,'par_rew'): 
            setattr(inp,key, self.par_rew[key])
        return inp 
    
    def GEN_existent_beam_dist_dpa_rad(self,in_p,flag_bd):
        from ocelot.adaptors.genesis import read_beam_file, read_edist_file, read_rad_file, read_dpa_file,read_dfl_file

        inp = in_p
        if flag_bd == 'beam':
            bfd = read_beam_file(getattr(self,'file_beam')) 
            setattr(bfd,'filePath',getattr(self,'file_pout')+'existent_beam')
            setattr(inp,'beam',bfd)
            return inp
        elif flag_bd == 'edist':
            bfd = read_edist_file(getattr(self,'file_edist'))
            setattr(bfd,'filePath',getattr(self,'file_pout')+'existent_edist')
            setattr(inp,'edist',bfd)
            return inp
        elif flag_bd == 'rad':
            bfd = read_rad_file(getattr(self,'file_rad')) 
            setattr(inp,'rad',bfd)
            return inp
        elif flag_bd == 'dpa':
            bfd = read_dpa_file(getattr(self,'file_dpa'),nbins = getattr(inp,'nbins'),npart = getattr(inp,'npart'))
            setattr(inp,'dpa',bfd)
            return inp 
        elif flag_bd == 'dfl':
            if getattr(inp,'dgrid')==0:
                rbeam = sqrt((getattr(inp,'rxbeam')**2)+(getattr(inp,'rybeam')**2))
                ray = sqrt((getattr(inp,'zrayl')*getattr(inp,'xlamds')/np.pi)*(1+((getattr(inp,'zwaist')**2)/(getattr(inp,'zrayl')**2))))
                d_len = float(getattr(inp,'rmax0'))*(ray+rbeam)
            else:
                d_len = float(2*getattr(inp,'dgrid'))
            bfd = read_dfl_file(getattr(self,'file_dfl'),Nxy=getattr(inp,'ncar'),Lxy = d_len,zsep=getattr(inp,'zsep'),xlamds=getattr(inp,'xlamds'))
            setattr(inp,'dfl',bfd)
            return inp
        else:
            pass

    def GEN_scan(self,n_p ,A_inp,A_und,in_p): 
        from ocelot.cpbd.beam import Twiss
        from ocelot.gui.genesis_plot import plot_edist  
        inp = in_p
        A_inpt = A_inp
        inp_file = inp.run_dir + 'run.' + str(inp.runid) + '.inp'
        inp_file = filename_from_path(inp_file)

        A_simul = ['npart','ncar','zsep','delz','dmpfld','fbess0','dgrid','rmax0','xkx','xky','iwityp']
        A_undl = ['quadd', 'quadf','fl','dl','drl','nsec','nwig','aw0', 'awd']

        if (self.parameter in A_simul):
            setattr(inp,self.parameter,n_p)
            inp.lat = A_und['Magnetic Lattice']
            setattr(inp,'magin',1)
            print(' ++++++++++ Scan {0} of the parameter {1}'.format(n_p, self.parameter))
        elif ((self.parameter in A_undl) or (self.parameter == 'xlamds') or self.parameter =='aw0'): 
            print(' ++++++++++ Steady State Scan {0} of the parameter {1} Quad optimisation'.format(n_p, self.parameter))
            setattr(inp,'type','steady')
            setattr(inp,'itdp',0)
            setattr(inp,'shotnoise',0)
            setattr(inp,'prad0',10)
            setattr(inp,'betamatch',True)
            if (self.parameter in A_undl):
                if self.parameter == 'aw0':
                    setattr(A_inpt,'aw0',n_p)
                    setattr(A_inpt,'awd',n_p)
                    setattr(inp,'aw0',n_p)
                    setattr(inp,'awd',n_p)
                else:
                    n_par = int(n_p)
                    if self.parameter in A_undl[0:2]:
                        for j_quad in A_undl[0:2]:
                            setattr(A_inpt,str(j_quad),n_par)
                            setattr(inp,str(j_quad),n_par)
                    else:
                        setattr(A_inpt,self.parameter,n_par)
                        setattr(inp,self.parameter,n_par)
            else:
                n_par =float(n_p)
                setattr(A_inpt,self.parameter,n_p)
                setattr(inp,self.parameter,n_p)
            setattr(inp,'magin',0)
            setattr(inp,'lat',None)
            #A_undulator=self.undulator_design(A_inpt)
            #setattr(inp,'lat',A_undulator['Magnetic Lattice'])
            #setattr(inp,'latticefile',inp_file+'.lat')
            #setattr(inp,'magin',1)
            #if self.parameter !='aw0':
            inp = self.beta_matching(inp,inp.run_dir)
            if inp.edist!=None:
                inp0=inp
                for i_tw in getattr(self,'tw_match'):
                    self.tw_match[i_tw]=getattr(inp0,i_tw.replace('_',''))
                inp=self.rematch_edist(inp0)
            inp.latticefile=None                       
        else:
            n_par = float(n_p) 
            setattr(A_inpt,self.parameter,n_p)
            setattr(inp,self.parameter,n_p)
            #setattr(inp,'lat',A_undulator['Magnetic Lattice'])
            #setattr(inp,'magin',1)
        return inp
        
    def gen_outplot_single(self,run_inp= [], itdp = True,savefig=True):
        import matplotlib.pyplot as plt
        from ocelot.gui.genesis_plot import subfig_rad_pow, subfig_rad_pow_evo
        from ocelot.adaptors.genesis import read_out_file
    
        dict_name={'p_int':'radiation power','energy': 'radiation pulse energy','el_e_spread': 'el.beam energy spread','el_energy': 'el.beam energy average','bunching': 'el.beam bunching','spec': 'radiation on-axis spectral density','dfl_spec':'total radiation spectral density','r_size':'radiation transv size','r_size_weighted':'radiation transv size (weighted)','xrms':'el.beam x size','yrms':'el.beam y size','error':'genesis simulation error','p_mid':'radiation power on-axis','phi_mid':'radiation phase on-axis','increment':'radiation power increment'}
        dict_unit={'p_int':'[W]','energy': '[J]','el_e_spread': '(gamma)','el_energy': '(gamma)','bunching': '','spec': '[arb.units]','dfl_spec': '[arb.units]','r_size':'[m]','xrms':'[m]','yrms':'[m]','error':''}
        figsize = (14,7)
        
        #if proj_dir[-1]!='/':
        #    self.file_pout+='/'
        if run_inp==[]:
            run_range=xrange(1000)
        #elif run_inp ==1:
         #   run_range=xrange(1)
        else:
            run_range=run_inp
        if itdp ==True:
            param_inp = ['energy','spec','bunching','el_energy','el_e_spread','xrms','yrms','p_int','dfl_spec','p_mid']
            s_inp=['max','mean']
            z_inp=[0,'end'] 
            n_seeds = len([self.file_pout+'scan_'+str(run_range[0])+'/'+files for files in os.listdir(self.file_pout+'scan_'+str(run_range[0])) if files.startswith('ip_s')])
            n_totl=int(n_seeds)*int(len(run_range))
        else:
            param_inp=['el_e_spread','bunching','el_energy','xrms','yrms','p_int']
            z_inp=[0,'end']
            s_inp=['max','mean'] 
            n_totl = len(run_inp)
        param_range = param_inp
        
        outlist=[]  
        run_range_good=[]
        text_l = []
        plt.rc('font',**{'size':8})
        plt.rcParams['text.usetex']=False
        plt.rcParams['text.latex.unicode']=False
        plt.rcParams.update({'font.size':8})  
    
        n_seeds = len([self.file_pout+'scan_'+str(run_range[0])+'/'+files for files in os.listdir(self.file_pout+'scan_'+str(run_range[0])) if files.startswith('ip_s')])
        colour=plt.cm.Set1(np.linspace(0.1,0.999999,n_totl))
        colour = map(lambda rgb:'#%02x%02x%02x' %(rgb[0]*255,rgb[1]*255,rgb[2]*255),
                     tuple(colour[:,0:-1]))
        fig,ax = plt.subplots(1,2)
        for i_run,runr in enumerate(run_range):
            ip_s=[]
            ip_seed = [self.file_pout+'scan_'+str(runr)+'/'+files for files in os.listdir(self.file_pout+'scan_'+str(runr)) if files.startswith('ip_s')]
            colourr = colour[i_run]
            for i_seed, ipseed in enumerate(ip_seed):
                out_file = [ipseed+'/'+files for files  in os.listdir(ipseed) if ((files.startswith('run.')) and (files.endswith('.gout')))]
                if os.path.isfile(out_file[0]):
                    g=read_out_file(str(out_file[0]),read_level=2)
                    ip_s.append(g)
                    
                text_leg = 'Scan ='+str(runr)+' ip_seed = '+str(ipseed[ipseed.find('ip_seed_')+8:])
                text_l.append(text_leg)
               # colourr = colour[i_seed]
                fig=subfig_rad_pow(ax[0],g,text_leg,colour[i_run],log=1)
                fig = subfig_rad_pow_evo(ax[1],g,text_leg,norm=1)
            run_range_good.append(i_run)
            outlist.append(ip_s)
        run_range=run_range_good
        if savefig==True:
            savefig='png'
            saving_path=self.file_pout+'results/'
            if not os.path.isdir(saving_path):
                os.makedirs(saving_path)
            print('      saving to '+saving_path)        
        art = []
        for plot_ax in ax:
            for tick in plot_ax.yaxis.get_major_ticks():
                tick.label.set_fontsize('8')
                for tick in plot_ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize('8')
                plot_ax.grid(True)
        lgd = ax[0].legend(loc=9,bbox_to_anchor=(0.98,0.4),prop={'size':4.5})
        art.append(lgd)
    
        if savefig!=False:
            print('      saving '+'power_plot'+'.'+savefig)
            plt.savefig(saving_path+'power_plot'+'.'+savefig,format=savefig,dpi=1200,figsize=(1,7), additional_artists=art,bbox_inches='tight')
    
    
        for param in param_range: 
            for s_ind in s_inp:
                s_fig_name='Z__'+'stage_1__'+dict_name.get(param,param).replace(' ','_').replace('.','_')+'__'+str(s_ind)
                fig=plt.figure(s_fig_name)
                fig.clf()
                fig.set_size_inches(figsize,forward=True)
                ax=fig.gca()
                for irun,runr in enumerate(run_range):
                    colourr = colour[irun]
                    text_leg = text_l[irun]
                    for i_seed in xrange(len(outlist[0])):
                        #colourr = colour[i_seed]
                        #text_leg = text_l[i_seed]
                        if not hasattr(outlist[irun][i_seed],param):
                            continue
                        else:
                            param_matrix=deepcopy(getattr(outlist[irun][i_seed],param))
    
                        if len(param_matrix) == len(outlist[irun][i_seed].z):
                            s_value=param_matrix
                        else:
                            if s_ind=='max':
                                s_value=np.amax(param_matrix,axis=0)
                            elif s_ind=='max_cur':
                                s_value=param_matrix[outlist[irun][i_seed].sn_Imax,:]
                            elif s_ind=='mean':
                                 s_value=np.average(param_matrix,axis=0,weights=outlist[irun][i_seed].I)#HMCC doing average current
                            #elif s_ind=='mean':#HMCC comment
                            #    s_value=np.mean(param_matrix,axis=0) #HMCC comment
                            else:
                                si=np.where(outlist[irun][i_seed].s<=s_ind)[-1][-1]
                                s_value=param_matrix[si,:]
                        fig=plt.plot(outlist[irun][i_seed].z,s_value,'0.8',color = colourr, linewidth=1,label=text_leg)
                        
                
                plt.xlabel('z [m]')
                plt.ylabel(dict_name.get(param,param)+' '+dict_unit.get(param,''))
                art = []
                lgd = ax.legend(loc=9,bbox_to_anchor=(1.11,1),prop={'size':4.5})
                art.append(lgd)
                if savefig!=False:
                    print('      saving '+s_fig_name+'.'+savefig)
                    plt.savefig(saving_path+s_fig_name+'.'+savefig,format=savefig,dpi=1200,figsize=(1,7), additional_artists=art,bbox_inches='tight')
    
        for param in param_range:
            for z_ind in z_inp:
                z_value=[]
                z_fig_name='S__'+'stage_1__'+dict_name.get(param,param).replace(' ','_').replace('.','_')+'__'+str(z_ind)+'__m'
                fig=plt.figure(z_fig_name)
                fig.clf()
                fig.set_size_inches(figsize,forward=True)
                ax = fig.gca()
                for irun,runr in enumerate(run_range):
                    colourr = colour[irun]
                    text_leg = text_l[irun]
                    for i_seed in xrange(len(outlist[0])):
                        #colourr = colour[i_seed]
                        #text_leg =text_l[i_seed]
                        if not hasattr(outlist[irun][i_seed],'s') or (itdp==False):
                               print('Steady state run')
                               break
                        elif not hasattr(outlist[irun][i_seed],param):#HMCC
                            break
                        else:
                            param_matrix=deepcopy(getattr(outlist[irun][i_seed],param))#HMCC
    
                        if len(param_matrix) == len(outlist[irun][i_seed].z): #case if the array is 1D (no s/z matrix presented) HMCC
                            break
                        else:
                            if z_ind=='end' or z_ind==inf:
                                z_value=param_matrix[:,-1] #after undulator
                            elif z_ind=='start':
                                z_value=param_matrix[:,0] #before undulator
                            else:
                                zi=np.where(outlist[irun][i_seed].z<=z_ind)[-1][-1]
                                z_value=param_matrix[:,zi]
                            if param=='spec':
                                freq_scale=outlist[irun][i_seed].freq_lamd#*1e9#HMCC
                                fig=plt.plot(freq_scale,z_value,'0.8',color = colourr,label=text_leg)
                                plt.xlabel('$\lambda$ [nm]')
                            else:
                                s_scale=outlist[irun][i_seed].s*1e6#HMCC
                                fig=plt.plot(s_scale,z_value,'0.8',color = colourr,label=text_leg)
                                plt.xlabel('s [um]')
                                plt.ylabel(dict_name.get(param,param)+' '+dict_unit.get(param,''))
                            if savefig!=False:
                                print('      saving '+z_fig_name+'.'+savefig)
                                art = []
                                lgd = ax.legend(loc=9,bbox_to_anchor=(1.11,1),prop={'size':4.5})
                                art.append(lgd)
                                plt.savefig(saving_path+z_fig_name+'.'+savefig,format=savefig,dpi=1200,figsize=(1,7), additional_artists=art,bbox_inches='tight')
        return fig
        
    def post_processing(self,inp,s_scan_inp):  
        from ocelot.adaptors.genesis import edist2beam
        plt.rcParams['text.usetex']=False
        plt.rcParams['text.latex.unicode']=False
        #Plot distribution
        if getattr(inp,'edist')!=None:
            plot_edist(getattr(inp,'edist'),figsize=20,savefig=True,showfig=False,plot_x_y=True,plot_xy_s=False,bins=(250,250,250,250)) 
            plot_edist(getattr(inp,'edist'),figsize=20,savefig=True,showfig=False,plot_x_y=False,plot_xy_s=True,bins=(250,250,250,250))
            bfd_b = edist2beam(getattr(inp,'edist'),step=float(10.0*getattr(inp,'xlamds')))
            setattr(bfd_b,'filePath',self.file_pout+'slice_edist')
            plot_beam(bfd_b,savefig=True,showfig=False)
        elif getattr(inp,'beam')!=None and hasattr(inp.edist,'I') and hasattr(inp.beam,'eloss'):
            plot_beam(getattr(inp,'beam'),figsize=20,savefig=True,showfig=False)
        # Power , electron and radiation plots

        if getattr(self,'i_scan')==0:
                if getattr(self,'stat_run')>1:
                    plot_gen_stat(self.file_pout, run_inp=s_scan_inp, stage_inp=xrange(1), param_inp=[], s_param_inp=['p_int', 'energy','xrms','yrms','bunching','r_size_weighted','spec', 'error','p_mid'], z_param_inp=['p_int', 'phi_mid_disp','energy','el_e_spread','spec','xrms','yrms', 'spec', 'bunching', 'wigner'], dfl_param_inp=[], run_param_inp=[], s_inp=['max','mean'], z_inp=[0,'end'], run_s_inp=[], run_z_inp=[], savefig=True, saveval=False,debug=0)
                elif getattr(self,'stat_run')==1:
                    if getattr(inp,'itdp')==0:
                        i_tdp = False
                    else:
                        i_tdp = True
                    self.gen_outplot_single(run_inp= s_scan_inp, itdp = i_tdp,savefig=True)
        else:
            self.gen_outplot_single(run_inp= s_scan_inp, itdp = False,savefig=True)
        plt.close("all")
        
    def GEN_simul_preproc(self,A_input,i_aft=0,i_fs=0):
        if not self.file_pout.endswith('/'):
            self.file_pout=self.file_pout+'/'
        
        inp_arr = []
        A_bbeam = ['rxbeam','rybeam','emitx','emity','alphax','alphay','xbeam','ybeam','pxbeam','pybeam']
        A_simul = ['alignradf','npart','ncar','zsep','delz','dmpfld','fbess0','dgrid','rmax0','xkx','xky','iwityp',
                   'nptr','lbc','zrayl','zstop','zwaist','delgam','xlamd','nscr','nscz','curpeak',
                   'iphsty','nharm','curlen','nbins','gamma0','isravg','isrsig','eloss','version',
                   'multconv','imagl','convharm','idril','ffspec','ndcut','ibfield','nslice','ntail',
                   'ippart','ispart','ipradi','isradi']
        A_td = ['itdp','prad0','shotnoise']
        A_und = ['quadd', 'quadf','fl','dl','drl','nsec','nwig','aw0', 'awd']

        print('++++ Output Path {0} ++++++'.format(self.file_pout))
    
        # Setting the number of noise realisations and scan (over quads or different parameters)
        
        if (self.i_scan ==0):
            s_scan = range(1)  
            num = self.stat_run
            run_ids = xrange(0,num)
            print('++++++++ No scan ++++++++++')
        elif (self.i_scan !=0):
            if (self.parameter in A_und):
                run_ids= range(1)
                if self.parameter !='aw0':
                    s_scan = range(int(self.init),int(self.end),int(np.ceil((self.end-self.init)/(self.n_scan))))
                else:
                    s_scan = np.linspace(self.init,self.end,self.n_scan) 
                print('++++ Quad scan, parameter  {0} ++++++'.format(self.parameter))
            elif (self.parameter=='xlamds'):
                run_ids= range(1)
                s_scan = np.linspace(self.init,self.end,self.n_scan)
                print('++++ Quad scan, parameter  {0} ++++++'.format(self.parameter))
            else:        
                s_scan = np.linspace(self.init,self.end,self.n_scan)
                num = self.stat_run
                run_ids = xrange(0,num)
                print('++++ Number of noise realisations {0} ++++++'.format(num))
    
            # setting the undulator design( Magnetic Lattice)
        A_undulator = self.undulator_design(A_input,i_fs)
        
            # Fill in the beam object
        A_beam = self.BeamDefinition(A_input)
        if (getattr(A_input,'itdp')==0):
            print('++++ Steady State run +++++')
            i_tdp = False
        elif (getattr(A_input,'itdp')==1):
            i_tdp = True
         
            # Generate input object
        inp = generate_input(A_undulator['Undulator Parameters'],A_beam,itdp=i_tdp)
        
        # Overwrite the simulation attributes of the input object with the ones defined in the input file
        for key in A_input.__dict__:
            if (key in A_simul) or (key in A_und) or (key in A_td) or (key =='xlamds') or (key == 'f1st') or (key == 'nslice') or (key == 'ntail'):
                setattr(inp,key, getattr(A_input,key))
      
        for key in ['edist','beam','dfl']:
            if getattr(A_input,key)!=None:
                setattr(inp,key,getattr(A_input,key))
        
        # Set up some input parameters
        if getattr(inp,'itdp')==0:
            setattr(inp,'type','steady')
        else:
            setattr(inp,'type','tdp')
        setattr(inp, 'awd', float(getattr(inp, 'aw0')))

        # idump attribute
        if (getattr(self,'idump')) == 1:
            setattr(inp,'idump',1)
            setattr(inp,'idmpfld',1)

        # Existent dist or beam file (if exists)
        if (getattr(self,'i_edist') == 1) and (hasattr(self,'file_edist')):
            inp=self.GEN_existent_beam_dist_dpa_rad(inp,'edist')
        elif (getattr(self,'i_beam') == 1) and (hasattr(self,'file_beam')):
            inp=self.GEN_existent_beam_dist_dpa_rad(inp,'beam')
        elif  (getattr(self,'i_rad') == 1) and (hasattr(self,'file_rad')):
            inp=self.GEN_existent_beam_dist_dpa_rad(inp,'rad')
        elif  (getattr(self,'i_dpa') == 1) and (hasattr(self,'file_dpa')):
            inp=self.GEN_existent_beam_dist_dpa_rad(inp,'dpa')
        elif  (getattr(self,'i_dfl') == 1) and (hasattr(self,'file_dfl')):
            inp=self.GEN_existent_beam_dist_dpa_rad(inp,'dfl')
        else:
            print('++++ No edist or beam or dpa or rad file available ++++++')
        
        
        # Read ASTRA file.
        if hasattr(self,'i_astra') and getattr(self,'i_astra')==1 and hasattr(self,'astra_file'):
            inp=self.convert_ASTRA_edist(inp)
            setattr(inp,'beam',None)
        elif  (hasattr(self,'i_astra') and getattr(self,'i_astra')==1) and not (hasattr(self,'astra_file')):
            print('Path of  the ASTRA file not provided')
            return
        else:
            print('No need to read ASTRA file')
        
        # Rematch beam (if the flag has been set within the data dictionary)
        if (getattr(inp,'edist')!= None) and hasattr(self,'i_match') and (getattr(self,'i_match')==1):
            inp = self.rematch_edist(inp)  

        # Setting up the time window of the distribution
        if getattr(self,'i_beam')==0:
            setattr(inp,'nslice',getattr(A_input,'nslice'))
        else:
            beamf = getattr(inp,'beam')
            if hasattr(beamf,'I'):
                setattr(inp,'curlen', getattr(beamf,tpulse) * speed_of_light / 1e15)
                setattr(inp,'nslice',8 * int(inp.curlen / inp.zsep / inp.xlamds))
            else:
                setattr(inp,'nslice',getattr(A_input,'nslice'))

        if (getattr(self,'i_edist')==1) or (getattr(inp,'edist')!=None) or  (getattr(inp,'beam')!=None) :
            setattr(inp,'ntail',0)
        else:
            if (getattr(self,'i_edist')==0) and getattr(A_input,'ntail')!=0 :
                setattr(inp,'ntail',int(getattr(A_input,'ntail')))
            else:
                setattr(inp,'ntail',-int(np.floor(getattr(inp,'nslice')/2)))
        # Overwrite the simulation attributes if the user has new values for them defined in the input data structure
        if (hasattr(self, 'i_rewrite')) and (hasattr(self, 'par_rew')) and (getattr(self, 'i_rewrite') == 1):
            inp = self.GEN_rewrite_par(inp)
        else:
            pass

        # Running over noise realisations and/or scan parameters
        for n_par in s_scan:
            for run_id in run_ids:           
                inp.runid = run_id
                #inp.lout =  [1,1,1,1,1,0,1,1,1,1,1,0,0,1,0,0,0,0,0]
                if ((self.stat_run==1) or (self.i_scan==1)):
                    setattr(inp,'ipseed',-1)
                else:
                    ipseed = np.random.randint(9999)
                    setattr(inp,'ipseed', ipseed)
                inp.run_dir = getattr(self,'file_pout')+'scan_'+str(n_par)+'/ip_seed_'+str(inp.ipseed)+'/' 
                
                try:
                    os.makedirs(inp.run_dir)
                    if self.stat_run>1:
                        while ipseed in  [int(files[8:]) for files in os.listdir(getattr(self,'file_pout')+'scan_'+str(n_par)) if os.path.isdir(getattr(self,'file_pout')+'scan_'+str(n_par)+'/'+files)]:
                            ipseed = np.random.randint(9999)
                            shutil.rmtree(inp.run_dir)
                        else:
                            setattr(inp,'ipseed',ipseed)
                            inp.run_dir = getattr(self,'file_pout')+'scan_'+str(n_par)+'/ip_seed_'+str(inp.ipseed)+'/'
                            os.makedirs(inp.run_dir)
                except OSError as exc:
                    if (exc.errno == errno.EEXIST) and os.path.isdir(self.file_pout+'scan'+str(n_par)+'/run_'+str(run_id)):
                        pass
                    else: 
                        raise
                if (A_input.latticefile != None or A_input.latticefile != 0) and A_input.magin ==1:
                    shutil.copyfile(os.path.dirname(self.gen_file)+'/'+A_input.latticefile,inp.run_dir+A_input.latticefile) 
                 
                if self.i_scan==1 and A_input.latticefile==None:
                    inp= self.GEN_scan(n_par ,A_input,A_undulator,inp)                      
                elif self.i_scan==0 and inp.f1st==1 and A_input.latticefile ==None:
                    inp.lat = A_undulator['Magnetic Lattice']
                    setattr(inp,'magin',1)
                elif A_input.latticefile !=None:
                    setattr(inp,'magin',1)
                    setattr(inp,'latticefile',getattr(A_input,'latticefile'))
                elif self.i_scan!=0 and inp.f1st ==0 and A_input.latticefile ==None:
                    inp.lat =None
                    setattr(inp,'magin',0)
                if inp.wcoefz[2]!= 0 :
                    inp.lat =None
                    inp.magin =0
                inp_arr.append(inp)                   
                launcher=get_genesis_launcher(self.gen_launch)
                print('+++++ Starting simulation of noise realisation {0}'.format(run_id))
                g = run_genesis(inp,launcher,i_aft=i_aft)
                setattr(g,'filePath',str(inp.run_dir))
                if (inp.itdp==1):
                    plot_gen_out_all(handle=g, savefig=True, showfig=False, choice=(1, 1, 1, 1, 3.05, 0, 0, 0, 0, 0, 0, 0, 0), vartype_dfl=complex128, debug=1)
                inp.latticefile=None
                inp.outputfile=None
                inp.edistfile=None
                inp.beamfile=None
                inp.fieldfile=None
                inp.radfile=None
        print('++++ postprocessing (results saved in {0}results++++'.format(self.file_pout))
        if i_aft!='xls':
            self.post_processing(inp,s_scan)
        plt.close("all")
        return inp_arr
