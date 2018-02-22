'''S2E_afterburner class
HMCC: 14-02-18 Implementation of the daughter class S2E_afterburner(creation of class with corresponding methods)
HMCC: 22-02-18 Correction to the modulator_amplifier method in order to do betamatching of twiss parameters,
               before the amplification stage starts.
'''
#################################################

from ocelot.adaptors.genesis import *
from ocelot.gui.genesis_plot import *
from ocelot.S2E_STFC.FEL_simulation_block import *
from ocelot.S2E_STFC.ming_xie import *
from ocelot.cpbd.beam import Twiss
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
from shutil import *
from copy import copy,deepcopy


class s2e_afterburner(FEL_simulation_block):
    __pathgf = '/home/qfi29231/riley_S2E/'
    def __init__(self,*initial_data,**kwargs):
        super(s2e_afterburner,self).__init__(*initial_data,**kwargs)
        self.__seedlambda = 3e-6
        self.__outlambda = 1e-7
        self.__nslice = 3e3
        self.__npart = 2048
        self.__gamma = 490.2
        self.__modules = 10
        self.__offsetradf=-16
        self.__nslicemod = 200
        self.__zsepmod = 1
     
    @property
    def pathgf(self):
        return self.__pathgf 

    @property
    def seedlambda(self):
        return self.__seedlambda 

    @property
    def outlambda(self):
        return self.__outlambda

    @property
    def nslice(self):
        return int(self.__nslice)

    @property
    def npart(self):
        return int(self.__npart)

    @property
    def gamma(self):
        return self.__gamma
    
    @property
    def modules(self):
        return int(self.__modules)
    
    @property
    def offsetradf(self):
        return self.__offsetradf
    
    @property
    def nslicemod(self):
        return int(self.__nslicemod) 
    
    @property
    def zsepmod(self):
        return int(self.__nslicemod)
    
    @pathgf.setter
    def pathgf(self,pgf):
        self.__pathgf = pgf
 
    @seedlambda.setter
    def seedlambda(self,slambda):
        self.__seedlambda = slambda 

    @outlambda.setter
    def outlambda(self,olambda):
        self.__outlambda = olambda

    @nslicemod.setter
    def nslicemod(self,n_sl):
        self.__nslicemod = int(n_sl)

    @npart.setter
    def npart(self,n_p):
        self.__npart = n_p 

    @gamma.setter
    def gamma(self,gmma):
        self.__gamma = gmma 

    @modules.setter
    def modules(self,mod):
        self.__modules = mod
    
    @offsetradf.setter
    def offsetradf(self,off_set):
        self.__offsetradf = int(off_set)
    
    @nslice.setter
    def nslice(self,n_sl):
        self.__nslice = n_sl  
    
    @zsepmod.setter
    def zsepmod(self,zsep):
        self.__zsepmod = int(zsep)

    def search_input(self,path,ext):
        return [path+files for files in os.listdir(path) if files.endswith(ext)][0]         
     
    def generation_beam(self,inp):
        beam=GenesisBeam()
        e_s = []
        e_g = []
        e_dg =[]
        scalefactor = int(self.seedlambda/self.outlambda)
        slicespace = float(self.seedlambda*self.zsepmod)
        slicewidth = float(self.seedlambda)
        numslice = int(self.nslicemod)
        
        dpa_f = getattr(inp,'dpa')
        count1=0
    
        g_file = (dpa_f.e).reshape(self.nslicemod,self.npart)
        ph_file = (dpa_f.ph).reshape(self.nslicemod,self.npart)
        data00 = np.column_stack((ph_file,g_file))
        data00[:,0]= data00[:,0]+4
        partnum =data00.shape[0]
        hist, bin_edges = np.histogram(data00[:,0],scalefactor)
	
        for sliceindex in range (0,int(numslice)):
            for i in range(0,scalefactor-1):
                chunkdata=np.zeros(3)
                count0 = count1
                count1 = count1 + hist[i]
                tempdata=data00[count0:count1,:]
                e_s.append(sliceindex*slicespace + slicewidth*((i+1.0)/scalefactor))
                e_g.append(np.mean(tempdata[:,1]))
                e_dg.append(np.std(tempdata[:,1]))
        setattr(beam,'z',np.array(e_s))
        setattr(beam,'g0',np.array(e_g))
        setattr(beam,'dg',np.array(e_dg))
        setattr(inp,'beam',beam)
        setattr(inp,'dpa',None)
        return inp   
    
    def prep_mod_amp_aft(self, path_mod,inp):
        inp0=inp
        f_out = read_out_file(self.search_input(path_mod,'.gout'))
        dpa_f = read_dpa_file_out(f_out)
        setattr(dpa_f,'filePath',path_mod+'mod'+'_dpa') 
        setattr(inp,'dpa',dpa_f)
        edist_dpa = dpa2edist(f_out,dpa_f,num_part = int(int(self.npart)*int(f_out('nslice'))),smear=0)
        #ebeam = edist2beam(edist_dpa,step=float(self.outlambda),i_aft=1)
        #setattr(ebeam,'filePath',path_mod+'mod'+'_beam') 
        setattr(edist_dpa,'filePath',path_mod+'mod'+'_edist') 
        plot_edist(edist_dpa,figsize=20,savefig=True,showfig=False,scatter=True,plot_x_y=True,plot_xy_s=False) 
        plot_edist(edist_dpa,figsize=20,savefig=True,showfig=False,scatter=True,plot_x_y=False,plot_xy_s=True)
        plot_dpa_bucket(dpa_f, slice_num=None, repeat=1, GeV=0, figsize=4, cmap=def_cmap, scatter=False,
                        energy_mean=None, legend=True, fig_name=None, savefig=True, showfig=False, suffix='mod',
                        bins=(f_out('nbins'),f_out('nbins')), debug=1, return_mode_gamma=0)
        setattr(inp0,'beam',None)
        setattr(inp0,'edist',edist_dpa)
        setattr(inp0,'dpa',None)
        setattr(inp0,'beamfile',None)
        setattr(inp0,'partfile',None)
        setattr(inp0,'edistfile',None)
        return inp0 

    def prep_afterburner(self,path_out,input_ext):
        gen_f = [self.pathgf+files for files in os.listdir(self.pathgf) if files.endswith('.in') and files.startswith(input_ext)][0]
        setattr(self,'gen_file',gen_f) 
        inp0 = super(s2e_afterburner,self).read_GEN_input_file()
        f_out = read_out_file(self.search_input(path_out,'.gout'))
        f_out.filePath = self.search_input(path_out,'.gout')
        setattr(inp0,'dpa',read_dpa_file_out(f_out,f_out.filePath+'.dpa')) 
        setattr(inp0,'dfl',read_dfl_file_out(f_out,f_out.filePath+'.dfl'))
        setattr(inp0,'lout',[1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0])
        setattr(inp0,'outputfile',None)
        setattr(inp0,'partfile',None)
        setattr(inp0,'radfile',None)
        setattr(inp0,'fieldfile',None) 
        setattr(inp0,'beamfile',None)
        setattr(inp0,'edistfile',None) 
        return inp0

    def modulator_amplifier(self,i_f_mamp):
        f_old = getattr(self,'file_pout')
        if i_f_mamp == 'mod' or i_f_mamp == 'amp':
            setattr(self,'file_pout',f_old+i_f_mamp+'/')
            file_in = 'demo_'+i_f_mamp
            gen_f = [self.pathgf+files for files in os.listdir(self.pathgf) if files.endswith('.in') and files.startswith(file_in)][0]
            setattr(self,'gen_file',gen_f)
            f_inp = super(s2e_afterburner,self).read_GEN_input_file()
            setattr(self,'idump',1)
            if i_f_mamp =='mod': 
                setattr(f_inp,'ippart',int(1))
                setattr(f_inp,'ispart',int(1))
                setattr(f_inp,'idmppar',int(1))
                f_inp2 = super(s2e_afterburner,self).GEN_simul_preproc(f_inp) 
                f_out = read_out_file(self.search_input(f_inp2[0].run_dir,'.gout'))
                f_out.filePath = self.search_input(f_inp2[0].run_dir,'.gout')
                f_inpamp = f_inp2[0]
            elif i_f_mamp == 'amp':
                tw = Twiss()
                f_inp2 =  self.prep_mod_amp_aft(f_old+'mod/scan_0/ip_seed_-1/' ,f_inp)
                f_inp2 = super(s2e_afterburner,self).beta_matching(f_inp2,f_old)
                dict_tw = {'beta_x':getattr(f_inp2,'gamma0')*(getattr(f_inp2,'rxbeam')**2)/(getattr(f_inp2,'emitx')), 
                           'beta_y':getattr(f_inp2,'gamma0')*(getattr(f_inp2,'rybeam')**2)/getattr(f_inp2,'emity'), 
                           'alpha_x':getattr(f_inp2,'alphax'),
                           'alpha_y':getattr(f_inp2,'alphay')}                 
                for key,value in dict_tw.items():
                    setattr(tw,key,value)                
                setattr(f_inp2,'edist',rematch_edist(getattr(f_inp2,'edist'),tw))
                setattr(f_inp2,'beam',edist2beam(rematch_edist(getattr(f_inp2,'edist'),tw),step=float(self.outlambda),i_aft=1))
                setattr(f_inp2.beam,'filePath',f_old+'mod/scan_0/ip_seed_-1/'+'mod'+'_beam')
                setattr(f_inp2,'edist',None)
                f_inp0 = self.GEN_simul_preproc(f_inp2)
                f_inpamp = f_inp0[0]
                setattr(f_inpamp,'beam',None)
                setattr(f_inpamp,'dpa',None)
                setattr(f_inpamp,'beamfile',None)
                setattr(f_inpamp,'partfile',None)
                setattr(f_inpamp,'edistfile',None)
                f_oamp = read_out_file(self.search_input(f_inp0[0].run_dir,'.gout'))
                f_oamp.filePath = self.search_input(f_inp0[0].run_dir,'.gout')
                setattr(f_inpamp,'fieldfile',None)
            setattr(self,'file_pout',f_old)
            return f_inpamp                
        else:
            print(' No modulator or amplifier ')
            return

    def after_burner(self,inp_ampl):
        inp0 = inp_ampl   
        g_out=[]
        file_old = self.gen_file
        fout_old = self.file_pout
        setattr(self,'file_pout',fout_old+'afterburner/')
        f_inp = self.prep_afterburner(inp_ampl.run_dir,'Pass1_original')
        f_inp.latticefile='odd.lat' 
        setattr(self,'file_pout',fout_old+'afterburner/n_mod0/')
        g_out.append(self.run_GENESIS(f_inp)) 

        f_inp2 = self.prep_afterburner(self.file_pout+'/scan_0/ip_seed_-1/','Pass_original') 
        f_inp2.run_dir=None
        for n_mod in range(1,getattr(self,'modules')): 
            nslice = getattr(f_inp2,'nslice')
            fold = getattr(self,'file_pout')+'scan_0/ip_seed_-1/'
            setattr(self,'file_pout',fout_old+'afterburner/n_mod'+str(n_mod)+'/')  
            if  (n_mod<(self,'modules')):
                if np.remainder(n_mod,16)==0:
                    base_n = 'odd'
                else:
                    if np.remainder(n_mod,8)==0:
                        base_n = 'even'
                    else:
                        base_n = 'noquad'
            else:
                break
            f_inp2.magin = 1
            f_inp2.latticefile = base_n+'.lat'  
            setattr(f_inp2,'nslice',int(nslice+getattr(self,'offsetradf')))          
            g0 = self.run_GENESIS(f_inp2)
            g0.filePath = self.search_input(f_inp2.run_dir,'.gout')
            g_out.append(g0)
            dpa_f = read_dpa_file_out(g0,g0.filePath+'.dpa')
            dfl = read_dfl_file_out(g0,g0.filePath+'.dfl') 
            setattr(f_inp2,'dpa',dpa_f)
            setattr(f_inp2,'dfl',dfl)
            setattr(f_inp2,'beam',None)
            setattr(f_inp2,'edist',None)
            setattr(f_inp2,'beamfile',None)
            setattr(f_inp2,'edistfile',None)
            setattr(f_inp2,'outputfile',None)
            setattr(f_inp2,'radfile',None)
            setattr(f_inp2,'latticefile',None)
        return g_out

    def run_GENESIS(self,inp):
        s_scan = range(1)  
        num = self.stat_run
        run_ids = xrange(0,num)
        for n_par in s_scan:
            for run_id in run_ids: 
                inp.idump=1
                setattr(inp,'itdp',1)
                setattr(inp,'type','tdp')
                inp.runid = run_id
                inp.ipseed = -1
    
                inp.run_dir = getattr(self,'file_pout')+'scan_'+str(n_par)+'/ip_seed_'+str(inp.ipseed)+'/' 
                try:
                    os.makedirs(inp.run_dir) 
                    copyfile(getattr(self,'pathgf')+inp.latticefile,inp.run_dir+inp.latticefile)
                except OSError as exc:
                    if (exc.errno == errno.EEXIST) and os.path.isdir(self.file_pout+'scan'+str(n_par)+'/run_'+str(run_id)) and (inp.latticefile !=None):
                        pass
                    else:
                        raise
                    
                launcher=get_genesis_launcher(self.gen_launch)
                print('+++++ Starting simulation of noise realisation {0}'.format(run_id))
                g = run_genesis(inp,launcher,i_aft=1)
                setattr(g,'filePath',str(inp.run_dir))
                inp.latticefile=None
                inp.outputfile=None
                inp.edistfile=None
                inp.beamfile=None
                inp.partfile=None
                inp.radfile=None
                inp.fieldfile=None 
        super(s2e_afterburner,self).post_processing(inp,s_scan)
        plt.close("all")
        return g 
    
    def GEN_simul_preproc(self,f_inp):
        f_inp2 = super(s2e_afterburner,self).GEN_simul_preproc(f_inp,i_aft=1)
        return f_inp2
