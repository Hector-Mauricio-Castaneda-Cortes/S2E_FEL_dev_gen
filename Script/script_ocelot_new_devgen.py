#!/usr/bin/python
'''Script to run the genesis simulation in OCELOT. It calls the different
elements that are required to make instances of objects and run genesis.
 It requires to be run in the directory where the input file is (it looks for a .in file
with the specifications of lattice design, simulation specs and hard coded values of 
Beam parameters.  It is necessary to append the path where OCELOT is being installed.

HMCC: 16-02-17 Initial prototype. Including Scan over parameters (from input file)
                and noise realisations. Things still to be implemented: Postprocessing,
                matching of the beam (via betamatch or using the other methods within the 
                OCELOT framework), tapering, wakefields.
HMCC: 21-02-17: Adding the function to extract the information from an existent GENESIS input file. The GENESIS input file should be in the root directory where the script input file and the script are located.
HMCC: 03-03-17: Fixing bug in the read_GEN_input_file function
                Adding the function to call the executable of betamatch. The executable needs to be in the directory where the script is being run
HMCC: 06-03-17: Adding the betamatch functionality. When scanning over quads, the scan is performed in steady state rather than time dependent. Modifications in oder to allow the user to perform the 
                scan over quads
HMCC: 08-03-17: Adding the postprocessing function from OCELOT to plot the results of statistical runs (OCELOT modifications are required for this).
HMCC: 09-03-17: Adding the function gen_outplot_single in order to plot scan and single runs. Calling the function gen_outplot_statistics to plot statistical runs. Add support to import and read existent
                GENESIS input files by using the implemented function GEN_input_file. At the moment, for single runs, the gen_outplot_single can only plot results at the end of the undulator, taking the 
                average over noise realisations and the maximum (peak power)
HMCC: 21-03-17  Fixing post-processing and calling new functions for plotting that are part of the new implementation of OCELOT.
HMCC: 22-03-17 Fixing nslice bug after the creation of the input file. Add support to the flat-top current profile by adding a flag i_profile to the input file.


'''
#################################################
### import of all modules that are required.
from __future__ import print_function
from copy import deepcopy,copy
import sys, os

from ocelot import *
from ocelot.utils.xfel_utils import *
from ocelot.gui.accelerator import *
from ocelot.gui.genesis_plot import *
from ocelot.optics.elements import Filter_freq
from ocelot.adaptors.genesis import *
import ocelot.cpbd.elements
import numpy as np
import matplotlib.pyplot as plt


### Path where OCELOT is being installed
sys.path.append('/home/qfi29231/.local/lib/python2.7/site-packages/Ocelot-16.8rc0-py2.7.egg/ocelot/')
####  Methods #################################

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir,name))]

def read_GEN_input_file(filename):
    from ocelot.adaptors.genesis import GenesisInput
    A_input = GenesisInput()
    with open(filename, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            else:
                splitLine = line.rsplit('=') 
                if splitLine[0].startswith('alphax'):
                    splitLine[0]=splitLine[0].replace('=','').rstrip()
                else:
                    splitLine[0]=splitLine[0][1:].replace('=','').rstrip()         
                splitLine[-1] = splitLine[-1].replace('D+','e').replace('D-','e-').replace('\n','').replace('\r','').rstrip()
                splitLine[-1] = splitLine[-1].replace('E+','e').replace('E-','e-').replace('\n','').replace('\r','').rstrip()
                if(splitLine[0][:].endswith('file') and splitLine[0][:].startswith('mag')):
                    val_attr = str(splitLine[-1])
                elif (splitLine[0].startswith('$n')) or (splitLine[0].find('filetype')!=-1):
                    continue
                elif(splitLine[0].startswith('$end')):
                    break
                elif (str(splitLine[0]).startswith('i') or str(splitLine[0]).startswith('n') 
                      or str(splitLine[0]).startswith('lbc') or str(splitLine[0]).startswith('magin') or str(splitLine[0]).startswith('magout')   
                      or str(splitLine[0]).startswith('ffspec') or str(splitLine[0]).startswith('convharm') 
                      or str(splitLine[0]).startswith('multconv')):
                    val_attr = int(float(splitLine[-1].replace('=',"")))
                elif str(splitLine[0]).startswith('xk') :
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
                else:
                    val_attr = float(splitLine[-1].replace('=',""))        
                setattr(A_input,str(splitLine[0]),val_attr)                
    f.close()
    return A_input

def read_input_file(f_path):
    if (f_path.endswith('/') is not 'True'):
        f_path=f_path+'/'
    else:
        f_path=f_path
    A_dir = [files for files in os.listdir(f_path) if files.endswith('input_file.in')]
    A_arg_int = ['stat_run',',ipseed','nsec','nwig','zsep','fl','dl','drl',
                 'npart','nslice','shotnoise','delz','dmpfld',
                 'itdp', 'delz','ncar','iwityp','ntail','fbess0','gauss_fl',
                 'rmax0','i_scan','n_scan','in_gen']
    A_arg_float = ['xlamd','xlamds','gamma0','aw0','awd','quadf','quadd',
                   'curlen','curpeak','delgam','rxbeam','rybeam',
                   'emitx','emity','alphax','alphay','prad0','dgrid',
                   'xkx','xky','init','end']
    A_content={}
    with open(f_path+str(A_dir[0]),'r') as f:
        for indx,line in enumerate(f.readlines()):
            splitLine = (line.replace('\n','')).rsplit()
            if splitLine[0].startswith('#'):
                continue
            elif splitLine[0] in A_arg_int:
                splitLine[-1]=int(splitLine[-1])
            elif splitLine[0] in A_arg_float:
                splitLine[-1] = float(splitLine[-1])
            elif splitLine[0].startswith('exp_dir'):
                splitLine[-1]=splitLine[-1][splitLine[-1].find('/')-1:-1]
            A_content[splitLine[0]]=splitLine[-1]
            if splitLine[0].startswith('parameter'):
                A_content['parameter']=str(splitLine[-1][1:])
            elif splitLine[0].startswith('genesis_launcher'):
                A_content['genesis_launcher']= str(splitLine[-1]) 
            elif splitLine[0].startswith('gen_filename'):
                A_content['gen_filename']= str(splitLine[-1])
    f.close() 
    if A_content['in_gen']==0:
        print('++++++++++ No GENESIS input file ++++++++++++')
        pass
    else:
        print('++++++++++ GENESIS input file {} ++++++++++++'.format(A_content['gen_filename']))
        A_inp = read_GEN_input_file(A_content['gen_filename'])
        for attr in A_inp.__dict__.keys():
            if (attr in A_arg_int) or (attr in A_arg_float): 
                A_content[attr]=getattr(A_inp,attr)
            else:
                continue
    return A_content
    
def undulator_design(A_contents):
    from ocelot.cpbd.elements import Drift, Quadrupole, Undulator
    from ocelot.cpbd.magnetic_lattice import MagneticLattice
    from ocelot.common import globals
    from ocelot.common.globals import m_e_GeV, speed_of_light, h_eV_s
    from ocelot.rad.undulator_params import Ephoton2K
    from ocelot.rad.undulator_params import UndulatorParameters
    import numpy as np

    ## Taking it from the Notebook
    xlamd = A_contents['xlamd']
    nwig=A_contents['nwig']
    E_beam=A_contents['gamma0']*m_e_GeV
    E_photon=h_eV_s*speed_of_light/A_contents['xlamds']
    p_beam = np.sqrt(E_beam**2-m_e_GeV**2)
    fl=A_contents['fl']
    dl=A_contents['dl']
    drl=A_contents['drl']
    quadf=A_contents['quadf']
    quadd=A_contents['quadd']
    nsec=A_contents['nsec']
    drl = int((fl+drl-nwig)/2)-1
    if A_contents['iwityp']==0:
        kw0=float(A_contents['aw0'])
    elif A_contents['iwityp']==1:
        kw0=float(A_contents['aw0']*sqrt(2))
   
   # Instance of the Undulator object
    und= Undulator(lperiod=xlamd, nperiods=nwig, Kx=kw0)
   
   # Calculation of the Undulator parameter from the Photon and Beam Energies)
    if A_contents['i_scan']==1 and A_contents['parameter']=='aw0':
        und.Kx = kw0
    else:
        und.Kx=Ephoton2K(E_photon,und.lperiod,E_beam)    
   
   # Drift sections (if they are the same)
    d_rift = Drift(l=drl*und.lperiod)
    d_rift2 = Drift(l=drl*und.lperiod)
   
   # Definition of Quads
 
    qf= Quadrupole(l=fl*und.lperiod,k1=0.3*quadf/p_beam)
    qd = Quadrupole(l=dl*und.lperiod,k1=0.3*quadd/p_beam)
    qdh=deepcopy(qd)
    qdh.l/=2
   
   # Creating the cell
 
    extra_fodo = (und,d_rift,qdh)
    cell_ps = (und, d_rift,qf,d_rift2, und, d_rift, qd, d_rift2)
    l_fodo = MagneticLattice(cell_ps).totalLen/2 ##Length of fodo cell
    sase3= MagneticLattice((und,d_rift, qd,d_rift2)+int(nsec/2)*cell_ps) # Lattice with nsec modules
    up = UndulatorParameters(und,E_beam) # Instance of the Class UndulatorParameters
    print('++++ Undulator Parameters +++')
    up.printParameters()
    return {'Undulator Parameters':up,'Magnetic Lattice':sase3}
    
def BeamDefinition(A_contents):
    from ocelot.common.globals import m_e_GeV, speed_of_light
    from ocelot.cpbd.beam import Beam   
    beamf = Beam()
    A_dict= {'E':A_contents['gamma0']*m_e_GeV,'sigma_E':A_contents['delgam']*m_e_GeV,'beta_x':A_contents['gamma0']*(A_contents['rxbeam']**2)/(A_contents['emitx']), 
        'beta_y':A_contents['gamma0']*(A_contents['rybeam']**2)/A_contents['emity'], 'alpha_x':A_contents['alphax'],'alpha_y':A_contents['alphay'],
        'emit_x':A_contents['emitx']/A_contents['gamma0'],'emit_y' : A_contents['emity']/A_contents['gamma0'],'emit_xn':A_contents['emitx'],'emit_yn':A_contents['emity'],
         'x' :  0.000000e+00,
        'y' : 0.000000e+00,'px':0,'py':0,'I':A_contents['curpeak'],'tpulse':1e15*(A_contents['curlen']/speed_of_light)}
    for item in A_dict:
        setattr(beamf, item,A_dict[item])
    return beamf

def beta_matching(inp,f_path):
    import os
    import shutil
    from ocelot.adaptors.genesis import filename_from_path
    from ocelot.common.globals import m_e_GeV

    A_params = ['rxbeam', 'rybeam','alphax','alphay','emitx','emity']
    inp0 = inp
    os.system('cp /scratch2b/qfi29231/betamatch_dir/betamatch %s' %(f_path+'/betamatch'))
    os.chdir(f_path)
    with open(f_path+'/mod_file.in','w') as f:
        f.write(inp0.input())
    f.close()
    with open(f_path+'/beta_input_file.in','w') as f:
        f.write('mod_file.in')
    f.close()
    with open(f_path+'/'+inp0.latticefile,'w') as f:
        f.write(generate_lattice(inp0.lat,unit=inp0.xlamd*inp0.delz,energy=inp0.gamma0*m_e_GeV))   
    f.close()
    try:
        os.system('./betamatch < beta_input_file.in')
    except:
        print('Betamatch did not work')
        
    inp2 = read_GEN_input_file(f_path+'/TEMPLATE.IN')
    for params in A_params:
        if getattr(inp0,params)==getattr(inp2,params):
            print('Betamatch did not work')
        else:
            print('Betamatch worked')
            setattr(inp0,params,getattr(inp2,params))
    os.remove(f_path+'/beta_input_file.in')
    os.remove(f_path+'/TEMPLATE.IN')
    os.remove(f_path+'/mod_file.in')
    return inp0

def gen_outplot_single(proj_dir,run_inp= [], itdp = True,savefig=True):
    import matplotlib.pyplot as plt
    from ocelot.gui.genesis_plot import subfig_rad_pow, subfig_rad_pow_evo
    from ocelot.adaptors.genesis import read_out_file
    import copy, os
    from numpy import swapaxes,mean

    dict_name={'p_int':'radiation power','energy': 'radiation pulse energy','el_e_spread': 'el.beam energy spread','el_energy': 'el.beam energy average','bunching': 'el.beam bunching','spec': 'radiation on-axis spectral density','dfl_spec':'total radiation spectral density','r_size':'radiation transv size','r_size_weighted':'radiation transv size (weighted)','xrms':'el.beam x size','yrms':'el.beam y size','error':'genesis simulation error','p_mid':'radiation power on-axis','phi_mid':'radiation phase on-axis','increment':'radiation power increment'}
    dict_unit={'p_int':'[W]','energy': '[J]','el_e_spread': '(gamma)','el_energy': '(gamma)','bunching': '','spec': '[arb.units]','dfl_spec': '[arb.units]','r_size':'[m]','xrms':'[m]','yrms':'[m]','error':''}
    figsize = (14,7)
    
    if proj_dir[-1]!='/':
        proj_dir+='/'
    if run_inp==[]:
        run_range=xrange(1000)
    #elif run_inp ==1:
     #   run_range=xrange(1)
    else:
        run_range=run_inp
    if itdp ==True:
        param_inp = ['energy','spec','bunching','el_e_spread','xrms','yrms']
        s_inp=['max','mean']
        z_inp=[0,'end'] 
        n_seeds = len([proj_dir+'scan_'+str(run_range[0])+'/'+files for files in os.listdir(proj_dir+'scan_'+str(run_range[0])) if files.startswith('ip_s')])
        n_totl=int(n_seeds)*int(len(run_range))
    else:
        param_inp=['el_e_spread','bunching','el_e_spread','xrms','yrms']
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

    n_seeds = len([proj_dir+'scan_'+str(run_range[0])+'/'+files for files in os.listdir(proj_dir+'scan_'+str(run_range[0])) if files.startswith('ip_s')])
    colour=plt.cm.Set1(np.linspace(0.1,0.999999,n_totl))
    colour = map(lambda rgb:'#%02x%02x%02x' %(rgb[0]*255,rgb[1]*255,rgb[2]*255),
                 tuple(colour[:,0:-1]))
    fig,ax = plt.subplots()
    for i_run,runr in enumerate(run_range):
        ip_s=[]
        ip_seed = [proj_dir+'scan_'+str(runr)+'/'+files for files in os.listdir(proj_dir+'scan_'+str(runr)) if files.startswith('ip_s')]
        colourr = colour[i_run]
        for i_seed, ipseed in enumerate(ip_seed):
            out_file = [ipseed+'/'+files for files  in os.listdir(ipseed) if ((files.startswith('run.')) and (files.endswith('.gout')))]
            if os.path.isfile(out_file[0]):
                g=read_out_file(str(out_file[0]),read_level=2)
                ip_s.append(g)
                
            text_leg = 'Scan ='+str(runr)+' ip_seed = '+str(ipseed[ipseed.find('ip_seed_')+8:])
            text_l.append(text_leg)
           # colourr = colour[i_seed]
            fig=subfig_rad_pow(ax,g,text_leg,colour[i_run],log=1)
            fig = subfig_rad_pow_evo(ax,g,text_leg,norm=1)
        run_range_good.append(i_run)
        outlist.append(ip_s)
    run_range=run_range_good
    if savefig==True:
        savefig='png'
        saving_path=proj_dir+'results/'
        if not os.path.isdir(saving_path):
            os.makedirs(saving_path)
        print('      saving to '+saving_path)        
    art = []
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize('8')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize('8')
    ax.grid(True)
    lgd = plt.legend(loc=9,bbox_to_anchor=(1.01,0.9),prop={'size':4.5})
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
                        param_matrix=copy.deepcopy(getattr(outlist[irun][i_seed],param))

                    if len(param_matrix) == len(outlist[irun][i_seed].z):
                        s_value=param_matrix
                    else:
                        if s_ind=='max':
                            s_value=np.amax(param_matrix,axis=0)
                        elif s_ind=='max_cur':
                            s_value=param_matrix[outlist[irun][i_seed].sn_Imax,:]
                        elif s_ind=='mean':
                            s_value=np.mean(param_matrix,axis=0)
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
                        param_matrix=copy.deepcopy(getattr(outlist[irun][i_seed],param))#HMCC

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
               
                    

##########  MAIN ###############################################      
def main():
    import os,sys
    import numpy as np
    from ocelot.common import globals
    from ocelot.adaptors.genesis import generate_input,generate_lattice,get_genesis_launcher,run_genesis,filename_from_path
    from ocelot.gui.genesis_plot import plot_gen_out_all,plot_gen_stat
   
    # Reading the input file

    A_input = read_input_file(os.getcwd()+'/')
    A_bbeam = ['gamma0','curlen','curpeak','delgam','rxbeam','rybeam','emitx','emity','alphax','alphay','xbeam','ybeam','pxbeam','pybeam']
    A_simul = ['npart','ncar','zsep','delz','dmpfld','fbess0','dgrid','rmax0','xkx','xky','iwityp']
    A_td = ['itdp','prad0','shotnoise']
    A_und = ['quadd', 'quadf','fl','dl','drl','nsec','nwig','aw0']

    #set simulation parameters by calling the read_input_file function
    exp_dir=A_input['exp_dir'][1:-1]+'/'
    print('++++ Output Path {} ++++++'.format(exp_dir))

    # Setting the number of noise realisations and scan (over quads or different parameters)
    A_input['parameter']= str(A_input['parameter'][:-1])
    if (A_input['i_scan'] ==0):
        s_scan = range(1)  
        num = A_input['stat_run']
        run_ids = xrange(0,num)
        print('++++++++ No scan ++++++++++')
    elif (A_input['i_scan']!=0):
        if (A_input['parameter'] in A_und):
            run_ids= range(1)
            if A_input['parameter'] !='aw0':
                s_scan = range(int(A_input['init']),int(A_input['end']),int(np.ceil((A_input['end']-A_input['init'])/(A_input['n_scan']-1))))
            else:
                s_scan = np.linspace(A_input['init'],A_input['end'],A_input['n_scan']) 
            print('++++ Quad scan, parameter  {} ++++++'.format(A_input['parameter']))
        elif (A_input['parameter']=='xlamds'):
            run_ids= range(1)
            s_scan = np.linspace(A_input['init'],A_input['end'],A_input['n_scan'])
            print('++++ Quad scan, parameter  {} ++++++'.format(A_input['parameter']))
        else:        
            s_scan = np.linspace(A_input['init'],A_input['end'],A_input['n_scan'])
            num = A_input['stat_run']
            run_ids = xrange(0,num)
            print('++++ Number of noise realisations {} ++++++'.format(num))  

        # setting the undulator design( Magnetic Lattice)
    A_undulator = undulator_design(A_input)
    
        # Fill in the beam object
    A_beam = BeamDefinition(A_input)
    if (A_input['itdp']==0):
        print('++++ Steady State run +++++')
        i_tdp = False
    elif (A_input['itdp']==1):
        i_tdp = True
     
        # Generate input object
    inp = generate_input(A_undulator['Undulator Parameters'],A_beam,itdp=i_tdp)

      # Overwrite the simulation attributes of the input object with the ones defined in the input file
    for in_put in A_input:
        if (in_put in A_simul) or (in_put in A_und) or (in_put in A_td):
            setattr(inp,in_put, A_input[in_put])
        else:
            continue

        # Set the object for a flat-top or Gaussian Current distribution
    if int(A_input['i_profile'])==1:
        setattr(inp,'nslice',8*int(inp.curlen/inp.zsep/inp.xlamds))
        setattr(inp,'ntail',-getattr(inp,'nslice')/2)
    elif int(A_input['i_profile'])==0:
        setattr(inp,'nslice',int(inp.curlen/inp.zsep/inp.xlamds))
        setattr(inp,'curlen',-getattr(inp,'curlen'))
        setattr(inp,'ntail',0)
    else:
        print('++++ Profile not set to be square (i_profile = 0) or Gaussian (i_profile=1) +++++++++++')
        return

    # Running over noise realisations and/or scan parameters
    for n_par in s_scan:
        for run_id in run_ids:           
            inp.runid = run_id
            #inp.ipseed = 17111*(int(inp.runid)+1)
            inp.ipseed = np.random.randint(9999)
            inp.run_dir = exp_dir+'scan_'+str(n_par)+'/ip_seed_'+str(inp.ipseed)+'/' 
            inp_path = inp.run_dir + 'run.' + str(inp.runid) + '.inp'
            inp_file = filename_from_path(inp_path)
            try:
                os.makedirs(inp.run_dir)
            except OSError as exc:
                if (exc.errno == errno.EEXIST) and os.path.isdir(exp_dir+'scan'+str(n_par)+'/run_'+str(run_id)):
                    pass
                else: 
                    raise

            if A_input['i_scan']==1:
                if (A_input['parameter'] in A_simul):
                    setattr(inp,A_input['parameter'][1:-1],n_par)
                    inp.lat = A_undulator['Magnetic Lattice']
                    print(' ++++++++++ Scan {} of the parameter {}'.format(n_par, A_input['parameter']))
                elif ((A_input['parameter'] in A_und) or (A_input['parameter'] == 'xlamds')): 
                    print(' ++++++++++ Steady State Scan {} of the parameter {} Quad optimisation'.format(n_par, A_input['parameter']))
                    setattr(inp,'type','steady')
                    setattr(inp,'itdp',0)
                    setattr(inp,'shotnoise',0)
                    setattr(inp,'prad0',10)
                    setattr(inp,'betamatch',True)
                    if (A_input['parameter'] in A_und):
                        if A_input['parameter'] == 'aw0':
                            A_input[str(A_input['parameter'])]=n_par                            
                        else:
                            n_par = int(n_par)
                            if A_input['parameter'] in A_und[0:2]:
                                for j_quad in A_und[0:2]:
                                    if j_quad =='quadd':
                                        n_par = -n_par
                                    else:
                                        n_par = n_par
                                    A_input[str(j_quad)]=n_par                  
                            else:
                                A_input[str(A_input['parameter'])]=n_par
                    else:
                        n_par =float(n_par)
                        A_input[str(A_input['parameter'])]=n_par
                    A_undulator=undulator_design(A_input)
                    inp.lat =A_undulator['Magnetic Lattice'] 
                    inp.latticefile = inp_file+'.lat'
                    inp = beta_matching(inp,inp.run_dir) 
                    inp.latticefile=None                       
                else:
                    n_par = float(n_par) 
                    A_input[str(A_input['parameter'])]=n_par
                    setattr(inp,str(A_input['parameter']),n_par)
                    inp.lat = A_undulator['Magnetic Lattice']                        
            else:
                inp.lat = A_undulator['Magnetic Lattice']                   
               
            launcher = get_genesis_launcher(A_input['genesis_launcher'][1:-1])
            print('+++++ Starting simulation of noise realisation {}'.format(run_id))
            g = run_genesis(inp,launcher)
            setattr(g,'filePath',str(inp.run_dir))
            if (inp.itdp==1):
                fig = plot_gen_out_all(handle=g, savefig=True, showfig=False, choice=(1, 1, 1, 1, 6.05, 0, 0, 0, 0, 0, 0, 0, 0), vartype_dfl=complex128, debug=1)
            inp.latticefile=None
            inp.outputfile=None

#Post- processing
# Setting off LaTeX in order to avoid the error when a plot is generated.
    print('+++++ Post-processing +++++++++++++++++++')
    plt.rcParams['text.usetex']=False
    plt.rcParams['text.latex.unicode']=False
    if (A_input['i_scan']==0):
        if (A_input['stat_run'] >1):
            print('++++++++Statistical runs (results saved in {}results) ++++++'.format(exp_dir))
            fig2=plot_gen_stat(exp_dir, run_inp=s_scan, stage_inp=xrange(1), param_inp=[], s_param_inp=['p_int', 'energy','xrms','yrms','bunching','r_size_weighted','spec', 'error'], z_param_inp=['p_int', 'phi_mid_disp','xrms','#yrms', 'spec', 'bunching', 'wigner'], dfl_param_inp=[], run_param_inp=[], s_inp=['max','mean'], z_inp=[0,'end'], run_s_inp=[], run_z_inp=[], savefig=True, saveval=False,debug=0)
        elif (A_input['stat_run'] ==1):
            print('++++++++ Single run (results saved in {}results) ++++++'.format(exp_dir))
            if (A_input['itdp']==0):
                i_tdp = False
            else:
                i_tdp = True
            fig=gen_outplot_single(exp_dir,run_inp= s_scan, itdp = i_tdp,savefig=True)
    else:
        print('++++++++scan (results saved in {}results) ++++++'.format(exp_dir))
        fig=gen_outplot_single(exp_dir,run_inp= s_scan, itdp = False,savefig=True)        
    
                
################## END OF METHODS ############################################
if __name__=="__main__":
    print('++++++ Script to run GENESIS ++++++++')
    main()
    print('+++++++++ End of simulation +++++++++')
