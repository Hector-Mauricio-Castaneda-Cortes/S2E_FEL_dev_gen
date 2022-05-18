#!/usr/bin/env python
# coding: utf-8

# # DEMO OF FEL SIMULATION BLOCK AND FEL CALCULATIONS
# 
# #### Hector Mauricio Castaneda Cortes
# #### 17/05/2022
# 
# 
# 
# ## Table of contents
# 1. [Setting up packages and installation](#introduction)
# 2. [Import the FEL simulation block class](#FEL_simulation_block)
#     1. [Setting up the data dictionary](#data_dict)
#     2. [Make an instance of the FEL simulation block](#instance)
#     3. [Read Genesis input file](#read_genesis_file)
#     4. [FEL parameters calculated by OCELOT](#fel_parameters)
#     5. [Run simulation](#run_simulation)
#     6. [Run several simulations changing parameters - TD](#loop_over_parameter)
#     
# 3. [Reading output file](#reading_output_file)
# 4. [Brightness and BW calculation across the undulator](#brightness_bw)
# 5. [Calculation of the coherence length](#coherence_length)

# ## Setting up the packages that are required for the FEL simulation block (Lynx)
#  <a name="introduction"></a>
# #### Clone the repository from [Lynx repository](https://github.com/Hector-Mauricio-Castaneda-Cortes/S2E_FEL_dev_gen)
# 
# 
# Firstly, the path where OCELOT is installed needs to be appended. That is why os and sys are imported. If you are running on apsv2 or apclara1, you can also export the PATH (EXPORT PATH=$PATH:path_ocelot). In addition, we are importing numpy, matplotlib for plotting purposes and deepcopy from copy

# In[1]:


import os
import sys
sys.path.append("C:////Users/qfi29231/Documents/S2E_FEL_dev_gen/ocelot/")
import numpy as np
import matplotlib.pyplot as plt
from copy import copy, deepcopy


# In addition, we are going to import the methods within the ocelot.adaptor.genesis (including the input, output, edist, dpa, dfl classes and methods to read the output files, dpa, dfl files, convert edist to beam or ASTRA to edist files)
# 
# 

# In[19]:


from ocelot.adaptors.genesis import * ##Not very pythonic, but it will import all the methods and classes within the ocelot.adaptors.genesis package
from ocelot.rad.fel import *


# ## FEL simulation block <a name="FEL_simulation_block"></a>
# 
# ### Import the Class
# 
# To import the FEL simulation block class (and run the simulations), it is necessary to import the class from the ocelot.S2E_STFC package. THe implementation goes as follows:

# In[20]:


from ocelot.S2E_STFC.FEL_simulation_block import FEL_simulation_block


# ### Setting up the data dictionary <a name="data_dict"></a>
# 
# The data dictionary is the input that the initialise the instance of the FEL simulation block object(let's call it Genesis simulation object). It will initialise some of the attributes that the instance will have once called.
# 
# - gen_file: Path of the Genesis input file (v2 or v3)
# - file_pout: Path where the input files and output files are generated. If there is only one statistical run, the tree directory will be file_pout/scan_0/ip_seed_-1/  (by default, the ipseed parameter for one noise realisation is -1).
# - stat_run: Number of statistical runs.
# - gen_launch: Genesis executable. It could be genesis2.0 (in apsv2) or genesis_mpi (in apclara1). I copied the genesis_mpi executable to my home directory in apclara1, and modified the code such that when the simulation launches, the executable in my home directory can be called (if you need to modify the number of processes or t he executable, you need to look for the method get_genesis_launcher in ocelot/adaptors/genesis.py).
# - i_edist, 'file_edist': Flag and path where an external edist is located. It uses the function read_edist from ocelot/adaptors/genesis.py If the flag is 0 (default), it does not look for an edist file.
# - i_beam, 'file_file_beam': Flag and path where an external beam is located. It uses the function read_beam from ocelot/adaptors/genesis.py (the size of the bin is 100 wavelengths by default. You can modify it in the script ocelot/S2E_STFC/FEL_simulation_block.py, GEN_read_beam() method).  If the flag is 0 (default), it does not look for a beam file.
# - i_rad, 'file_rad': Flag and path where an external rad is located.   If the flag is 0 (default), it does not look for an rad file.
# - i_dpa, 'file_dpa': Flag and path where an external dpa file is located. It uses the function read_dpa from ocelot/adaptors/genesis.py (you can set the number of particles in the dpa file). You can modify it in the script ocelot/S2E_STFC/FEL_simulation_block.py).  If the flag is 0 (default), it does not look for a dpa file.
# - i_dfl, 'file_dfl': Flag and path where an external dfl file is located. It uses the function read_dfl from ocelot/adaptors/genesis.py.  If the flag is 0 (default), it does not look for a dfl file.
# - i_astra, 'astra_file': Flag and  path where an external ASTRA file is located. It converts the ASTRA file into an edist object. If the flag is 0 (default), it does not look for an ASTRA file.
# - i_match, tw_match dictionary: If the flag is set to 1, it will match the external distribution to the twiss parameters which are defined within the tw_match dictionary (keys: alpha_x,alpha_y,beta_x,beta_y)
# - i_rewrite, par_rew: If the flat is set to 1, it will rewrite any parameter of the original input file, replacing it by a new value defined within the dictionary par_rew.
# - i_scan: Flag to perform a scan over a parameter in the input file (if it is set to one):
#     * Parameter: Parameter of the input file to be scanned (if you do a scan, the tree directory will look like file_pout/scan_value-parameter/ip_seed_-1/). The scan is done in steady state by definition. <ins>If you want to do a time dependent scan over a variable, just call the simulation within a loop and replace the parameter in the input object.</ins>
#     * Init: starting value of the parameter
#     * End: Ending value of the parameter for the scan
#     * n_scan: Number of scanned values
# 
# 

# In[21]:


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


# ### Making an instance of the FEL_simulation_block class <a name="instance"></a>
# 
# Once you have defined the data dictionary, you can make the instance of the FEL simulation block (simulation run object) with the attributes as defined within the input data dictionary

# In[22]:


f=FEL_simulation_block(data)


# ### Read GENESIS input file (calling the method read_GEN_input_file within the class) <a name="read_genesis_file"></a>
# 
# The method read_GEN_input_file() is used to read the input file and generate the input object. THe attributes of this object can be overwritten (if you need to change a parameter of the simulation in the input file, just type A_inp.aw0=0.1, and it will change the aw0 to 0.1 in the input object). Just in case, to avoid errors, once you have created the instance of the input object, change fbess0, trama, to zero. Also, the input object has as attributes dpa, dfl, beam and edist objects. The idea is that you can use the methods in the adaptors/genesis.py package to read external dpa, dfl, beam and edist files. And these objects will become dpa, dfl,beam,edist files generated before the simulation starts.

# In[23]:


A_inp=f.read_GEN_input_file()
A_inp.fbess0=int(0) # Change the parameter fbess0 of the input object
A_inp.trama=int(0)
## A_inp.edist= read_edist_file(path_to_edist_file)  If you want to read an external edist file (you need to set A_input.distfile=None)
## A_inp.beam= read_beam_file(path_to_beam_file,step=step in m)  If you want to read an external beam file (you need to set A_input.beamfile=None)
## A_inp.dfl= read_dfl_file(path_to_dfl_file)  If you want to read an external dfl file (you need to set A_input.fieldfile=None)
## A_inp.dpa = read_dpa_file(path_to_dpa_file)  If you want to read an external beam file (you need to set A_input.partfile=None)

# The dfl, dpa, beam, and edist objects have attributes beam.filePath. If you define those attributes and plot the beam or the edist object, you can
#select where the results (plots) are to be saved.

# It fails because the input file is not in the gen_file path. Change it in the data dictionary
# For a steady state simulation
#           A_inp.type = 'steady'
#           A_inp.shotnoise=0
#           A_inp.prad0=1e3
#           A_inp.itdp=0

# For a TD
#           A_inp.type='tdp'
#           A_inp.shotnoise=1
#           A_inp.prad0=1e-12
#           A_inp.itdp=1


# ### FEL parameters <a name="fel_parameters"></a>
# 
# There is a function within the ocelot package that allows you to display some FEL properties calculated from the input file. It is a method called CalculateFelParameters. It gets the input object as input of the function. You can print out the FEL parameters calculated from the input file , and then, once you create the instance by calling the calculate FEL Parameters, you can access the value of the parameter listed

# In[32]:


printFelParameters(A_inp)

fel= calculateFelParameters(A_inp) # it creates an instance of the FEL Parameters class
print('rho={},etad={},lscale={}'.format(fel.rho,fel.xie_etad,fel.xie_lscale))


# ### Calling the method to run the simulation and making the postprocessing (GEN_simul_preproc). Input: the Input object which is the outcome of the read_GEN_input_file <a name="run_simulation"></a>
# 
# Once you have set up everything related to your input object, you can start the simulation, calling the method GEN_simu_preproc from the FEL_simulation_block object. The input of that method is the input object you have defined earlier.
# 
# The output is an array of output objects (output objects are objects that are generated when you read the output file of Genesis. WHen you run the simulation, it will read the output file at the end. To access the output file, you just need to do g_out[0]. You read it and plot it as you do it in MATLAB, I think

# In[13]:


g_out=f.GEN_simul_preproc(A_inp)

#If you want to have access the genesis output object (at the end of the simulation), and you have only one noise realisation,
# you need just to type g_out= g_out[0]. If you have more than one noise realisation, all the output objects will be appended
# in an array (in this case, the name of the array is g_out). If you want to know to which simulation does a particular output
#object belongs, you can ask Python for the filepath of the output object:

# gout= g_out[0]
# print(gout.filePath)

# If you want to get the pulse energy, just type gout.energy (length of this array is the same as gout.z, 
#position along the undulator)

#If you want to get the power, gout.p_int (dimentions of the matrix, number of slices times number of positions along the undulator)

# If you want to get the spectrum, gout.spec (it has the dimentions of the length of gout.freq_eV or gout.freq_nm , frequencies
# in eV or nm and positions along the undulator, so that you can get the spectra at some particular positions of the undulator) . 
# Therefore, you can plot the spectra


# ### Loop scan over parameter <a name="loop_over_parameter"></a>
# 
# If you want to do a loop over a parameter, time dependent, you can do the following (for example scan of aw0 in time dependent)

# In[ ]:


start_point=0.
end_point=1.e-6
num = 10
array_input = np.linspace(start_point,end_point,num=num)
g_out=[]

f=FEL_simulation_block(data)
file_old = f.file_pout # directory where to save the results
for aw00 in array_input:
    f.file_pout = os.path.join(file_old,'aw0_'+str(aw0)) # generates subfolders with the results of each step of the scan
    A_inp=f.read_GEN_input_file()
    A_inp.fbess0=int(0) # Change the parameter fbess0 of the input object
    A_inp.trama=int(0)
    A_inp.tdp='tdp'
    A_inp.shotnoise=1
    A_inp.prad0=1e-12
    A_inp.itdp=1
    g_out.append(f.GEN_simul_preproc(A_inp)) # it appends an output object at the end of each simulation. 
    


# ## Reading of the output file <a name="reading_output_file"></a>
# 
# If you have run a simulation and you have a genesisv2 file, you can read it using the method read_file_out from adaptors/genesis.py. The input will be  the path of the output file (change the path of the output file in the next cell to an available output file).

# In[14]:


g_out = read_out_file('/scratch2b/qfi29231/XLS/run.0.gout')


# ## Calculation of Brightness, BW <a name="brightness_bw"></a>
# 
# If you have the output object, you can calculate the brightness across the undulator and the BW. These methods are defined within the FEL_simulation_block object (simulation object). So, you can define a dummy input data dictionary (in case you are not running a previous simulation) and the method will receive an output object as input

# In[16]:


data_dummy={} ## Dummy data dictionary, g_out output object
f_obj_two = FEL_simulation_block(data_dummy)
brightness_bw = f_obj_two.brightness_calculation(g_out)


# The brightness_bw is a dictionary and has the results of the brightness and BW calculation . It has three keys: 
# 
# - 'brightness' :Brightness across the undulator (length the same as g_out.z)
# - 'bw_std': BW across the undulator (as plotted by xgenesis)
# - 'z': Position across the undulator
# 
# if you need to call the brightness, you can access it via brightness_bw['brightness']
# 

# ## Calculation of coherence length <a name="coherence_length"></a>
# 
# Similarly to the case of the Brightness, if you have an output object (after reading an output file), you can create a dummy dictionary to have access to the method to calculate the coherence length of a pulse as outcome of a FEL simulation. The calculation of the coherence length follows a script implemented by Neil in MATLAB, in which you determine from the output the amplitude and phase of the electric field, and use the usual definition of the degree of coherence to integrate the field and its conjugate at different times and transverse positions. 

# In[ ]:


data_dummy={} ## Dummy data dictionary, g_out output object
f_obj_two = FEL_simulation_block(data_dummy)
coh_length =f_obj_two.get_lcoh(g_out)




