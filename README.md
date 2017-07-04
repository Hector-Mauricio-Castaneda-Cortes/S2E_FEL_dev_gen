# S2E_FEL_dev_gen
New version of OCELOT (including new script to run it)

 Update
 

23-03-17

After the migration of OCELOT to a new branch for development, I created a new repository which has the modifications in order to run OCELOT 
in our local server. I also have come with a new version of the Script with the pre-processing and post-processing routines 
updated accordingly.

In order to install this version of OCELOT in your apsv2 account just do:

git clone https://github.com/Hector-Mauricio-Castaneda-Cortes/S2E_FEL_dev_gen.git

and then overwrite the Ocelot folder which is located in 
/home/user_name/.local/lib/python2.7/site-packagesOcelot-16.8rc0-py2.7.egg/ocelot

I recommend you to rename this directory to 
/home/user_name/.local/lib/python2.7/site-packagesOcelot-16.8rc0-py2.7.egg/ocelot_old

In case something goes wrong. And afterwards, overwrite the folder by the one you cloned from the repository.

The repository consists of two directories:
Ocelot/
Scripts/

The Script directory has the script_ocelot_new_devgen.py on it and the input_file.in file. Those are the ones that are needed to run OCELOT. 
You can modify the input_file.in as you wish. The input_file.in needs to be in the same folder that the script_new_devgen.py

To run the script you must only type:
python script_ocelot_new_devgen.py
in the command line

If you are reading from an existent Genesis input file, you must have the GENESIS input file in the directory where the input_file.in and 
the script are located. Within the input file, there is a flag (in_gen) which is set to 0 if you do not have a GENESIS input file available 
and 1 if there is any. If so, the name of the GENESIS input file must be supplied in front of the gen_file field (the in_gen flag must 
be set to 1).

Please note that after every key within the input file there is a tab. Do not erase the tab. After the key, there must be a tab space and 
then the value of the parameter (I will fix that later).

Even if you have a genesis input file, if you want to do noise realisations or scan over quads you will need to modify the input_file.in

If you want to do noise realisations, there is a flag called stat_run. Set the number of noise realisations that are going to be run.
Here is a description of the input_file.in:


###### Output path #################################
exp_dir	      '/scratch2b/qfi29231/results_dev_single/'   (Output path where the results and output files are going to be stored)
###### Type of Launcher (Genesis v2 or v3) #########
genesis_launcher	'genesis2.0'                          (GENESIS launcher, could be genesis2.0 or genesis3.2.1)
###### Statistical runs ############################
stat_run	   1                                          (Number of statistical runs, larger than 0) 
###### GENESIS input file ##########################
in_gen 0                                                  (Flag to set if there is a GENESIS input file available to read from, 1, or not, 0).
gen_filename	gen_input.in                              (Name of the GENESIS input file. Even if there is one, noise realisations and scans are controlled via the input_file.in)  
###### Time dependent flags       ##################      (Time dependent flags (controls steady state or time dependent simulations)  
shotnoise   1
itdp	    1
prad0	    1e-12
###### Beam Energy and wavelength ##################      (Beam Energy and wavelength values) (if there is an available GENESIS input file, these values are overwritten)
xlamds	1.35076e-8
gamma0	1.46871e3
###### Undulator design and lattice ################      (Undulator parameters: They are used to fill in the attributes of the Undulator Lattice object. The lattice file is generated automatically from the object in run time)
xlamd	0.036
nsec	23
nwig	58
aw0	0.7852
fl	2
dl	2
drl	78
quadf	16
quadd	-16
iwityp	1
xkx	0.5
xky	0.5
##### Beam parameters #############################         (Beam parameters: hard coded for the moment. Eventually, the idea is to read these parameters from a dist file or a beam file and generate a beam file in run time)
i_profile  1                                                 (i_profile: Flag to determine if the beam is gaussian (1) or flat-top(0)) If it is gaussian, OCELOT calculates the time window as 8 sigma. 
curlen	2.0944e-5
curpeak	430
delgam	1.357
rxbeam	6.491585e-5
rybeam	2.95641e-5
emitx	6.5e-7
emity	6.5e-7
alphax	-2.464459e-2
alphay	-8.971086e-3
nslice	927
ntail	-463
gauss_fl	1
##### Simulation parameters #######################        (Simulation parameters: If there is an available GENESIS input file, these are overwritten)
zsep	10
npart	8192
ncar	151
delz	1
dmpfld	1
fbess0	0
dgrid	0
rmax0	9
##### Scan parameters (flag 0 if there is no scan) #        (Scan parameters)
i_scan	   0                                                (Flag for scan. If it is 0, no scan is considered. If it is 1, then the scan is performed)   
parameter  'aw0'                                             (Name of the parameter to do scan. If it is a quad scan (quadf, quadd, fl,dl,drl, nsec, nwig, aw0), the scan is done in Steady State and rematching is done with betamatch)
init	   12                                                (Initial value of the scan parameter)
end	   25                                                    (End value of the scan parameter_)
n_scan	   10                                                (Number of scans)


30-06-17

Implementation of the new class FEL_simulation_block as package of the modified version of OCELOT. Therefore, in order to run a simulation, there is no need for the dummy input , but an existent GENESIS input file which can be read by a method of the class. A script was added to the Script folder in order to show the usage of the FEL_simulation_block class

HMCC






