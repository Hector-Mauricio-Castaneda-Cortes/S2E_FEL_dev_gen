# S2E_FEL_dev_gen
New version of OCELOT (including new script to run it)


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












