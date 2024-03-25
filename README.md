# mixing-enthalpy
Contributors: Joshua Cheng, Pravan Omprakash <br>
M-CUBE, WUSTL<br> <br>
Binary mixing enthalpy calculation using SQS, DFT and functional models.

To pull this repository:
```commandline
git@github.com/Pravanop/mixing-enthalpy.git
cd mixing-enthalpy
```
Create a virtual environment, and then run
```commandline
python3.9 -m venv env
source env/bin/activate
pip install -r ./requirements.txt
```
You might have to configure your interpreter settings in the IDE to point to the right env location
## Step 1

Modifying input.yaml  <br> <br>
The input.yaml file contains all the required information from the user to create binary SQS alloys and optionally 
vasp input runs (for internal use mostly).
<br><br>
The important tags to update are:
- element_list : A list of strings that contain the element space you are interested in.
- lattice : The crystal system the binaries will be created for. Options: 'BCC', 'FCC', 'HCP', 'SC'
- doping_percent: The mole fraction in the binary. Set to 50 for most purposes.
- coord_type: The creation of rndstr.in with fractional or cartesian coords seems to change the convergence of the 
  effect. Incase, the sqs doesn't retrieve fast outputs, please change the type and try again.
- supercell: Depending on the crystal system, different sizes needed to be provided. Recommended : BCC: [3, 2, 2], 
  FCC: [3, 2, 1]
- corrdump_path and mcsqs_path: This code works alongwith the ATAT code for creating SQS alloys. You can download 
  ATAT from here: . Follow the instructions to install ATAT. Then update the executable files for 'corrdump' and 
  'mcsqs' for these tags.
- run_vasp : **For Internal Group Use**  Set to true, and it will create vaspruns based on the metadata in .
  /vasp_metadata
<br> <br>

### If run_vasp is True, then please follow these instructions:

In the /vasp_metadata folder, there are two files of importance 
<br> <br>
- INCAR.json : For this project, you must update all the INCAR tags you are interested in. Please change this once 
  only and leave it for each project. 
- metadata.yaml : Update the paths in the metadata file, corresponding to base.yaml, INCAR.json and PP. Additionally 
  you might want to change the xc_functional, kpoints and runjob files. Instructions for each tag are given in the 
  metadata.yaml file
- **The PP folder for VASP, must be added into this folder, and its path updated in metadata.yaml.** 


## Step 2

Once all these input files have been modified, the DFT_input.py file can be run. An 'Outputs' folder will be created 
with all the intermediary files as well as the finished 'vasp_runs'. These vasp_run files can be sent to the cluster 
for further steps. <br>

```
python DFT_input.py
```

For further steps, the energies from VASP runs have to be added manually. **[WIP]**

## PART 2
Scripts to create phase diagrams, intermetallics, higher order enthalpy values have been created in calcEnthalpy. 
They need a mixing enthalpy dataset got from Part 1 and DFT calculations of this project.

A UI is now ready!
To run the UI, ensure that streamlit is installed in python (Check requirements.txt)
```commandline
streamlit run app.py
```
This will open up in your browser, and you can play around with it. 
This script is not bug free, especially in terms of absolute paths. It is recommended to check if any of my paths 
are laying dormant in the code, and you can change it. The UX as of itself, is quite robust with error handling 
taken care of. <br>
THE UI IS A **WIP**.