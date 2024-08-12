<h1> A Fast and Robust Method for  Predicting the Phase Stability of High Entropy Alloys using Pairwise Mixing Enthalpy</h1>


<br>
<br>

This package helps to calculate the mixing enthalpy for high entropy alloys, create phase diagrams using convex hulls to predict stability and then find favourable deposition pathways for synthesis.

To begin, 
1) Pull the repo into your local machine.
2) Create a 'Data' folder with the following structure:<br>
```
Data
│      
└───Input Data
│   └───$Source$
│       │   $Lattice$_$Source$.json
│       │   element_list_$Lattice$_$Source$.txt
│       │   ...
│   
└───Output Data
    |
```
  i) Replace $Source$ with a unique identifier of your choice. Replace $Lattice$ with 'bcc', 'hcp', 'fcc' etc. 
  ii) The .json file should have the form:
         ```
         {'El1-'El2' : mixing enthalpy value in ev/atom,
           ....}
         ```
     <br> Note: **The string **'El1-El2'** must be alphabetically sorted!**
  iii) The element text file that is of the form 'El1,El2,EL3,El4...' <br>
      **Note: There must be no space between the commas and element names. The element list does not have to contain all the elements in the json file. The computed enthalpies will be governed by the element list not the json file.**
3) Add an Api key from materials Project in the 'callMPI' folder if you want to extract intermetallics data. Currently only materials project is supported.

<br> Once the input data has been added, move into the calculateEnthalpy dataset and run calculateEnthalpyDataset.py. 
<br> Ternary Phase diagrams can be plotted using higherOrderPhaseDiagrams.
<br> Reaction pathways can be computed and plotted using reactionPathways. 

Instructions to run each module is provided in each script, with example runs.
