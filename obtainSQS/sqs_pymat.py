from pymatgen.transformations.advanced_transformations import SQSTransformation
from pymatgen.core.structure import Structure

sqs_transform = SQSTransformation(scaling = 4)
s = Structure.from_file("/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/Outputs_HCP/ele_POSCARs/Cr_HCP.vasp")
for i in s.sites:
	i.species = {"Fe":0.5, "Cr":0.5}
print(s)
print(sqs_transform.apply_transformation(s))