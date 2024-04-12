from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.structure import Structure

malnac_struct = Structure.from_file("../Outputs_HCP_old/ele_POSCARs/Cr_HCP.vasp")
malnac02_struct = Structure.from_file("../Outputs_HCP_old/ele_POSCARs/Cr_HCP_trans.vasp")
matcher = StructureMatcher(primitive_cell=False)
print(matcher.fit(malnac_struct, malnac02_struct))
print(matcher.get_transformation(malnac_struct, malnac02_struct)[0].T)
