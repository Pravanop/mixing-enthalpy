from emmet.core.thermo import ThermoType
from pymatgen.core.composition import Composition
from calcEnthalpy_package.phase_diagram.pdEntry_local import PDEntryLocal
from pymatgen.ext.matproj import MPRester


class IntermetallicExtractions:
	"""
    Class responsible for retrieving intermetallic data.
    """
	
	@staticmethod
	def get_MP_intermetallic(alloy_list: list[str],
							 api_key: str) -> list[PDEntryLocal]:
		"""
        Retrieves intermetallic data for given alloy compositions using MPRester.

        Parameters
        ----------
        alloy_list : list of str
            List of alloy compositions.
        api_key : str
            API key for MPRester.

        Returns
        -------
        list of PDEntry
            List of phase diagram entries with intermetallic data.

        """
		pd_entries_list = []
		with MPRester(api_key, mute_progress_bars=True) as mpr:
			gga = mpr.materials.thermo.search(chemsys=alloy_list,
											  fields=['composition', 'formation_energy_per_atom'],
											  thermo_types=[ThermoType.GGA_GGA_U])
			for i in gga:
				energy = i.formation_energy_per_atom
				name = Composition(i.composition)
				pd_entries_list.append(PDEntryLocal(composition=name,
											   energy=energy * name.num_atoms,
											   name=f'{name.alphabetical_formula}_MP'))
		return pd_entries_list
