from emmet.core.thermo import ThermoType
from pymatgen.core.composition import Composition
from ..phase_diagram.pdEntry_local import PDEntryLocal
from pymatgen.ext.matproj import MPRester
from typing import List


class IntermetallicExtractions:
    """
    A class to extract intermetallic phase data from the Materials Project (MP) database.

    This class provides methods to retrieve thermodynamic data, such as formation energies,
    for intermetallic phases in a given alloy system using the Materials Project API. The data
    can be used to generate phase diagrams or perform further thermodynamic analysis.

    Methods:
            get_MP_intermetallic(alloy_list, api_key):
                    Retrieves formation energies and compositions for the specified alloy system
                    from the Materials Project API.
    """

    @staticmethod
    def get_MP_intermetallic(alloy_list: List[str], api_key: str) -> List[PDEntryLocal]:
        """
        Retrieve intermetallic phase data from the Materials Project (MP) API for a given alloy system.

        Args:
                alloy_list (List[str]): A list of alloys (strings) that define the chemical system,
                                                                e.g., ['Fe-Ni', 'Cu-Al'].
                api_key (str): Your API key for accessing the Materials Project database.

        Returns:
                List[PDEntryLocal]: A list of phase diagram entries (PDEntryLocal objects)
                                                        containing compositions and formation energies.

        Example::
        
                alloy_list = ['Fe-Ni', 'Cu-Al']
                api_key = 'your_mp_api_key'

                # Retrieves intermetallic data for the Fe-Ni and Cu-Al systems
                entries = IntermetallicExtractions.get_MP_intermetallic(alloy_list, api_key)
                # entries will contain a list of PDEntryLocal objects with composition and formation energy.
        """

        pd_entries_list = []

        # Using MPRester to access the Materials Project data
        with MPRester(api_key, mute_progress_bars=True) as mpr:
            # Retrieve GGA formation energies and compositions for the specified alloy system
            gga = mpr.materials.thermo.search(
                chemsys=alloy_list,
                fields=["composition", "formation_energy_per_atom"],
                thermo_types=[ThermoType.GGA_GGA_U],
            )

            # Process and store the retrieved data in PDEntryLocal format
            for i in gga:
                energy = i.formation_energy_per_atom
                name = Composition(i.composition)
                pd_entries_list.append(
                    PDEntryLocal(
                        composition=name,
                        energy=energy * name.num_atoms,
                        name=f"{name.alphabetical_formula}_MP",
                    )
                )

        return pd_entries_list
