from typing import Union
import numpy as np
import pandas as pd
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Element

from calcEnthalpy_package.math_operations.thermo_calculations import ThermoMaths
from convex_hull import ConvexHull


class StabilityDecomposition:
	
	def __init__(self,
				 data: dict,
				 end_member: dict,
				 tm: ThermoMaths,
				 flags: dict,
				 grid_size: int,
				 api_key: str):
		
		self.im_flag = flags['im_flag']
		self.equi_flag = flags['equi_flag']
		self.correction = flags['correction']
		self.tm = tm
		self.data = data
		self.end_member = end_member
		
	
	
