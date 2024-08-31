import numpy as np
import json

class thermo_calc():
	def __init__(self, binary):

		omegas_dict = json.load(open('./bokas_omegas_processed.json'))
		end_members_dict = json.load(open('./bokas_end_members_dict.json'))
