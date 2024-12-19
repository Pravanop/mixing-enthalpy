
import pickle

with open('../database/intermetallic_database_bokas.pickle', 'rb') as handle:
	im_list = pickle.load(handle)
	
im_names = [str(i.name).split('_')[0] for i in im_list]\

#take out numbers using regex
import re

im_names = [i for i in im_names for c in i if not c.isdigit()]

for i in im_names:
	if 'Nb' in i and 'V' in i:
		print(i)
# print(im_names)