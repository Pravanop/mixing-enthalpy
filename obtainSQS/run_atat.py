"""
Python routine to run certain ATAT functions originally written in C++ and distributed from
Link:
Citation:
Contributer : Pravan Omprakash
"""

import subprocess
import os
import signal
import time
import re

def run_corrdump(
		corrdump_path: str,
		path_dir: str,
		pair_interaction: float = 2.4 ,
		triplet_interaction: float = 3.0
		) :
	"""
	
	:param corrdump_path:
	:param pair_interaction:
	:param triplet_interaction:
	:return:
	"""
	corrdump = corrdump_path
	os.chdir(path_dir)
	command = f"{corrdump} -ro -noe -nop -clus -rnd -2={pair_interaction} -3={triplet_interaction} -l rndstr.in"
	stdout = subprocess.run([command] , shell = True , check = True)
	return stdout

def run_mcsqs(
		mcsqs_path: str,
		path_dir: str, ) -> subprocess.Popen :
	"""
	
	:param mcsqs_path:
	:return:
	"""
	mcsqs = mcsqs_path
	os.chdir(path_dir)
	command = f"{mcsqs} -rc"
	process = subprocess.Popen(command , shell = True , stdout = subprocess.PIPE)
	time.sleep(5)
	check = True
	start = time.time()
	while check :
		stop = time.time()
		if stop - start >= 60 :
			check = False
			if process.poll() is None :
				os.killpg(os.getpgid(process.pid) , signal.SIGTERM)
		if check_mcsqs(path_dir = path_dir) :
			check = False
			if process.poll() is None :
				os.killpg(os.getpgid(process.pid) , signal.SIGTERM)
	
	return process

def check_mcsqs(path_dir: str) :
	"""
	
	:param path_dir:
	:return:
	"""
	os.chdir(path_dir)
	with open(f"{path_dir}/mcsqs.log" , 'r') as f :
		lines = f.read().splitlines()
	
	if "Objective_function= Perfect_match" in lines :
		return True
	
	else :
		best_obj = lines[-1]
		obj = re.search(pattern = '.*(-[0-9]\.[0-9]+)' , string = best_obj).group(1)
		if obj < -1 :
			return True
		else :
			return False

def runsqs(
		pair_interaction: float ,
		triplet_interaction: float ,
		corrdump_path: str ,
		mcsqs_path: str ,
		path_dir: str
		) :
	"""
	
	:param pair_interaction:
	:param triplet_interaction:
	:return:
	"""
	stdout_corrdump = run_corrdump(
		path_dir = path_dir,
		corrdump_path = corrdump_path ,
		pair_interaction = pair_interaction ,
		triplet_interaction = triplet_interaction
		)
	
	assert stdout_corrdump.returncode == 0
	stdout_mcsqs = run_mcsqs(
		mcsqs_path = mcsqs_path ,
		path_dir = path_dir
		)
	assert stdout_mcsqs.returncode == 0
	
	return True
