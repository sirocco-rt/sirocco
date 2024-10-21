#!/usr/bin/env python 
'''
Test the masterfiles specified in the arrays all point to the right files.

needs python to be compiled.

Usage:
	test_masterfiles.py [SIROCCO VERSION]
'''
import py_read_output as rd 
import subprocess, os, sys

# set env variable
# Do not call this when we're on ReadTheDocs
if not os.environ.get('READTHEDOCS'):
	SIROCCO = os.environ["SIROCCO"]

# change these if you want to test different files. These are all in data/ as of October 2019
macro_files = ['h20', 'h10_hetop_lohe1_standard80', 'h10_standard80', 'h10_hetop_standard80', 'h20_hetop_standard80']
std_files = ['standard80_sn_kurucz', 'standard80']

def run_file(pf_template, masterfile, version=""):
	'''
	runs python version using a given pf_template but with masterfile modified.
	Runs the file as _tmp.pf, then deletes it and all files created.
	'''

	# copy pf and modify data string in OrderedDict
	pf = pf_template
	pf["Atomic_data"] = "data/{}".format(masterfile) 

	# write to file 
	rd.write_pf("_tmp.pf", pf)

	# run process
	process = subprocess.run("{}/bin/py{} -i _tmp.pf".format(SIROCCO, version), shell=True, stdout=subprocess.PIPE)

	subprocess.run("/bin/rm -f _tmp.sig _tmp.out.pf logfile _tmp.out.pf.old _tmp.pf", shell=True, stdout=subprocess.PIPE)
	subprocess.run("/bin/rm -rf diag__tmp", shell=True, stdout=subprocess.PIPE)

	return (process)

def check_run (process, m):
	'''
	check a process object returned 0 for masterfile string m.
	prints info to screen.
	'''
	if (process.returncode != 0): 
		print ("Error for masterfile {} with simple atom model".format(m))	
		print ("Last Line of stdout:", process.stdout.decode('UTF-8').splitlines()[-1])
	else:
		print ("All good.")

	print ("--------")



def run_test(VERSION):
	'''
	run the test. 
	'''

	macro_template = rd.read_pf("{}/examples/extended/cv_macro_benchmark.pf".format(SIROCCO))
	std_template = rd.read_pf("{}/examples/basic/cv_standard.pf".format(SIROCCO))

	subprocess.run("Setup_Py_Dir")


	# check macro atom masterfiles with just the macro atom model
	for i, m in enumerate(macro_files):
		process = run_file (macro_template, m, VERSION)
		print (m, process.returncode, end=": ")
		check_run (process, m)

	# for the standard data sets, check both types of models, since python can always read simple data
	for i, m in enumerate(std_files):
		process = run_file (std_template, m, VERSION)
		print (m, process.returncode, end=": ")
		check_run (process, m)

		process = run_file (macro_template, m, VERSION)
		print (m, process.returncode, end=": ")
		check_run (process, m)

	return (0)

if __name__ == "__main__":
	if len (sys.argv) > 1:
		VERSION = sys.argv[1]
	else: 
		VERSION = ""

	print ("Testing atomic masterfiles...")
	run_test(VERSION)
