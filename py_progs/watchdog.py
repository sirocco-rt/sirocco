#!/usr/bin/env python 

'''
During-run parser for error logs from parallel simulation runs.
	
This is the during-processing watchdog code to deal with the error logs in 
SIROCCO's parallel mode. Very similar to py_error.py

usage and arguments:

	python watchdog.py [-f folder] root
	root is the parameter files
	-f folder is the folder the diagfiles are in, defaults to diag_root

returns:

	Error summary of ongoing run printed to screen

130722	JM	Coding began -- initial tests conducted successfully

'''

	
## import modules
import sys, subprocess


#help string for user
def help(flag):
	''' help function for user'''
	if flag=='help' or flag =='-h':
		print ('''usage and arguments:

		python py_error.py [-f folder] root
		root is the parameter files
		-f folder is the folder the diagfiles are in, defaults to diag_root''')
		sys.exit(0)
	else: print ('Watchdog.py: for help type use -h flag')

#strip function for comparing error strings (removes actual numbers from error message)
def strip_error(s):
	'''little routine which strips numerical digits from errors'''
	result_nodigits = ''.join(i for i in s if not i.isdigit())
	result =' '.join( result_nodigits.split() )
	return result


if __name__ == "__main__":		# allows one to run from command line without running automatically with write_docs.py

	# read root pf file from command line
	if len(sys.argv) > 1:
		root = sys.argv[-1]
	else:
		print ('Error: no arguments provided. Exiting.')
		exit (0)

	#print help if needed
	help(root)
		


	# set defaults
	nprocessors = 1
	diagfolder_name = 'diag_' + root


	# check for folder args
	for i in range(len(sys.argv)):
		if sys.argv[i]=="-f":
			diagfolder_name = sys.argv[i+1]


	# list files in folder, put in array 'files'
	command = "ls " + diagfolder_name + "/" + root + "_*.diag &"
	output = subprocess.check_output ( command , shell=True )
	files = output.split("\n")
	if files[-1] == "": files = files[0:-1]


	# number of processors and therefore diag files provided as argument
	nfiles = len(files)
	nprocessors = nfiles



	# create some arrays to store the error logs and increment the counters
	# note that we also count the number of threads an error occurs in

	#2d arrays with index i,j, which store the jth error in thread i, plus the counting array
	error_log_thread = []
	error_count_thread  = []

	#errors that have exceeded error counts for each thread are stored here
	error_nolonger=[]
	error_nolonger_thread=[]

	errors=[]


	#everything now initialised, can collate errors

	for i in range(nfiles):

		filepath = files [i]

		#two temporary arrays for storing error logs in this thread
		#these get appended to the arrays above, error_log_thread
		temp_error_log=[]
		temp_error_count=[]
		temp_errors=[]
		

		# grep for error count errors in the diag file
		command = "grep 'Error: error_count:' " + filepath +"&"
		output = subprocess.check_output ( command , shell=True )

		# split the grep output up into lines
		data = output.split("\n")
		if data[-1] == "": data = data[0:-1] #this is needed because the array has an extra element

		# now cycle through the lines in output and append these message to an array
		for j in range(len(data)):
			error_nolonger.append(data[j])
			error_nolonger_thread.append(i)






		# grep for errors themselves in the diag file
		command = "grep 'Error:' " + filepath +"&"
		output = subprocess.check_output ( command , shell=True )
		#print output, filepath

		# split the grep output up into lines
		data = output.split("\n")
		if data[-1] == "": data = data[0:-1]

		# now cycle through the lines in output
		for j in range(len(data)):
			error_string=strip_error(data[j])
			temp_errors.append(error_string)
			No_log=True
			
			for k in range(len(temp_error_log)):
				if error_string==temp_error_log[k]:
					No_log = False	# we've found a match, increment the counts
					temp_error_count[k]+=1
			if No_log:
				temp_error_log.append ( error_string )
				temp_error_count.append (1)	

		error_log_thread.append(temp_error_log)
		error_count_thread.append(temp_error_count)	
		errors.append(temp_errors)
		

	#1d arrays with index i,j, which store the jth error in thread i, plus the counting array
	error_log_total = []
	error_count_total  = []
	thread_count = []



	#now we cycle over each error log for each thread and count the totals
	for i in range(len(error_log_thread)):
		#cycle other error_log array for 
		for j in range(len(error_log_thread[i])):
			error_string=error_log_thread[i][j]
			error_count=error_count_thread[i][j]
			No_log=True
			for k in range(len(error_log_total)):
				if error_string==error_log_total[k]:
					No_log = False	# we've found a match, increment the counts
					error_count_total[k]+=error_count
					thread_count[k] += 1
			if No_log:
				error_log_total.append ( error_string )
				error_count_total.append (error_count)
				thread_count.append(1)







	n_logs=len(error_log_total)
	n_nolonger=len(error_nolonger)

	# print the output to screen

	# first we print the total errors across all scripts
	print ("WATCHDOG Error summary: "+root)
	print ("Collated errors for " + str(nprocessors) + " processors")
	for i in range(n_logs):
		print ("\t%d -- %s -- %d" % (error_count_total[i], error_log_total[i], thread_count[i]))


	# now the errors for each thread 
	print ('\n-------------------------------\n')
	print ('Errors which will no longer be logged')
	print ('Thread -- Error')
	for i in range(n_nolonger):
		print  ('%d -- %s %s' % (error_nolonger_thread[i], 'Error count:', error_nolonger[i][57:]))


	print ('\n-------------------------------\n')
	for i in range(nprocessors):
		print ('\n')
		print ('Thread %d errors' % i)
		print ('Recurrences --  Description')
		for j in range(len(error_log_thread[i])):
			print ("\t%d -- %s" % (error_count_thread[i][j], error_log_thread[i][j]))

	# all done.
	
				


