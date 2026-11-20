#!/usr/bin/env python
'''
	emissivity_plot.py

	generate a plot of the total matom emissivities, with the last level being the k-> emissivity, 
	for a series of runs. 

	Usage: python emissivity_plot.py root1 [root2 root3...]
'''
import os 
import numpy as np 
import sys
import matplotlib.pyplot as plt

def get_summed_emissivities(root):
	'''
	quick hacky way to extract emissivities from a diag file
	uses grep and sed
	'''

	# get emissivity info out of diag file 
	os.system(r"grep _emiss diag_{}/{}_0.diag > total_emiss_{}.dat".format(root, root, root))

	# remove some text and get columns to be the same length
	os.system(r"sed -i -e 's/Macro Atom level emissivities (summed over cells)://g' total_emiss_{}.dat".format(root))
	os.system(r"sed -i -e 's/(summed over cells)://g' total_emiss_{}.dat".format(root))

	# put into arrays
	absorbed, emissivity = np.genfromtxt("total_emiss_{}.dat".format(root), usecols=(3,5), unpack=True)
	return (absorbed, emissivity)

def make_plot(roots, savename="emiss.png"):
	'''
	make an emissivity plot for a list of root python runs.
	save the plot as savename.
	'''
	plt.figure()

	# loop over root filenames, get emissivities, and plot
	for root in roots:
		absorbed, emissivity = get_summed_emissivities(root)
		Nlevels = len(absorbed)
		n = np.arange(1, Nlevels + 1, 1)
		plt.plot(n, emissivity, "-o", label=root)

	plt.xlabel("$n$", fontsize=18)
	plt.ylabel("Level emissivity, total", fontsize=18) # really a luminosity or energy
	plt.semilogy()
	plt.savefig(savename)

if __name__ == "__main__":
	roots = sys.argv[1:]
	make_plot(roots)
