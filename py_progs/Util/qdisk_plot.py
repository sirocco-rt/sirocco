#!/usr/bin/env python 
'''
    University of Southampton -- JM -- July 2016

                qdisk_plot.py

Synopsis:
    make a plot of the quantities in the disk.diag file

Usage:
    Either import as a module in a python session e.g.
    import qdisk_plot 

    or run from the command line e.g. 

    qdisk_plot.py root 

 
Arguments:
    root 
        root filename to analyse
'''

import pylab as p
import os, sys 
import numpy as np
from astropy.io import ascii

def qdisk_plot(root):

	# some labels
	ylabels = ["Heating", r"$N_{\mathrm{hit}}$", r"$N_{\mathrm{hit}}/N_{\mathrm{tot}}$", 
	r"$T_{\mathrm{heat}}$", r"$T_{\mathrm{irrad}}$", r"$W_{\mathrm{irrad}}$"]

	log_lin = [1,0,0,1,1,1]


	p.figure(figsize=(9,10))

	disk_diag = "diag_%s/%s.disk.diag" % (root, root)

	# read the disk_diag file
	a = ascii.read(disk_diag)

	# cyce through the physical quantities and plot for each annulus
	for j, name in enumerate(a.colnames[3:]):
		
		p.subplot(3,2,j+1)
		p.plot(a[name], ls="steps", c="k", linewidth=2)
		p.ylabel(ylabels[j])
		p.xlabel("Annulus")

		if log_lin[j]:
			p.semilogy()

	p.savefig("qdisk_%s.png" % root, dpi=300)

if __name__ == "__main__":

	root = sys.argv[1]
	qdisk_plot(root)