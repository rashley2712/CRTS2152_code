#!/usr/bin/env python3
import argparse, sys, numpy, copy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads UCAM data (JSON) and plots them.')
	parser.add_argument('inputFiles', type=str, nargs='+', help='Name of the input files.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	arg = parser.parse_args()

	blocking = False
	if arg.interactive: blocking = True

	
	# Set up the matplotlib environment
	generallib.setMatplotlibDefaults()
	params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
			}
	matplotlib.rcParams.update(params)
	plotWidth = 8
	plotHeight = plotWidth/1.62
	
	ccds = []
	for index, f in enumerate(arg.inputFiles):
		print("Loading from %s"%f)
		ccd = photometrylib.target(f)
		ccd.loadFromJSON(f)
		ccd.filterName = f.split('-')[-1].split('.')[0]
		ccd.night = f.split('-')[1]
		ccd.sigmaClip(15, 3)
		print(ccd.night)
		ccd.loadEphemeris(filename='ephem.dat')
		ccds.append(ccd)
	
	for c in ccds:
		print(c.id)
		matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		
		dates = c.getColumn(c.dateColumn)
		flux = c.getColumn(c.fluxColumn)
		flux_err = c.getColumn(c.fluxErrorColumn)
	
		# Plot the light curve	
		matplotlib.pyplot.errorbar(dates, flux,color='k', yerr=flux_err, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3)
		matplotlib.pyplot.title(c.filterName)
		matplotlib.pyplot.xlabel("HJD")
		matplotlib.pyplot.ylabel("counts")
		matplotlib.pyplot.gca().set_ylim(bottom = 0)
		matplotlib.pyplot.draw()

	phasePlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	colours = ['r', 'g', 'b']
	for index, c in enumerate(ccds):
		print(c.id)
		
		dates = c.getColumn(c.dateColumn)
		c.calcPhase()
		phases = c.getColumn('phase')
		flux = c.getColumn(c.fluxColumn)
		flux_err = c.getColumn(c.fluxErrorColumn)
		colour = colours[index % len(colours)]
		# Plot the light curve	
		matplotlib.pyplot.errorbar(phases, flux,color=colour, yerr=flux_err, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3)
		matplotlib.pyplot.errorbar([p+1 for p in phases], flux,color=colour, yerr=flux_err, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3)
	matplotlib.pyplot.xlabel("HJD")
	matplotlib.pyplot.ylabel("counts")
	matplotlib.pyplot.gca().set_ylim(bottom = 0)
	matplotlib.pyplot.draw()
	

	#matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.show()

	sys.exit()
	