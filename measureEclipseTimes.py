#!/usr/bin/env python3
import argparse, sys, numpy, copy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads UCAM tries to measure the eclipse times.')
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

	lightCurves = []	
	for index, f in enumerate(arg.inputFiles):
		print("Loading from %s"%f)
		ccd = photometrylib.target(f)
		ccd.loadFromJSON(f)
		ccd.filterName = f.split('-')[-1].split('.')[0]
		ccd.night = f.split('-')[1]
		ccd.debug = False
		ccd.sigmaClip(15, 3)
		print("Filter %s, Night %s"%(ccd.filterName, ccd.night)
		lightCurves.append(ccd)

	
	
	for l in lightCurves:
		print(l.id)
		matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		
		dates = c.getColumn(c.dateColumn)
		flux = c.getColumn(c.fluxColumn)
		flux_err = c.getColumn(c.fluxErrorColumn)
	
		# Plot the light curve	
		matplotlib.pyplot.errorbar(dates, flux,color='k', yerr=flux_err, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3)
		matplotlib.pyplot.title(c.filterName + "'")
		matplotlib.pyplot.xlabel("HJD")
		matplotlib.pyplot.ylabel("counts")
		matplotlib.pyplot.gca().set_ylim(bottom = 0)
		matplotlib.pyplot.draw()



	matplotlib.pyplot.show()

	sys.exit()
	