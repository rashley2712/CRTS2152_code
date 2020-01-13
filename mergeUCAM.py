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
		ccd.debug = False
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
		matplotlib.pyplot.title(c.filterName + "'")
		matplotlib.pyplot.xlabel("HJD")
		matplotlib.pyplot.ylabel("counts")
		matplotlib.pyplot.gca().set_ylim(bottom = 0)
		matplotlib.pyplot.draw()

	# Split the data into 2 nights
	nights = []
	for c in ccds:
		night = c.night
		found = False
		for n in nights:
			if n==night: found = True
		if not found: nights.append(night)
	print("All nights:", nights)

	phasePlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight*1.4))
	colours = ['g', 'r', 'purple', 'brown']
	filterColours = { 	'u' : 'blue',
						'g' : 'green',
						'r' : 'red', 
						'i' : 'purple', 
						'z' : 'brown'}
	offset = 0.1
	midDates = {}
	for index, c in enumerate(ccds):
		print(c.id)
		dates = c.getColumn(c.dateColumn)
		midDate = (min(dates) + max(dates)) / 2
		for i, n in enumerate(nights):
			if n == c.night: 
				offsetIndex = offset*i
				midDates[n] = "%.1f"%midDate
		c.calcPhase()
		c.phaseBin(numBins=500)
		phases = c.phases
		flux = c.bins
		flux_err = c.errors
		colour = filterColours[c.filterName]
		if c.filterName  == 'u': 
			flux = [f/4 for f in flux]
			flux_err = [fe/4 for fe in flux_err]
		# Plot the light curve	
		matplotlib.pyplot.errorbar(phases, [f + offsetIndex for f in flux], color=colour, yerr=flux_err, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3, alpha=0.4)
		matplotlib.pyplot.errorbar([p+1 for p in phases], [f + offsetIndex for f in flux], color=colour, yerr=flux_err, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3, alpha=0.4)
	for i, n in enumerate(nights):
		matplotlib.pyplot.text(0.25, offset*i+0.005, "HJD: " + midDates[n], horizontalalignment='center', fontsize='large')
		
	matplotlib.pyplot.xlabel("Orbital phase")
	matplotlib.pyplot.ylabel("Relative counts")
	matplotlib.pyplot.gca().set_ylim(bottom = 0)
	matplotlib.pyplot.gca().set_xlim(0, 2)
	matplotlib.pyplot.plot([1, 1], matplotlib.pyplot.gca().get_ylim(), color='gray', ls='--')
	for i, n in enumerate(nights):
		matplotlib.pyplot.plot([0, 2], [offset*i, offset*i], ls=':', color='k')
	filterList = ['u', 'g', 'r', 'i', 'z']
	for f in filterList:
		matplotlib.pyplot.errorbar([0], [0], yerr=0, label=f, color = filterColours[f], fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=4, alpha=1.0)
	
	axes = matplotlib.pyplot.gca()
	matplotlib.rcParams['legend.loc'] = 'upper right'
	axes.legend()
	
	matplotlib.pyplot.draw()
	

	if arg.save is not None:
		savename = generallib.addSuffixToFilename(arg.save, "all")
		print("Writing to file: %s"%savename)
		matplotlib.pyplot.savefig(savename)


	zoomPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight*1.4))
	colours = ['g', 'r', 'purple', 'brown']
	filterColours = { 	'u' : 'blue',
						'g' : 'green',
						'r' : 'red', 
						'i' : 'purple', 
						'z' : 'brown'}
	offset = 0.1
	midDates = {}
	startPhase = 0.9
	endPhase = 1.1
	for index, c in enumerate(ccds):
		print(c.id)
		dates = c.getColumn(c.dateColumn)
		midDate = (min(dates) + max(dates)) / 2
		for i, n in enumerate(nights):
			if n == c.night: 
				offsetIndex = offset*i
				midDates[n] = "%.1f"%midDate
		c.calcPhase()
		#c.phaseBin(numBins=500)
		#phases = c.phases
		phase = c.getColumn('phase')
		flux = c.getColumn(c.fluxColumn)
		flux_err = c.getColumn(c.fluxErrorColumn)
		newPhases, newFlux, newFluxErr = [],  [],  []
		# Filter out phases
		for p, f, fe in zip(phase, flux, flux_err):
			ps = p
			if p<0.5: ps = ps+1
			if ps > startPhase and ps < endPhase:
				newPhases.append(ps)
				newFlux.append(f)
				newFluxErr.append(fe)
		phases = newPhases
		flux = newFlux
		flux_err = newFluxErr
		colour = filterColours[c.filterName]
		if c.filterName  == 'u': 
			flux = [f/4 for f in flux]
			flux_err = [fe/4 for fe in flux_err]
		# Plot the light curve	
		matplotlib.pyplot.errorbar(phases, [f + offsetIndex for f in flux], color=colour, yerr=flux_err, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3, alpha=0.4)
		# matplotlib.pyplot.errorbar([p+1 for p in phases], [f + offsetIndex for f in flux], color=colour, yerr=flux_err, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3, alpha=0.4)
	for i, n in enumerate(nights):
		matplotlib.pyplot.text(0.925, offset*i+0.005, "HJD: " + midDates[n], horizontalalignment='center', fontsize='large')
		
	matplotlib.pyplot.xlabel("Orbital phase")
	matplotlib.pyplot.ylabel("Relative counts")
	matplotlib.pyplot.gca().set_ylim(bottom = 0)
	matplotlib.pyplot.gca().set_xlim(startPhase, endPhase)
	matplotlib.pyplot.plot([1, 1], matplotlib.pyplot.gca().get_ylim(), color='gray', ls='--')
	for i, n in enumerate(nights):
		matplotlib.pyplot.plot([0, 2], [offset*i, offset*i], ls=':', color='k')
	filterList = ['u', 'g', 'r', 'i', 'z']
	for f in filterList:
		matplotlib.pyplot.errorbar([0], [0], yerr=0, label=f, color = filterColours[f], fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=4, alpha=1.0)
	
	axes = matplotlib.pyplot.gca()
	matplotlib.rcParams['legend.loc'] = 'upper right'
	axes.legend()
	
	matplotlib.pyplot.draw()
	

	if arg.save is not None:
		savename = generallib.addSuffixToFilename(arg.save, "eclipse")
		print("Writing to file: %s"%savename)
		matplotlib.pyplot.savefig(savename)


	#matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.show()

	sys.exit()
	