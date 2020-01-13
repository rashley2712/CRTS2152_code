#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads ASAS-SN data (csv) and plots it.')
	parser.add_argument('inputFile', type=str, help='Name of the input file.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	parser.add_argument('-j', '--tojson', type=str, help="Dump object to a JSON file.")
	arg = parser.parse_args()

	blocking = False
	if arg.interactive: blocking = True


	photometryLoader = photometrylib.loadPhotometry(debug=False)

	object = photometryLoader.loadFromASASSNJSON(arg.inputFile)
	object.id = "CRTS2152"
	object.duplicateColumn('hjd', 'HJD')
	object.dateColumn = 'HJD'
	if not object.loadEphemeris(filename='ephem.dat'):
		print("Need ephemeris to plot phase")
		noEphem = True
	else: 
		noEphem = False
	
	# Check for 'zero' errors
	dates = object.getColumn(object.dateColumn)
	flux = object.getColumn(object.fluxColumn)
	flux_err = object.getColumn(object.fluxErrorColumn)
	for d in object.data: 
		if d['mag_err']==0: 
			mag_err = 2.5 * (d['flux_err']/numpy.log(10)/d['flux'])
			print(d," flux error!")
			print("new err", mag_err)
			d['mag_err'] = mag_err
		
		
	
	
	dates = object.getColumn(object.dateColumn)
	mags = object.getColumn(object.fluxColumn)
	mag_errors = object.getColumn(object.fluxErrorColumn)
	filters = object.getColumn('filter')
	
	if arg.tojson is not None:
		object.writeToJSON(arg.tojson)

	# Set up the matplotlib environment
	generallib.setMatplotlibDefaults()
	params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
			}
	matplotlib.rcParams.update(params)
	plotWidth = 8
	plotHeight = plotWidth/1.62
	photometryPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	
	# Plot the light curve	
	# Split the data from the different filters
	filterNames = ['V', 'g']
	colours = ['g', 'b']
	markers = ['^', '.']
	for filter , colour, marker in zip(filterNames, colours, markers):
		xValues = []
		yValues = []
		yErrors = []
		for (x, y, e, f) in zip(dates, mags, mag_errors, filters):
			if f == filter:
				xValues.append(x)
				yValues.append(y)
				yErrors.append(e)
		print(filter, colour, len(xValues))
		# Filter upper and value measurements
		xValues, yValues, yErrors = numpy.array(xValues), numpy.array(yValues), numpy.array(yErrors)
		matplotlib.pyplot.errorbar(xValues, yValues, color=colour, yerr=yErrors, fmt = '.', ecolor='gray', capsize=0, marker = marker)
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.xlabel("HJD")
	matplotlib.pyplot.ylabel("ASAS-SN magnitude")
	matplotlib.pyplot.draw()
	if arg.save is not None:
		savename = "long-" + arg.save
		print("Writing to file: %s"%savename)
		matplotlib.pyplot.savefig(savename)

	object.calcPhase()
	
	# Check for 'zero' errors
	dates = object.getColumn(object.dateColumn)
	phases = object.getColumn('phase')
	flux = object.getColumn(object.fluxColumn)
	flux_err = object.getColumn(object.fluxErrorColumn)
	for d in object.data: 
		if d['mag_err']==0: 
			mag_err = 2.5 * (d['flux_err']/numpy.log(10)/d['flux'])
			print(d," flux error!")
			print("new err", mag_err)
			d['mag_err'] = mag_err
		
		
	

	if not noEphem:
		object.debug = False
		#object.phaseBin(numBins = 500)

		binnedPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		
		dates = object.getColumn(object.dateColumn)
		phases = object.getColumn('phase')
		flux = object.getColumn(object.fluxColumn)
		flux_err = object.getColumn(object.fluxErrorColumn)

		matplotlib.pyplot.errorbar(phases, flux, color='k', yerr= flux_err, fmt = '.', ecolor='gray', capsize=0)
		matplotlib.pyplot.errorbar([p+1 for p in phases], flux, color='k', yerr= flux_err, fmt = '.', ecolor='gray', capsize=0)
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.draw()
		if arg.save is not None:
			savename = "phase-" + arg.save
			print("Writing to file: %s"%savename)
			matplotlib.pyplot.savefig(savename)
	matplotlib.pyplot.xlabel("Orbital phase")
	matplotlib.pyplot.ylabel("ASAS-SN magnitude")
		
	matplotlib.pyplot.show()

	# Find potential outbursts


	for HJD, mag, mag_error in zip(dates, mags, mag_errors):
		if(mag<14.5):
			print(HJD, mag, mag_error)

	print("earliest date: ", min(dates))
	print("latest date: ", max(dates))
	totalDuration = max(dates) - min(dates)
	print("Total duration %f days or %f years"%(totalDuration, totalDuration/365))
	
	

	
			
