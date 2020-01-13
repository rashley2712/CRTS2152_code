#!/usr/bin/env python3
import argparse, sys, numpy, copy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads UCAM data (logfile) and plots it.')
	parser.add_argument('inputFile', type=str, help='Name of the input file.')
	parser.add_argument('telescope', type=str, help='Name of the telescope.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	parser.add_argument('-j', '--tojson', type=str, help="Dump object to a JSON file.")
	arg = parser.parse_args()

	blocking = False
	if arg.interactive: blocking = True

	photometryLoader = photometrylib.loadPhotometry(debug=True)

	objects = photometryLoader.loadFromHCAM(arg.inputFile)
	
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
	
	for o in objects:
		print(o.id)
		o.loadEphemeris(filename='ephem.dat')
		o.telescope = arg.telescope
		errorCount = 0
		errorCodes = o.getColumn('errorflag')
		savedErrors = []
		for e in errorCodes:
			if e != 0: 
				errorCount+=1
				savedErrors.append(e)
		if errorCount>0:
			print("Errors:", errorCount)
			print(savedErrors)
			#sys.exit()
		o.computeJDfromMJD()
		o.computeHJDfromJD()
		o.dateColumn="HJD"
	
		dates = o.getColumn(o.dateColumn)
		flux = o.getColumn(o.fluxColumn)
		flux_err = o.getColumn(o.fluxErrorColumn)
	
		# Plot the light curve	
		matplotlib.pyplot.errorbar(dates, flux, color='g', yerr=flux_err, fmt = '.', ecolor='lightgreen', capsize=0, marker = '.')


	#matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.xlabel("HJD")
	matplotlib.pyplot.ylabel("counts")
	matplotlib.pyplot.draw()
	
	# Now merge the object and the comparison by dividing the flux
	newObjects = []
	for index in range(0, len(objects), 2):
		newTarget = copy.deepcopy(objects[index])
		comparison = objects[index+1]
		targetFlux = newTarget.getColumn(newTarget.fluxColumn)
		comparisonFlux = comparison.getColumn(comparison.fluxColumn)
		targetFluxError = newTarget.getColumn(newTarget.fluxErrorColumn)
		comparisonFluxError = comparison.getColumn(comparison.fluxErrorColumn)
		newFlux = [tf/cf for tf, cf in zip(targetFlux, comparisonFlux)]
		newFlux_err = [ nf * numpy.sqrt((tfe/tf)**2 + (cfe/cf)**2) for nf, tf, tfe, cf, cfe in zip(newFlux, targetFlux, targetFluxError, comparisonFlux, comparisonFluxError) ]
		
		newTarget.setColumn('comparisonCounts', comparisonFlux)
		newTarget.setColumn('comparisonCountsError', comparisonFluxError)
		newTarget.setColumn('relativeCounts', newFlux)
		newTarget.fluxColumn = 'relativeCounts'
		newTarget.setColumn('relativeCountsError', newFlux_err)
		newTarget.fluxErrorColumn = 'relativeCountsError'
		newTarget.id = "CRTS2152-%d"%(int(index/2))

		newObjects.append(newTarget)
		
		relativeFluxPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		dates = newTarget.getColumn(newTarget.dateColumn)
		matplotlib.pyplot.errorbar(dates, newFlux, color='k', yerr=newFlux_err, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3)
		matplotlib.pyplot.title(newTarget.id)
		matplotlib.pyplot.draw()

	if arg.save is not None:
		savename = arg.save
		print("Writing to file: %s"%savename)
		matplotlib.pyplot.savefig(savename)

	if arg.tojson is not None:
		for o in newObjects:
			filename = o.id + ".json" 
			o.writeToJSON(filename)

	# Some extra data
	for o in newObjects:
		print(o.id)
		dates = o.getColumn(o.dateColumn)
		start = min(dates)
		end = max(dates)
		print("Total length: %f days %f minutes."%((end-start), (end-start)*24*60))
		print("Mid run %f"%((start+end)/2))
		print("Exposures: %d  %fs per exposure."%(len(dates), (end-start)*86400/len(dates)))

	matplotlib.pyplot.show()
	