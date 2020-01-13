#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib
import astropy


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads JSON photometry and plots and tries to normalise them.')
	parser.add_argument('inputFiles', type=str, nargs="+", help='Names of the input files.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	parser.add_argument('-j', '--tojson', action="store_true", help="Dump object to a JSON file.")
	arg = parser.parse_args()

	blocking = False
	if arg.interactive: blocking = True

		# Set up the matplotlib environment
	generallib.setMatplotlibDefaults()
	params = {	'axes.labelsize': 'x-large',
				'xtick.labelsize': 'x-large',
				'ytick.labelsize': 'x-large',
				'legend.loc': 'upper left',
				'legend.fontsize': 'large'
			}
	matplotlib.rcParams.update(params)
	plotWidth = 8
	plotHeight = plotWidth/1.62
	

	# Load each set of photometry into an array
	eachSource = []
	for f in arg.inputFiles:
		newJSONdata = photometrylib.target('new')
		newJSONdata.loadFromJSON(f)
		print("%s file contains %d points from source: %s"%(f, len(newJSONdata.getColumn(newJSONdata.dateColumn)), newJSONdata.source))
		eachSource.append(newJSONdata)

	# Split the ASAS-SN data into its two separate filters
	newSources = []
	for s in eachSource:
		if s.source == 'ASAS-SN':
			asassndata = s
			filterList = []
			for d in asassndata.data:
				notFound = True
				for f in filterList:
					if f==d['filter']: notFound = False
				if notFound: filterList.append(d['filter'])
			for f in filterList:
				newTarget = photometrylib.target("from filter %s"%f)
				for d in asassndata.data:
					if d['filter'] == f:
						newTarget.data.append(d)
				newTarget.source="ASAS-SN %s"%f
				newTarget.dateColumn = asassndata.dateColumn
				newTarget.fluxColumn = asassndata.fluxColumn
				newTarget.fluxErrorColumn = asassndata.fluxErrorColumn
				newSources.append(newTarget)
		else:
			newSources.append(s)
		

	eachSource = newSources
	
	# convert the magnitudes to flux units
	for s in eachSource:
		for d in s.data:
			d['fluxFromMag'] = 10**((d[s.fluxColumn]) / -2.5)
			d['fluxErrorFromMag'] = d['fluxFromMag'] * numpy.log(10) * d[s.fluxErrorColumn]
			print(d['fluxFromMag'], d['fluxErrorFromMag'])
	#	s.fluxColumn = 'fluxFromMag'
	#	s.fluxErrorColumn = 'fluxErrorFromMag'
	
	# Plot each data set
	for s in eachSource:
		photometryPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		xValues = s.getColumn(s.dateColumn)
		yValues = s.getColumn(s.fluxColumn)
		yErrors = s.getColumn(s.fluxErrorColumn)
		colour = 'k'
		ecolour = 'gray'
		marker = '.'
		matplotlib.pyplot.errorbar(xValues, yValues, color=colour, yerr=yErrors, fmt = '.', ecolor=ecolour, capsize=0, marker = marker)
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.xlabel(s.dateColumn)
		matplotlib.pyplot.ylabel(s.fluxColumn)
		matplotlib.pyplot.title(s.source)

		matplotlib.pyplot.draw()
	

	# Plot the light curves together
	colours = ['r', 'g', 'b', 'k', 'purple']
	ecolours = ['pink', 'lightgreen', 'lightblue', 'gray']
	markers = ['.', 's', 'v', '^', '+']
	timeStart = 1E10
	timeStop = 0
	mergedPlot = matplotlib.pyplot.figure(figsize=(plotWidth*2, plotHeight))
	for index, s in enumerate(eachSource):
		xValues = s.getColumn(s.dateColumn)
		if max(xValues)>timeStop: timeStop = max(xValues)
		if min(xValues)<timeStart: timeStart = min(xValues)
		yValues = s.getColumn(s.fluxColumn)
		yErrors = s.getColumn(s.fluxErrorColumn)
		colour = colours[index % len(colours)]
		ecolour = ecolours[index % len(ecolours)]
		marker = markers[index % len(markers)]
		matplotlib.pyplot.errorbar(xValues, yValues, color=colour, yerr=yErrors, fmt = '.', ecolor=ecolour, capsize=0, marker = marker, alpha=0.5, label=s.source, ms=4)
	axes = matplotlib.pyplot.gca()
	(bottom, top) = axes.get_ylim()

	labels = []
	maxheight = 12.5
	length = -0.3
	offset = .1
	label = {'text': 'ULTRACAM', 'HJD':  2457917.8, 'height' : 13.0}
	labels.append(label)
	label = {'text': 'HiPERCAM', 'HJD':  2458046.4, 'height' : 14.0}
	labels.append(label)
	label = {'text': 'ULTRACAM', 'HJD':  2458407.6, 'height' : 13.0}
	labels.append(label)
	
	for l in labels:
		height = l['height']
		matplotlib.pyplot.plot( [l['HJD'], l['HJD']], [height - length, height], color='black')
		matplotlib.pyplot.text( l['HJD'], height - offset, l['text'], horizontalalignment='center', fontsize=12)
	
	

	#matplotlib.pyplot.plot([w1mstart, w1mstart], [bottom, top], color='gray', linestyle=':')
	#matplotlib.pyplot.plot([w1mstop, w1mstop], [bottom, top], color='gray', linestyle=':')
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.gca().set_xlim(left= timeStart, right= timeStop+160)
	matplotlib.pyplot.gca().set_ylim(top= maxheight)
	matplotlib.pyplot.xlabel("HJD")
	matplotlib.pyplot.ylabel("Survey magnitude")
	matplotlib.pyplot.legend()	
	if arg.save:
		filename = generallib.addSuffixToFilename(arg.save, "all")
		print("Writing file: %s"%filename)
		matplotlib.pyplot.savefig(filename)
	
	# Do a phase plot
	params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
				'legend.loc': 'upper left',
				'legend.fontsize': 'medium'
			}
	matplotlib.rcParams.update(params)
	plotWidth = 8
	plotHeight = plotWidth/1.62
	
	
	phasePlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	for index, s in enumerate(eachSource):
		s.loadEphemeris('ephem.dat')
		s.calcPhase()
		xValues = s.getColumn('phase')
		yValues = s.getColumn(s.fluxColumn)
		yErrors = s.getColumn(s.fluxErrorColumn)
		colour = colours[index % len(colours)]
		ecolour = ecolours[index % len(ecolours)]
		marker = markers[index % len(markers)]
		matplotlib.pyplot.errorbar(xValues, yValues, color=colour, yerr=yErrors, fmt = '.', ecolor=ecolour, capsize=0, marker = marker, alpha=0.5, label=s.source, ms=4)
		matplotlib.pyplot.errorbar([x+1 for x in xValues], yValues, color=colour, yerr=yErrors, fmt = '.', ecolor=ecolour, capsize=0, marker = marker, alpha=0.5, ms=4)
	axes = matplotlib.pyplot.gca()
	(bottom, top) = axes.get_ylim()
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.gca().set_xlim(left= 0, right= 2)
	#matplotlib.pyplot.gca().set_ylim(top= maxheight)
	matplotlib.pyplot.xlabel("Orbital phase")
	matplotlib.pyplot.ylabel("Survey magnitude")
	matplotlib.pyplot.legend()	
	matplotlib.pyplot.draw()
	if arg.save:
		filename = generallib.addSuffixToFilename(arg.save, "phase")
		print("Writing file: %s"%filename)
		matplotlib.pyplot.savefig(filename)
	matplotlib.pyplot.show()
	
	#sys.exit()
	# Remove 'outburst' data by looking for points x mags above/below the mean
	range = 10.75
	for s in eachSource:
		fluxValues = numpy.array(s.getColumn(s.fluxColumn))
		fluxErrors = numpy.array(s.getColumn(s.fluxErrorColumn))
		weights = 1. / (fluxErrors**2)
		mean = numpy.average(fluxValues, weights=weights)
		print("%s has a mean of %f"%(s.source, mean))
		newData = []
		for d in s.data:
			if abs(d[s.fluxColumn] - mean) > range: continue
			newData.append(d)
		print('Removed %d points outside the range of %f'%(len(s.data)-len(newData), range))
		s.data = newData
	

	# Now find the mean of each source and subtract it to 'normalise' in mags
	for s in eachSource:
		fluxValues = s.getColumn(s.fluxColumn)
		fluxErrors = s.getColumn(s.fluxErrorColumn)
		fluxValues = numpy.array(s.getColumn(s.fluxColumn))
		fluxErrors = numpy.array(s.getColumn(s.fluxErrorColumn))
		weights = 1. / (fluxErrors**2)
		mean = numpy.average(fluxValues, weights=weights)
		print(s.source, mean)
		for d in s.data:
			d[s.fluxColumn] = d[s.fluxColumn] - mean
			

	# Plot again
	colours = ['r', 'g', 'b', 'k']
	ecolours = ['pink', 'lightgreen', 'lightblue', 'gray']
	mergedPlot = matplotlib.pyplot.figure(figsize=(plotWidth*2, plotHeight))
	for index, s in enumerate(eachSource):
		xValues = s.getColumn(s.dateColumn)
		yValues = s.getColumn(s.fluxColumn)
		yErrors = s.getColumn(s.fluxErrorColumn)
		colour = colours[index % len(colours)]
		ecolour = ecolours[index % len(colours)]
		marker = '.'
		matplotlib.pyplot.errorbar(xValues, yValues, color=colour, yerr=yErrors, fmt = '.', ecolor=ecolour, capsize=0, marker = marker, alpha=0.5, label=s.source)
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.xlabel("HJD")
	matplotlib.pyplot.ylabel("mag")
	matplotlib.pyplot.legend()	
	
	# Now do the Lomb-Scargle

	allDates = []
	allFlux = []
	allFluxErrors = []
	for s in eachSource:
		xValues = s.getColumn(s.dateColumn)
		allDates = allDates + xValues
		yValues = s.getColumn(s.fluxColumn)
		allFlux = allFlux + yValues
		yErrors = s.getColumn(s.fluxErrorColumn)
		allFluxErrors = allFluxErrors + yErrors
		
	print("In total there are %d datapoints."%len(allDates))

	from astropy.timeseries import LombScargle
	frequency = numpy.linspace(0.1, 6, 1000)
	# power = LombScargle(allDates, allFlux, allFluxErrors).power(frequency)
	#frequency, power = LombScargle(allDates, allFlux, allFluxErrors).autopower()
	frequency, power = LombScargle(allDates, allFlux, allFluxErrors).autopower(minimum_frequency=0.1, maximum_frequency=20)
	
	lsPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	matplotlib.pyplot.plot(frequency, power) 

	best_frequency = frequency[numpy.argmax(power)]
	print("Best frequency is %f d^-1 or a period of %f days"%(best_frequency, 1/best_frequency))
	best_frequency/=2
	print("Half of that is: %f d^-1 or a period of %f days"%(best_frequency, 1/best_frequency))
	
	matplotlib.pyplot.show()

	# do a quick LS for each data set separately
	for s in eachSource:
		print()
		print("Data for: %s"%s.source)
		xValues = s.getColumn(s.dateColumn)
		yValues = s.getColumn(s.fluxColumn)
		yErrors = s.getColumn(s.fluxErrorColumn)
		frequency, power = LombScargle(xValues, yValues, yErrors).autopower()

		best_frequency = frequency[numpy.argmax(power)]
		print("Best frequency is %f d^-1 or a period of %f days"%(best_frequency, 1/best_frequency))
		best_frequency/=2
		print("Half of that is: %f d^-1 or a period of %f days"%(best_frequency, 1/best_frequency))
		print()
		
	

	firstDate = numpy.min(allDates)
	lastDate = numpy.max(allDates)
	print("HJDs run from %f to %f"%(firstDate, lastDate))
	print("Total time span is %d days or %f years."%(int(lastDate-firstDate), (lastDate - firstDate)/365))

	# Write all data to a text file for further analysis
	outputFile = open("merged.dat", "wt")
	outputFile.write("# HJD, mag, mag_err, source\n")
	for s in eachSource:
		for d in s.data:
			outputFile.write("%f, %f, %f, %s\n"%(d['HJD'], d[s.fluxColumn], d[s.fluxErrorColumn], s.source))

	outputFile.close()


	sys.exit()

