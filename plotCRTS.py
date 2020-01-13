#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib
from trm.aspec import amp_spec 
from astropy.timeseries import LombScargle
from scipy.optimize import leastsq
import scipy.optimize

	


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads CRTS data and plots it.')
	parser.add_argument('inputFile', type=str, help='Name of the input file.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	parser.add_argument('-j', '--json', type=str, help="Save to a JSON file (specify filename).")

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
    
	photometryLoader = photometrylib.loadPhotometry(debug=True)

	targets = photometryLoader.loadFromCRTS(arg.inputFile)
	for t in targets:
		t.id = "CRTS2152"
		if t.loadEphemeris('ephem.dat'): t.hasEphemeris = True
		print(t.ephemeris)
		t.computeHJDs()
		if (t.ephemeris.orbit): 
			t.hasOrbit = True
			t.calcPhase()	
		else: t.hasOrbit = False
		t.dateColumn = 'HJD'
		if arg.json is not None:
			t.writeToJSON(arg.json)
	print(len(targets))

	for object in targets:
		photometryPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		xValues = object.getColumn('HJD')
		yValues = object.getColumn('mag')
		yErrors = object.getColumn('err')
		matplotlib.pyplot.errorbar(xValues, yValues, color='k', yerr=yErrors, fmt = '.', ecolor='gray', capsize=0)
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.xlabel('HJD')
		matplotlib.pyplot.ylabel('CRTS magnitude')
		matplotlib.pyplot.gca().set_xlim(left=min(xValues), right=max(xValues))
		if arg.save is not None:
			filename = generallib.addSuffixToFilename(arg.save, "full")
			print("Writing to file: %s"%filename)
			matplotlib.pyplot.savefig(filename)
		
		matplotlib.pyplot.draw()
		

	t = targets[0]
	# Do a Periodogram
	startFrequency  = 1.1  # Cycles/day
	stopFrequency  = 20   # Cycles/day
	spacing = 0.0001
	numsamples = int((stopFrequency-startFrequency)/spacing)
	
	f0,df,Nf = startFrequency, spacing, numsamples
	startTime = numpy.min(xValues)
	times = numpy.array([ x - startTime for x in xValues ])
	yMean = numpy.mean(yValues)
	c = numpy.array([y - yMean for y in yValues])
	ce = yErrors

	fs, amps = amp_spec(times,c,ce,f0,df,Nf)
	fmax = fs[numpy.argmax(amps)]
	print('Frequency of maximum amplitude =',fmax,'cycles/day  period=',1/fmax)
	print('Period of maximum amplitude =',86400/fmax,'seconds or', 24/fmax, 'hours')
	
	frequency = numpy.linspace(0.1, 6, 1000)
	#frequency, power = LombScargle(allDates, allFlux, allFluxErrors).autopower()
	frequency, power = LombScargle(times, c, ce).autopower(minimum_frequency=0.1, maximum_frequency=20)
	fmax = frequency[numpy.argmax(power)]
	print('Frequency of maximum amplitude =',fmax,'cycles/day  period=',1/fmax)
	print('Period of maximum amplitude =',86400/fmax,'seconds or', 24/fmax, 'hours')
	periodogramPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	matplotlib.pyplot.plot(fs, amps)
	matplotlib.pyplot.plot(frequency, power, alpha=0.4)
	matplotlib.pyplot.draw()

	def sineFit(x, a, b, c, f):
		phase = 2*numpy.pi*f*x
		return c + a*numpy.cos(phase) + b*numpy.sin(phase)

	def model(p,t):
		c,f,a,b = p
		phase = 2*numpy.pi*f*t
		return c + a*numpy.cos(phase) + b*numpy.sin(phase)

	def res(p,t,c):
		return c-model(p,t)

	def cosHarmonics(x, c, a1, freq, t0):
		return c + a1 * numpy.cos( 2*numpy.pi*freq*(x-t0) )

	p0 = (1,fmax,0,0)
	ret,cov = leastsq(res,p0,(times,c))
	p0 = ret
	print("Covariance", cov)
	fmax = p0[1]
	Pbest = 1/fmax
	print('Fitted period =',86400*Pbest,'seconds')
	print('Fitted amplitude =',numpy.sqrt(p0[2]**2+p0[3]**2))

	guess = [0, 0, 1, fmax]
	result, covariance = scipy.optimize.curve_fit(sineFit, times, c, guess, ce)
	a, b, c, f = result
	print("Result", result)
	errors = numpy.sqrt(numpy.diag(covariance))
	frequency = result[3]
	frequencyError = errors[3]
	phase = numpy.arctan(-a / b) + numpy.pi/2.
	print("Phase (radians)", phase)
	print("Fit frequency: %f (%f)"%(frequency, frequencyError))
	periodDays = 1/frequency
	periodError = frequencyError/(frequency*frequency)
	print("Fit period: %.10f (%.10f)"%(periodDays, periodError))
	T0 = startTime + phase/2/numpy.pi * periodDays
	T0_error = phase/2/numpy.pi * periodError
	print("T0: %.10f (%.10f)"%(T0, T0_error))

	# Do the fit all over again

	x = t.getColumn('HJD')
	y = t.getColumn('mag')
	ye = t.getColumn('err')
	guess = (numpy.mean(y), -1, fmax, startTime)
	result, covariance = scipy.optimize.curve_fit(cosHarmonics, x, y, guess, ye)
	errors = numpy.sqrt(numpy.diag(covariance))
	print("New result:", result)
	constant, a1, frequency, T0 = result
	con_err, a1_err, frequencyError, T0_error = errors
	periodDays = 1/frequency
	periodError = frequencyError/(frequency*frequency)
	T0-= periodDays/2. 
	print("Fit period: %.10f (%.10f)"%(periodDays, periodError))
	print("T0: %.10f (%.10f)"%(T0, T0_error))
	

	if not t.hasOrbit: 
		sys.exit()
	
	# Update the ephemeris
	for t in targets:
		t.ephemeris.Period = periodDays
		t.ephemeris.Period_error = periodError
		t.ephemeris.T0 = T0
		t.ephemeris.T0_error = T0_error
		t.calcPhase()	
	
	
	for object in targets:
		dates = object.getColumn(object.dateColumn)
		totalDuration = max(dates) - min(dates)
		startDate = dates[0]
		phasePlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		xValues = object.getColumn('phase')
		yValues = object.getColumn('mag')
		yErrors = object.getColumn('err')
		
		#matplotlib.pyplot.scatter(xValues, yValues, c=colourLookup, cmap='gist_rainbow')
		#matplotlib.pyplot.scatter([x+1 for x in xValues], yValues, c=colourLookup, cmap='gist_rainbow')
		matplotlib.pyplot.errorbar(xValues, yValues, color='k', yerr=yErrors, fmt = '.', ecolor='gray', capsize=0)
		matplotlib.pyplot.errorbar([x+1 for x in xValues], yValues, color='k', yerr=yErrors, fmt = '.', ecolor='gray', capsize=0)
		matplotlib.pyplot.gca().invert_yaxis()
	
		def cosPhase(x, a1, c):
			return c + a1*numpy.cos(2*numpy.pi*(x+0.5)) 
		# Plot the fit
		xFit = numpy.arange(0, 2, 0.03)
		yFit = cosPhase(xFit, a1, constant)
		matplotlib.pyplot.plot(xFit, yFit)
		matplotlib.pyplot.plot([1, 1], matplotlib.pyplot.gca().get_ylim(),  ls=':')
	
		matplotlib.pyplot.xlabel('Orbital phase')
		matplotlib.pyplot.ylabel('CRTS magnitude')
		matplotlib.pyplot.gca().set_xlim(left=0, right=2)
		matplotlib.pyplot.draw()
		if arg.save is not None:
			filename = generallib.addSuffixToFilename(arg.save, "phase")
			print("Writing to file: %s"%filename)
			matplotlib.pyplot.savefig(filename)
		matplotlib.pyplot.show()

		dates = object.getColumn(object.dateColumn)
		totalDuration = max(dates) - min(dates)
		print("total length %f days or %f years."%(totalDuration, totalDuration/365))

			
