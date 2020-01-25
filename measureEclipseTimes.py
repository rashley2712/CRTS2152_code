#!/usr/bin/env python3
import argparse, sys, numpy, copy, scipy.optimize
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
		print("Filter %s, Night %s"%(ccd.filterName, ccd.night))
		print("Egress start %f, egress end %f"%(ccd.egressStart, ccd.egressEnd))
		lightCurves.append(ccd)

	
	
	for l in lightCurves:
		print(l.id)
		matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		dates = l.getColumn(l.dateColumn)
		flux = l.getColumn(l.fluxColumn)
		flux_err = l.getColumn(l.fluxErrorColumn)
		startDate = int(min(dates))
		dates = [d-startDate for d in dates]
		# Plot the light curve	
		matplotlib.pyplot.errorbar(dates, flux,color='k', yerr=flux_err, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3)
		matplotlib.pyplot.title(l.filterName + "'")
		matplotlib.pyplot.xlabel("HJD - %d"%startDate)
		matplotlib.pyplot.ylabel("counts")
		matplotlib.pyplot.gca().set_ylim(bottom = 0)
		matplotlib.pyplot.draw()

		ingressStart = l.ingressStart - startDate
		ingressEnd = l.ingressEnd - startDate
		# Filter out just the ingress
		filteredDates, filteredFlux, filteredFluxErr = [], [], []
		for d, f, fe in zip(dates, flux, flux_err):
			if d > ingressStart and d < ingressEnd:
				filteredDates.append(d)
				filteredFlux.append(f)
				filteredFluxErr.append(fe)
		matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		matplotlib.pyplot.errorbar(filteredDates, filteredFlux, color='k', yerr=filteredFluxErr, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3)
		matplotlib.pyplot.title(l.filterName + "' ingress")
		matplotlib.pyplot.xlabel("HJD - %d"%startDate)
		matplotlib.pyplot.ylabel("counts")
		
		def sigmoid(x, a1, a2, a3, a4, a5):
			y = a1 / (1 + numpy.exp(-a2*(x-a3))) + a4 + a5 * (x - a3)
			return y
	
		# Initial parameters
		drop = filteredFlux[-1] - filteredFlux[0]
		a1 = drop
		print("Drop is [a0]:", a1)
		a2 = 30000.
		print("Sharpness [a2] is", a2)
		a3 = numpy.median(filteredDates)
		print("Midpoint of drop is [a3]", a3)
		a4 = filteredFlux[0]
		print("Initial value is [a4]", a4)
		a5 = 0.
		print("Linear slope [a5] is", a5)
		
		xFit = numpy.arange(filteredDates[0], filteredDates[-1], 0.00001)
		yFit = sigmoid(xFit, a1, a2, a3, a4, a5)
		matplotlib.pyplot.plot(xFit, yFit)
		matplotlib.pyplot.show()
		guess = numpy.array([a1, a2, a3, a4, a5])
		result, covariance = scipy.optimize.curve_fit(sigmoid, filteredDates, filteredFlux, guess, filteredFluxErr)
		errors = numpy.sqrt(numpy.diag(covariance))
		(a1, a2, a3, a4, a5) = result
		yFit = sigmoid(xFit, a1, a2, a3, a4, a5)
		matplotlib.pyplot.plot(xFit, yFit, color ='g')
		egressTime = a3
		egressTimeErr = errors[2]
		#egressTimeErr = 0
		print("ingress time: %f (%f)"%(egressTime, egressTimeErr))
		yLims = matplotlib.pyplot.gca().get_ylim()
		matplotlib.pyplot.plot([egressTime, egressTime], yLims, color='gray', ls='--')
		matplotlib.pyplot.plot([egressTime - egressTimeErr, egressTime - egressTimeErr], yLims, color='gray', ls=':')
		matplotlib.pyplot.plot([egressTime + egressTimeErr, egressTime + egressTimeErr], yLims, color='gray', ls=':')
		matplotlib.pyplot.draw()

		# Filter out just the egress
		filteredDates, filteredFlux, filteredFluxErr = [], [], []
		egressStart = l.egressStart - startDate
		egressEnd = l.egressEnd - startDate
		for d, f, fe in zip(dates, flux, flux_err):
			if d > egressStart and d < egressEnd:
				filteredDates.append(d)
				filteredFlux.append(f)
				filteredFluxErr.append(fe)
		
		matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		matplotlib.pyplot.errorbar(filteredDates, filteredFlux, color='k', yerr=filteredFluxErr, fmt = '.', ecolor='gray', capsize=0, marker = '.', ms=3)
		matplotlib.pyplot.title(l.filterName + "' egress")
		matplotlib.pyplot.xlabel("HJD - %d"%startDate)
		matplotlib.pyplot.ylabel("counts")
	
		def sigmoid(x, a1, a2, a3, a4, a5):
			y = a1 / (1 + numpy.exp(-a2*(x-a3))) + a4 + a5 * (x - a3)
			return y
	
		# Initial parameters
		drop = filteredFlux[-1] - filteredFlux[0]
		a1 = drop
		print("Drop is [a0]:", a1)
		a2 = 30000.
		print("Sharpness [a2] is", a2)
		a3 = numpy.median(filteredDates)
		print("Midpoint of drop is [a3]", a3)
		a4 = filteredFlux[0]
		print("Initial value is [a4]", a4)
		a5 = 0.
		print("Linear slope [a5] is", a5)
		
		xFit = numpy.arange(filteredDates[0], filteredDates[-1], 0.00001)
		yFit = sigmoid(xFit, a1, a2, a3, a4, a5)
		matplotlib.pyplot.plot(xFit, yFit)

		guess = numpy.array([a1, a2, a3, a4, a5])
		result, covariance = scipy.optimize.curve_fit(sigmoid, filteredDates, filteredFlux, guess, filteredFluxErr)
		errors = numpy.sqrt(numpy.diag(covariance))
		(a1, a2, a3, a4, a5) = result
		yFit = sigmoid(xFit, a1, a2, a3, a4, a5)
		matplotlib.pyplot.plot(xFit, yFit, color ='g')
		egressTime = a3
		egressTimeErr = errors[2]
		#egressTimeErr = 0
		print("egress time: %f (%f)"%(egressTime, egressTimeErr))
		yLims = matplotlib.pyplot.gca().get_ylim()
		matplotlib.pyplot.plot([egressTime, egressTime], yLims, color='gray', ls='--')
		matplotlib.pyplot.plot([egressTime - egressTimeErr, egressTime - egressTimeErr], yLims, color='gray', ls=':')
		matplotlib.pyplot.plot([egressTime + egressTimeErr, egressTime + egressTimeErr], yLims, color='gray', ls=':')
		matplotlib.pyplot.draw()
	


	matplotlib.pyplot.show()

	sys.exit()
	