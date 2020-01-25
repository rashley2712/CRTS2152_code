#!/usr/bin/env python3
import sys, os, argparse, copy
import configlib, spectrumlib
import matplotlib.pyplot, numpy, scipy
import generallib, datetimelib


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads a RV .dat file and plots it.')
	parser.add_argument('inputFile', type=str, help='RV data.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	
	arg = parser.parse_args()

  	# Set up the matplotlib environment
	generallib.setMatplotlibDefaults()
	params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
			}
	matplotlib.rcParams.update(params)    
	plotWidth = 8
	plotHeight = plotWidth/1.62
    

	ephem = datetimelib.ephemerisObject()
	ephem.loadFromFile('ephem.dat')

	rvData = []
	inputFile = open(arg.inputFile, 'rt')
	for line in inputFile:
		line = line.strip()
		if line[0]=='#': continue
		fields = line.split(',')
		entry = {}
		entry['HJD'] = float(fields[0])
		entry['rv'] = float(fields[3])
		entry['rv_err'] = float(fields[4])
		good = int(fields[6])
		if good:
			rvData.append(entry)
		
	inputFile.close()

	phases = []
	rvs = []
	rv_errs = []
	for r in rvData:
		print(r)
		phases.append(ephem.getPhase(r['HJD']))
		rvs.append(r['rv'])
		rv_errs.append(r['rv_err'])
	phasePlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	phasePlot.canvas.set_window_title('RV - phase')
	
	matplotlib.pyplot.errorbar(phases, rvs, color='k', yerr=rv_errs, fmt = '.', ecolor='gray', capsize=0)
	matplotlib.pyplot.errorbar([1+p for p in phases], rvs, color='k', yerr=rv_errs, fmt = '.', ecolor='gray', capsize=0)
	
	def sine(x, K2, gamma):
		return gamma + K2 * numpy.sin(2*numpy.pi*x)

	guess = [400, 0]
	result, covariance = scipy.optimize.curve_fit(sine, phases, rvs, guess, rv_errs)
	errors = numpy.sqrt(numpy.diag(covariance))
	print(result)
	K2 = result[0]
	K2_err = errors[0]
	gamma = result[1]
	gamma_err = errors[1]
	xFit = numpy.arange(0, 2, 0.01)
	yFit = sine(xFit, K2, gamma)
	print("K2 %f (%f)   \t gamma %f (%f)"%(K2, K2_err, gamma, gamma_err))
	matplotlib.pyplot.plot(xFit, yFit, color='g', ls=':')
	matplotlib.pyplot.plot([0,2], [gamma, gamma], color='gray', ls='--')
	matplotlib.pyplot.gca().set_xlim(left = 0, right=2.0)
	matplotlib.pyplot.xlabel('Orbital phase')
	matplotlib.pyplot.ylabel('Radial velocity (km\,s$^{-1}$)')
	matplotlib.pyplot.draw()
	if arg.save is not None:
		matplotlib.pyplot.savefig(arg.save)
	matplotlib.pyplot.show()

	
	sys.exit()

	filenames = []
	if arg.list:
		# Load the list of files.
		if len(arg.inputFiles)>1:
			print("You can only give me one list of filenames.")
			sys.exit()
		filename = arg.inputFiles[0]
		fileList = open(filename, 'r')
		for line in fileList:
			filenames.append(str(line.strip()))
	else:
		filenames = arg.inputFiles

	spectra = []
	for fileIndex, f in enumerate(filenames):
		spectrum = spectrumlib.spectrumObject()
		spectrum.loadFromJSON(f)
		print("%d: \t%s, contains object: %s."%(fileIndex+1, f, spectrum.objectName))
		spectra.append(spectrum)

	numSpectra = len(spectra)

	# Trim out the important region of the spectrum
	print("Discarding all info outside of the range %f to %f Angstroms."%(config.lowerWavelength, config.upperWavelength))
	for s in spectra:
		print(s.HJD)
		s.trimWavelengthRange(config.lowerWavelength, config.upperWavelength)

	rvInfo = rvdata()
	plotWidth = 8
	plotHeight = plotWidth/1.62
	spectrumOverview = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	spectrumOverview.canvas.set_window_title('Overview')
	plotWidth = 6
	plotHeight = plotWidth/1.62
	fitZoom = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	fitZoom.canvas.set_window_title('Fit')
	velocity = 0
	velocityError = 0
	wavelength = 0
	wavelengthError = 0
	reducedChiSq = 0
	for s in spectra:
		wavelengthGuess = config.blueRestWavelength
		repeatFit = True
		while repeatFit:
			wavelengths = s.getWavelengths()
			fluxes = s.getFlux()
			fluxErrors = s.getFluxErrors()
			matplotlib.pyplot.figure(spectrumOverview.number)
			matplotlib.pyplot.title(s.HJD)
			matplotlib.pyplot.step(wavelengths, fluxes, where='mid', color='blue')
			matplotlib.pyplot.step(wavelengths, fluxErrors, color='red', where='mid')
			axes = matplotlib.pyplot.gca()
			(yLower, yUpper) = axes.get_ylim()
			matplotlib.pyplot.plot([config.fitLower, config.fitLower], [yLower, yUpper], linestyle=':', color='grey')
			matplotlib.pyplot.plot([config.fitUpper, config.fitUpper], [yLower, yUpper], linestyle=':', color='grey')
			matplotlib.pyplot.draw()

			# Get the mean value of the continuum outside of the fitting region
			continuumSpectrum = copy.deepcopy(s)
			continuumSpectrum.snipWavelengthRange(config.fitLower, config.fitUpper)
			continuumMean = numpy.mean(continuumSpectrum.flux)
			matplotlib.pyplot.plot([numpy.min(wavelengths), numpy.max(wavelengths)], [continuumMean, continuumMean], linestyle=':', color='green')

			# Extract just the small region of the spectrum required for the fit
			fitSpectrum = copy.deepcopy(s)
			fitSpectrum.trimWavelengthRange(config.fitLower, config.fitUpper)
			wavelengths = fitSpectrum.getWavelengths()
			fluxes = fitSpectrum.getFlux()
			fluxErrors = fitSpectrum.getFluxErrors()
			depthGuess = -1.0 * ((continuumMean - numpy.min(fluxes)) * 0.8)

			# Plot the spectrum
			matplotlib.pyplot.figure(fitZoom.number)
			matplotlib.pyplot.step(wavelengths, fluxes, where='mid', color='blue')

			# Fit the doubleGaussian
			a0 = continuumMean        			# Constant
			a1 = 0								# Slope
			a2 = depthGuess						# Depth
			a3 = wavelengthGuess				# Position of blue line
			bounds = ( [0, -0.2, -a0, config.fitLower], [2*a0, 0.2, 0, config.fitUpper] )
			# print("Bounds:", bounds)
			separation = config.separation
			width = config.width
			def doubleGaussian(x, a0, a1, a2, a3):
				global width, separation
				w = width
				s = separation
				y = a0 + a1 * (x-a3) + a2 * numpy.exp(-.5 * ((x-a3)/w)**2) + a2 * numpy.exp(-.5 * (((x-(a3+s))/w)**2) )
				return y

			goodValues = numpy.full(len(wavelengths), True, dtype=bool)
			rejectedPoints = 0
			newRejectedPoints = True
			iteration = 0
			while newRejectedPoints:
				# Plot the starting value
				xValues = numpy.arange(numpy.min(wavelengths), numpy.max(wavelengths), 0.1)
				yValues = doubleGaussian(xValues, a0, a1, a2, a3)
				matplotlib.pyplot.plot(xValues, yValues, color='green', linestyle=':')

				# Plot the fit and the spectrum
				matplotlib.pyplot.figure(fitZoom.number)
				matplotlib.pyplot.step(wavelengths, fluxes, where='mid', color='blue')

				(yLower, yUpper) = matplotlib.pyplot.gca().get_ylim()
				matplotlib.pyplot.step(wavelengths, [fe+yLower for fe in fluxErrors], color='red', where='mid')
				matplotlib.pyplot.show(block = False)

				guess = [a0, a1, a2, a3]
				x = numpy.array(wavelengths)
				y = numpy.array(fluxes)
				ye = numpy.array(fluxErrors)
				try:
					results, covariance = scipy.optimize.curve_fit(doubleGaussian, x[goodValues], y[goodValues], guess, ye[goodValues], absolute_sigma = True, bounds=bounds)
				except ValueError:
					print("Fit failed. Try to tweak")
					newRejectedPoints = False
					continue
				errors = numpy.sqrt(numpy.diag(covariance))
				(a0, a1, a2, a3) = results
				wavelength = a3
				wavelengthError = errors[3]
				print("Centroid blueward wavelength %f [%f] A"%(wavelength, wavelengthError))
				velocity = (wavelength - config.blueRestWavelength)/config.blueRestWavelength * 3E5
				velocityError = 3E5 / config.blueRestWavelength * wavelengthError
				print("Velocity %f [%f] km/s"%(velocity, velocityError))

				# Plot the fitted curve
				xValues = numpy.arange(numpy.min(wavelengths), numpy.max(wavelengths), 0.1)
				yValues = doubleGaussian(xValues, a0, a1, a2, a3)
				matplotlib.pyplot.plot(xValues, yValues, color='green', linestyle='-')

				# Calculate the chiSquared
				chiSq = 0
				sigma = 0
				for w, flux, fluxError in zip(x[goodValues], y[goodValues], ye[goodValues]):
					fittedFlux = doubleGaussian(w, a0, a1, a2, a3)
					chiSq+= ((flux - fittedFlux)/fluxError)**2
				print("Chi squared", chiSq)
				reducedChiSq = chiSq / (numpy.sum(goodValues) - 4)
				print("Reduced Chi squared", reducedChiSq)

				# Calculate the residuals

				residuals = numpy.array([abs(f - doubleGaussian(w, a0, a1, a2, a3)) for w, f, fe in zip(wavelengths, fluxes, fluxErrors)])
				scaledResiduals = [r/(fe*numpy.sqrt(reducedChiSq)) for r, fe in zip(residuals, fluxErrors)]
				for w, r, rs, fe in zip(wavelengths, residuals, scaledResiduals, fluxErrors):
					print(w, r, rs, fe)
				worstResidual = numpy.max(numpy.array(scaledResiduals)[goodValues])
				index = numpy.argmax(numpy.array(scaledResiduals)[goodValues])
				print("Worst residual: ", wavelengths[index], worstResidual, scaledResiduals[index])
				if worstResidual < config.sigma:
					newRejectedPoints = False
				else:
					goodValues[index] = False

				# goodValues = scaledResiduals < numpy.full(len(wavelengths), config.sigma)
				goodResiduals = residuals[goodValues]
				matplotlib.pyplot.scatter(x[goodValues], [r + yLower for r in goodResiduals], marker='+', color='green')
				badResiduals = residuals[numpy.logical_not(goodValues)]
				matplotlib.pyplot.scatter(x[numpy.logical_not(goodValues)], [r + yLower for r in badResiduals], marker='x', color='red')

				iteration+= 1
				print("Iteration: %d    %d rejected points lying outside %2.1f sigma."%(iteration, numpy.sum(numpy.logical_not(goodValues)), config.sigma))
				if iteration > 10:
					print("The fit is not converging, rejecting fit.")
					break
				# generalUtils.query_yes_no("Continue?")
			matplotlib.pyplot.draw()

			phase = ephem.getPhase(s.HJD)
			#matplotlib.pyplot.figure(phasePlot.number)
			#matplotlib.pyplot.gca().set_xlim(left=0, right=1)
			#matplotlib.pyplot.scatter(phase, velocity)
		
			
			# Pause for input
			print("Phase is %f, Radial velocity is %f (%f)"%(phase, velocity, velocityError))
			print("Happy with the fit? ([Y], [n] or retry with a tweak to the [l]eft or [r]ight. ")
			choice = input().lower()

			if choice=='n':
				goodFit = False
				repeatFit = False
			elif choice=='l':
				wavelengthGuess-= 4
				repeatFit = True
			elif choice=='r':
				wavelengthGuess+= 4
				repeatFit = True
			else:
				goodFit = True
				repeatFit = False


			# Clear the old plots
			spectrumOverview.clf()
			fitZoom.clf()
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show(block=False)
			matplotlib.pyplot.pause(0.01)

		saveRV = { 'HJD': s.HJD,
					'vel': velocity, 'velError': velocityError,
					'wavelength': wavelength, 'wavelengthError': wavelengthError,
					'redChiSq': reducedChiSq,
					'good': goodFit }
		rvInfo.addData(saveRV)
		rvInfo.setFilename(s.objectName)
		rvInfo.sort()
		rvInfo.saveDataTSV()
		rvInfo.saveData()
