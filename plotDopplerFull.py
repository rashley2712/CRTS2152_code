#!/usr/bin/env python3
import numpy, argparse, sys
import matplotlib.pyplot
from matplotlib import gridspec
from astropy.io import fits
import generallib, datetimelib
import trm.roche

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads Doppler maps and comdat and makes a publishable plot.')
	parser.add_argument('map', type=str, help='Name of the map file.')
	parser.add_argument('input', type=str, help='Name of the input file.')
	parser.add_argument('comdat', type=str, help='Name of the comdat file.')
	parser.add_argument('-b', '--bins', type=int, default=50, help='Number of phase bins to use (default 50).')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', action="store_true", help="block the plots")
	arg = parser.parse_args()

	percentileUpper = 98.5
	percentileLower = 2
	colourMap = 'jet'
	Msol = 2E30  # kg
	Rsol = 696340000 # m
	G = 6.67408E-11 # m^3/kg/s

	M1 = 0.6
	M2 = 0.30
	q = M2/M1
	Porb = 0.1629749071 * 86400
	a = numpy.cbrt(G * (M1+M2)*Msol * Porb* Porb / 4 / numpy.pi/numpy.pi)
	omega = 2 * numpy.pi / Porb
	print("a:",a)
	print("a [Rsol]:", a/Rsol)
	print("w:",omega)
	scale = a*omega/1000
	print("scale:", scale);
		
	# Set up the matplotlib environment
	generallib.setMatplotlibDefaults()
	params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
				'legend.loc': 'upper left',
				'legend.fontsize': 'large'
			}
	matplotlib.rcParams.update(params)
	plotWidth = 6
	plotHeight = plotWidth*1.62

	dopplerPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	gs = gridspec.GridSpec(3, 1, height_ratios=[2, 0.65, 0.65])
	f = matplotlib.pyplot.subplot(gs[0])
	hdulist = fits.open(arg.map)
	print(hdulist.info())
	length = len(hdulist)
	numImages = int((length-1)/2)
	print("HDULength is %d thefore we have %d images"%(length, numImages))
	
	mainImage = 1
	data = hdulist[mainImage].data	
	nx = int(hdulist[mainImage].header['NAXIS1'])
	ny = int(hdulist[mainImage].header['NAXIS2'])
	vxy = float(hdulist[mainImage].header['VXY'])
	hdulist.close()
	vmin, vmax = numpy.median(data) - data.std(), numpy.median(data) + 4 * data.std()

	contrastData = generallib.percentiles(data, percentileLower, percentileUpper)
			
	matplotlib.pyplot.imshow(contrastData, origin='lower', cmap=colourMap, aspect='equal', extent=(-nx*vxy,nx*vxy, -ny*vxy,ny*vxy))
	(vx, vy) = trm.roche.vlobe2(q, n=200)
	matplotlib.pyplot.plot(vx*scale, vy*scale, color='g')
	(vx, vy) = trm.roche.vlobe1(q, n=200)
	matplotlib.pyplot.plot(vx*scale, vy*scale, ls="--", color='g')
	(vx, vy) = trm.roche.vstream(q, step=0.01, vtype=1, n=70)
	matplotlib.pyplot.plot(vx*scale, vy*scale, color='g')
	matplotlib.pyplot.xlabel("$\mathrm{V}_\mathrm{x} (\mathrm{km}\,\mathrm{s}^{-1})$")
	matplotlib.pyplot.ylabel("$\mathrm{V}_\mathrm{y} (\mathrm{km}\,\mathrm{s}^{-1})$")
	axis = matplotlib.pyplot.gca()
	print("Scale:", axis.get_xscale())
	axis.set_autoscale_on(False)
	axis.xaxis.set_label_position('top')
	axis.xaxis.tick_top()
	matplotlib.pyplot.draw()
	
	# Now draw the source data
	matplotlib.pyplot.subplot(gs[1])
	hdulist = fits.open(arg.input)
	print("Loading %s"%arg.input)	
	print(hdulist.info())
	timeData = hdulist[4].data
	HJDs = [ t[0] + 2400000.5 for t in timeData]
	flux = hdulist[1].data
	# Calculate the phases
	ephemeris = datetimelib.ephemerisObject()
	ephemeris.loadFromFile("ephem.dat")
	print(ephemeris)
	print(flux)
	fluxRecords = []
	for index, hjd in enumerate(HJDs):
		fluxRecord={}
		fluxRecord['fluxData'] = flux[index]
		phase = ephemeris.getPhase(hjd)
		#print(hjd, phase)
		fluxRecord['hjd'] = hjd
		fluxRecord['phase'] = phase
		fluxRecords.append(fluxRecord)
	
	fluxRecords = sorted(fluxRecords, key=lambda object: object['phase'], reverse = False)
	wLength = numpy.shape(flux)[1]
	numBins = arg.bins
	binWidth = 1/numBins
	binnedFlux = numpy.zeros((numBins, wLength))
	for b in range(numBins):
		pStart = b*binWidth
		pEnd = (b+1)*binWidth
		fluxesToBin = []
		for i in range(len(fluxRecords)):
			if (fluxRecords[i]['phase'] >= pStart) and (fluxRecords[i]['phase'] < pEnd): fluxesToBin.append(fluxRecords[i]['fluxData'])
		#print("Phase bin: %d: %f, %f contains %d elements."%(b, pStart, pEnd, len(fluxesToBin)))
		fluxEntry = numpy.zeros(wLength)
		if len(fluxesToBin)>0:
			for i in range(len(fluxesToBin)):
				fluxEntry+=fluxesToBin[i]
			fluxEntry/= len(fluxesToBin)
		binnedFlux[b] = fluxEntry
	
	wavelengths = hdulist[3].data
	(startWavelength, endWavelength) = (wavelengths[0][0], wavelengths[0][-1])
	#sourcePlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	contrastFlux = generallib.percentiles(binnedFlux, percentileLower, percentileUpper)
	matplotlib.pyplot.imshow(contrastFlux, origin='lower', cmap=colourMap, aspect=20,  extent=[startWavelength, endWavelength, 0, 1])
	#matplotlib.pyplot.axis('auto')
	matplotlib.pyplot.ylabel("Orbital phase")
	matplotlib.pyplot.xticks([], [])
	hdulist.close()

	# Now draw the comdat data
	hdulist = fits.open(arg.comdat)
	matplotlib.pyplot.subplot(gs[2])
	print("Loading %s"%arg.comdat)	
	print(hdulist.info())
	timeData = hdulist[4].data
	HJDs = [ t[0] + 2400000.5 for t in timeData]
	flux = hdulist[1].data
	# Calculate the phases
	ephemeris = datetimelib.ephemerisObject()
	ephemeris.loadFromFile("ephem.dat")
	#print(ephemeris)
	#print(flux)
	fluxRecords = []
	for index, hjd in enumerate(HJDs):
		fluxRecord={}
		fluxRecord['fluxData'] = flux[index]
		phase = ephemeris.getPhase(hjd)
		#print(hjd, phase)
		fluxRecord['hjd'] = hjd
		fluxRecord['phase'] = phase
		fluxRecords.append(fluxRecord)
	
	fluxRecords = sorted(fluxRecords, key=lambda object: object['phase'], reverse = False)
	wLength = numpy.shape(flux)[1]
	numBins = arg.bins
	binWidth = 1/numBins
	binnedFlux = numpy.zeros((numBins, wLength))
	for b in range(numBins):
		pStart = b*binWidth
		pEnd = (b+1)*binWidth
		fluxesToBin = []
		for i in range(len(fluxRecords)):
			if (fluxRecords[i]['phase'] >= pStart) and (fluxRecords[i]['phase'] < pEnd): fluxesToBin.append(fluxRecords[i]['fluxData'])
		#print("Phase bin: %d: %f, %f contains %d elements."%(b, pStart, pEnd, len(fluxesToBin)))
		fluxEntry = numpy.zeros(wLength)
		if len(fluxesToBin)>0:
			for i in range(len(fluxesToBin)):
				fluxEntry+=fluxesToBin[i]
			fluxEntry/= len(fluxesToBin)
		binnedFlux[b] = fluxEntry
	
	wavelengths = hdulist[3].data
	(startWavelength, endWavelength) = (wavelengths[0][0], wavelengths[0][-1])
	contrastFlux = generallib.percentiles(binnedFlux, percentileLower, percentileUpper)
	matplotlib.pyplot.imshow(contrastFlux, origin='lower', cmap=colourMap, aspect=20, extent=[startWavelength, endWavelength, 0, 1])
	#matplotlib.pyplot.axis('auto')
	#dopplerPlot.tight_layout()
	matplotlib.pyplot.subplots_adjust(wspace=None, hspace=0.1)	
	matplotlib.pyplot.draw()	

	matplotlib.pyplot.ylabel("Orbital phase")
	matplotlib.pyplot.xlabel("Wavelength (\AA)")
		
	if arg.save:
		#filename = generallib.addSuffixToFilename(arg.save, "map")
		filename = arg.save		
		print("Writing file: %s"%filename)
		matplotlib.pyplot.savefig(filename)			
	
	matplotlib.pyplot.show(block=arg.p)


