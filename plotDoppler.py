#!/usr/bin/env python3
import numpy, argparse, sys
import matplotlib.pyplot
from astropy.io import fits
import generallib
import trm.roche

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads JSON photometry and plots and tries to normalise them.')
	parser.add_argument('inputFile', type=str, help='Names of the input files.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")

	arg = parser.parse_args()

	Msol = 2E30  # kg
	Rsol = 696340000 # m
	G = 6.67408E-11 # m^3/kg/s

	M1 = 0.7
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
	params = {	'axes.labelsize': 'x-large',
				'xtick.labelsize': 'x-large',
				'ytick.labelsize': 'x-large',
				'legend.loc': 'upper left',
				'legend.fontsize': 'large'
			}
	matplotlib.rcParams.update(params)
	plotWidth = 8
	plotHeight = plotWidth/1.62

	hdulist = fits.open(arg.inputFile)
	print(hdulist.info())
	length = len(hdulist)
	numImages = int((length-1)/2)
	print("HDULength is %d thefore we have %d images"%(length, numImages))
	
	for index in range(numImages):
		num = index*2 + 1
		data = hdulist[num].data	
		nx = int(hdulist[num].header['NAXIS1'])
		ny = int(hdulist[num].header['NAXIS2'])
		vxy = float(hdulist[num].header['VXY'])
		
		vmin, vmax = numpy.median(data) - data.std(), numpy.median(data) + 4 * data.std()

		contrastData = generallib.percentiles(data, 5, 98)
		#plt.imshow(data, cmap='viridis', aspect='equal', vmin=vmin, vmax=vmax, extent=(-nx*vxy,nx*vxy, -ny*vxy,ny*vxy))
		dopplerPlot = matplotlib.pyplot.figure(figsize=(plotHeight, plotHeight))
			
		matplotlib.pyplot.imshow(contrastData, origin='lower', cmap='viridis', aspect='equal', extent=(-nx*vxy,nx*vxy, -ny*vxy,ny*vxy))
		#matplotlib.pyplot.imshow(data, origin='lower', cmap='viridis', aspect='equal', vmin=vmin, vmax=vmax, extent=(-nx*vxy,nx*vxy, -ny*vxy,ny*vxy))
		(vx, vy) = trm.roche.vlobe2(q, n=200)
		matplotlib.pyplot.plot(vx*scale, vy*scale, color='g')
		(vx, vy) = trm.roche.vlobe1(q, n=200)
		matplotlib.pyplot.plot(vx*scale, vy*scale, ls="--", color='g')
		(vx, vy) = trm.roche.vstream(q, step=0.01, vtype=1, n=80)
		matplotlib.pyplot.plot(vx*scale, vy*scale, color='g')
		matplotlib.pyplot.draw()
		if arg.save:
			if numImages>1: filename = generallib.addSuffixToFilename(arg.save, str(index))
			else: filename = arg.save
			print("Writing file: %s"%filename)
			matplotlib.pyplot.savefig(filename)			
	matplotlib.pyplot.show(block=False)
