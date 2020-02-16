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
	plotWidth = 7
	plotHeight = plotWidth

	dopplerPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
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
	
	if arg.save:
		filename = generallib.addSuffixToFilename(arg.save, "map")
		print("Writing file: %s"%filename)
		matplotlib.pyplot.savefig(filename)			
	
	matplotlib.pyplot.show(block=arg.p)


