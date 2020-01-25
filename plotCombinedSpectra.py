#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib, configlib
import spectrumlib, generallib
from astropy.io import fits


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads JSON spectrum file and plots it.')
	parser.add_argument('inputFile', type=str, nargs="+", help='Name of the input file.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-c', '--config', type=str, default="plot.cfg", help="Name of the config file. Default is 'plot.cfg'.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	arg = parser.parse_args()

	config = configlib.configClass(debug=False)
	config.load(arg.config)

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
	#spectrumPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	spectrumPlot,(axes1, axes2)= matplotlib.pyplot.subplots(1,2,sharey=False, sharex=False, figsize=(plotWidth, plotHeight))

	
	spectrumPlot.canvas.set_window_title(config.title)
	minW = 1E8
	maxW = 0
	maxH = 0
	minH = 1E8

	spectra = []
	for filename in arg.inputFile:
		spectrum = spectrumlib.spectrumObject()
		spectrum.loadFromJSON(filename)
		spectrum.convertFluxes()
		if min(spectrum.flux) < minH: minH = min(spectrum.flux)
		if max(spectrum.flux) > maxH: maxH = max(spectrum.flux)
		print(spectrum)
		spectra.append(spectrum)
	
	maxH = maxH*1.1

	spectra[0].trimWavelengthRange(4231, 4980)
	spectra[1].trimWavelengthRange(8023,  8767)
	
	length = 0.3
	offset = .1
	spectrum = spectra[0]
	axes1.step(spectrum.wavelengths, spectrum.flux,  color = 'black', lw=0.5)
	labels = [	{'line': '', 'wavelength':  4686, 'height' : 3.5}, 
	         	{'text': '$\mathrm{He}\\textsc{ii}$', 'wavelength':  4686, 'height' : 3.5},
				{'line': '', 'wavelength':  4861, 'height' : 3.5},
			 	{'text': 'H$\\beta$', 'wavelength':  4861, 'height' : 3.5},
				{'line': '', 'wavelength':  4340, 'height' : 3.5},
				{'text': 'H$\\gamma$', 'wavelength':  4340, 'height' : 3.5},
				{'text': '$\mathrm{He}\\textsc{i}$', 'wavelength':  4921, 'height' : 1.5},
				{'line': '', 'wavelength':  4921, 'height' : 1.5},
				{'text': '$\mathrm{He}\\textsc{i}$', 'wavelength':  4714, 'height' : 1.5},
				{'line': '', 'wavelength':  4714, 'height' : 1.5},
				{'text': '$\mathrm{C}\\textsc{ii}$', 'wavelength':  4268, 'height' : 1.5},
				{'line': '', 'wavelength':  4268, 'height' : 1.5},
				{'text': '$\mathrm{He}\\textsc{i}$', 'wavelength':  4388, 'height' : 1.5},
				{'line': '', 'wavelength':  4388, 'height' : 1.5},
				{'text': '$\mathrm{He}\\textsc{i}$', 'wavelength':  4472, 'height' : 2.0},
				{'line': '', 'wavelength':  4472, 'height' : 2.0},
				{'text': '$\mathrm{He}\\textsc{ii}$', 'wavelength':  4542, 'height' : 1.5},
				{'line': '', 'wavelength':  4542, 'height' : 1.5},
				{'text': 'Bowen blend', 'wavelength':  4641, 'height' : 0.05},
				{'line': '', 'wavelength':  4641, 'height' : 0.65},
				{'text': '$\mathrm{Fe}\\textsc{i}$', 'wavelength':  4416, 'height' : 0.05},
				{'line': '', 'wavelength':  4416, 'height' : 0.65}
				]
	
	for l in labels:
		height = l['height']
		print(l)
		if 'line' in l:
			axes1.plot( [l['wavelength'], l['wavelength']], [(height - length) * 1E-15, height * 1E-15], color='gray')
		else:
			axes1.text( l['wavelength'], (height + offset)*1E-15, l['text'], horizontalalignment='center', fontsize=12)
	axes1.set_xlim(left=min(spectrum.wavelengths), right=max(spectrum.wavelengths))


	spectrum = spectra[1]
	doubledFlux = [ f * 1 for f in spectrum.flux]
	axes2.step(spectrum.wavelengths, doubledFlux,  color = 'black', lw=0.5)
	labels = [	{'line': '', 'wavelength':  8498, 'height' : 1.0}, 
				{'line': '', 'wavelength':  8238, 'height' : 1.0},
				{'text': '$\mathrm{Fe}\\textsc{ii}$', 'wavelength':  8238, 'height' : 1.0},
				{'line': '', 'wavelength':  8542, 'height' : 1.0},
				{'text': '$\mathrm{Ca}\\textsc{ii}$', 'wavelength':  8590, 'height' : 1.0},
				{'line': '', 'wavelength':  8662, 'height' : 1.0}
			]			
	
	for l in labels:
		height = l['height']
		if 'line' in l:
			axes2.plot( [l['wavelength'], l['wavelength']], [(height - length) * 1E-15, height * 1E-15], color='gray')
		else:
			axes2.text( l['wavelength'], (height + offset)*1E-15, l['text'], horizontalalignment='center', fontsize=12)
	
	
	axes2.set_xlim(left=min(spectrum.wavelengths), right=max(spectrum.wavelengths))
	
	axes1.spines['right'].set_visible(False)
	axes2.spines['left'].set_visible(False)
	axes1.yaxis.tick_left()
	#axes1.tick_params(labelright='off')
	#axes2.tick_params(labelleft='off')
	axes2.yaxis.tick_right()

	axes1.set_ylim(bottom=0, top=maxH)
	axes2.set_ylim(bottom=0, top=maxH)

		# matplotlib.pyplot.title(config.title)
	#	axes = matplotlib.pyplot.gca()
	#	axes.set_xlim(left=min(spectrum.wavelengths), right=max(spectrum.wavelengths))
	#	axes.set_ylim(bottom=0, top=maxH)
	#	if index==1:
	# axes2.get_yaxis().set_visible(False)
	d = .015 # how big to make the diagonal lines in axes coordinates
	# arguments to pass plot, just so we don't keep repeating them
	kwargs = dict(transform=axes1.transAxes, color='k', clip_on=False, lw=0.5)
	axes1.plot((1-d,1+d), (-d,+d), **kwargs)
	axes1.plot((1-d,1+d),(1-d,1+d), **kwargs)

	kwargs.update(transform=axes2.transAxes)  # switch to the bottom axes
	axes2.plot((-d,+d), (1-d,1+d), **kwargs)
	axes2.plot((-d,+d), (-d,+d), **kwargs)

	
	spectrumPlot.add_subplot(111, frameon=False)
	# hide tick and tick label of the big axis
	matplotlib.pyplot.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
	matplotlib.pyplot.xlabel(config.xlabel)
	matplotlib.pyplot.ylabel(config.ylabel)
	#spectrumPlot.text(0.5, 0.04, 'common X', ha='center')
	
	#matplotlib.pyplot.xlabel(config.xlabel)

	#axes = spectrumPlot.gca()
	#axes.set_xlim(left=minW, right=maxW)
	# axes.set_ylim(bottom=0)

	labels = []
	maxheight = 12.5
	length = -0.3
	offset = .1
	label = {'text': 'CaII', 'wavelength':  8498, 'height' : 2.0}
	labels.append(label)
	
	for l in labels:
		height = l['height']
		#matplotlib.pyplot.plot( [l['wavelength'], l['wavelength']], [height - length, height], color='black')
		#matplotlib.pyplot.text( l['HJD'], height - offset, l['text'], horizontalalignment='center', fontsize=12)
	
	

	matplotlib.pyplot.draw()
	if arg.save is not None:
		print("Writing to file: %s"%arg.save)
		matplotlib.pyplot.savefig(arg.save)
	
	matplotlib.pyplot.show(block=True)
	matplotlib.pyplot.pause(arg.pause)
			