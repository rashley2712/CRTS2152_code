#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Merges DR2 and DR3 data.')
	parser.add_argument('inputFile', nargs="+", type=str, help='Name of the input files.')
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

	dataSets = []
	for filename in arg.inputFile:
		d = photometrylib.target(filename)
		d.loadFromJSON(filename)
		print(len(d.data))
		dataSets.append(d)



	photometryPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	
	colours = ['r', 'g', 'b']
	offset = 2
	for index, object in enumerate(dataSets):
		xValues = object.getColumn('HJD')
		yValues = object.getColumn('mag')
		yErrors = object.getColumn('err')
		colour = colours[ index % len(colours)]
		matplotlib.pyplot.errorbar(xValues, [y+index*offset for y in yValues] , color=colour, yerr=yErrors, fmt = '.', ecolor=colour, capsize=0)
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.xlabel('HJD')
		matplotlib.pyplot.ylabel('CRTS magnitude')
		# matplotlib.pyplot.gca().set_xlim(left=min(xValues), right=max(xValues))
		if arg.save is not None:
			filename = generallib.addSuffixToFilename(arg.save, "full")
			print("Writing to file: %s"%filename)
			matplotlib.pyplot.savefig(filename)
		
		matplotlib.pyplot.draw()
		

		dates = object.getColumn(object.dateColumn)
		totalDuration = max(dates) - min(dates)
		print("total length %f days or %f years."%(totalDuration, totalDuration/365))

			
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show()