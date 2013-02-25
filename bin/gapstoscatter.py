#!/usr/bin/env python

import matplotlib.pyplot as plt
from argparse import ArgumentParser

class CSVGapFile(object):
    def __init__( self, filename ):
        self.filename = filename
        # Will be a list of start, stop next to each other
        self.xaxis = {'LowCoverage': [], 'Gap': []}
        # Will just be a list of numbers which coorespond to the lines in the gapfile
        self.yaxis = {'LowCoverage': [], 'Gap': []}
        # Line names. Stores the line number with the first column value in the gapfile
        self.names = []
        self.numsamples = 0

        self.colors = { 'LowCoverage': '#000055', 'Gap': '#BA0000' }

    def parse( self ):
        # What about incorrect lines
        with open( self.filename ) as fh:
            for line in fh:
                try:
                    line = line.strip().split( ',' )
                    self.names.append( line[0] )
                    lineno = len( self.names )
                    # Loop through each range section(start, stop, type)
                    for i in range( 1, len( line ), 3 ):
                        start, end, rtype = line[i:i+3]
                        # Effectively put a point on every spot between start and end(inclusive)
                        for x in range( int( start ), int( end ) + 1 ):
                            self.xaxis[rtype].append( x )
                            self.yaxis[rtype].append( lineno )
                    # Keep track of how many samples are processed
                    self.numsamples += 1
                except ValueError:
                    continue

    def makeScatter( self, outputfile ):
        # Init some useful variables
        xmin = 0
        xmax = 2300
        xstep = 100
        ymin = 0
        ymax = self.numsamples
        ystep = 10

        fig, ax = plt.subplots()

        # Plot out the points
        for rtype, values in self.xaxis.iteritems():
            ax.scatter( values, self.yaxis[rtype], s = 5, c = self.colors[rtype], marker = 'o', edgecolors='none', label=rtype )

        ax.set_title( 'Gap/Low Coverage for %s' % self.filename )
        ax.set_xlabel( 'Nucleotide Position' )
        ax.set_ylabel( 'Sample Index Number' )

        # Grab the figure object from the plot
        #fig = plt.gcf()
        fig.set_size_inches( 24, 12 )
        plt.tight_layout( )
        plt.minorticks_on()

        # Create the x and y tick locations
        xticks = [i for i in range( xmin - xstep, xmax + xstep, xstep )]
        yticks = [i for i in range( ymin - ystep, ymax + ystep, ystep )]
        plt.xticks( xticks, xticks )
        plt.yticks( yticks, yticks )

        # Get the 
        ca = fig.gca()
        ca.set_xlim( xmin - xstep, xmax + xstep )
        ca.set_ylim( ymin - ystep, ymax + ystep )

        # Place legend in lower left(loc = 3)
        plt.legend( loc=3 )
        plt.savefig( outputfile, dpi=600 )

def ops( ):
    parser = ArgumentParser()

    parser.add_argument( '--csv', dest='csvfile', required=True, help='CSV Gaps file to parse' )
    parser.add_argument( '-o', dest='outputfile', required=False, default='gaps.png', help='Filepath to put output image[Default: ./gaps.png]' )
    
    return parser.parse_args()

def main( ops ):
    g = CSVGapFile( ops.csvfile )
    g.parse()
    g.makeScatter( ops.outputfile )

if __name__ == '__main__':
    ops = ops()
    main( ops )
