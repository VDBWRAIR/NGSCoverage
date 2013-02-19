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

        self.colors = { 'LowCoverage': '#FDF5E6', 'Gap': '#FF0000' }

    def parse( self ):
        with open( self.filename ) as fh:
            for line in fh:
                line = line.strip().split( ',' )
                self.names.append( line[0] )
                lineno = len( self.names )
                for i in range( 1, len( line ), 3 ):
                    start, end, rtype = line[i:i+3]
                    self.xaxis[rtype].append( start )
                    self.xaxis[rtype].append( end )
                    self.yaxis[rtype].append( lineno )
                    self.yaxis[rtype].append( lineno )

    def makeScatter( self ):
        for rtype, values in self.xaxis.iteritems():
            plt.scatter( values, self.yaxis[rtype], s = 5, c = self.colors[rtype], marker = 'o' )
        fig = plt.gcf()
        fig.set_size_inches( 24, 12 )
        plt.savefig( 'gaps.png', dpi=100 )

def ops( ):
    parser = ArgumentParser()

    parser.add_argument( '--csv', dest='csvfile', required=True, help='CSV Gaps file to parse' )
    
    return parser.parse_args()

def main( ops ):
    g = CSVGapFile( ops.csvfile )
    g.parse()
    g.makeScatter()

if __name__ == '__main__':
    ops = ops()
    main( ops )
