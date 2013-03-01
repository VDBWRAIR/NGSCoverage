#!/usr/bin/env python

import matplotlib.pyplot as plt
from argparse import ArgumentParser
from wrairlib.primer import Primer
import sys

class CSVGapFile(object):
    def __init__( self, filename, primerfile=None ):
        self.filename = filename
        self.primerfile = primerfile

        # Will be a list of start, stop next to each other
        self.xaxis = {'LowCoverage': [], 'Gap': []}
        # Will just be a list of numbers which coorespond to the lines in the gapfile
        self.yaxis = {'LowCoverage': [], 'Gap': []}
        self.yaxislabels = []
        # Line names. Stores the line number with the first column value in the gapfile
        # These should be the References that were assembled against...probably
        self.names = []

        self.numsamples = 0
        #ymax is same as numsamples so use property instead
        self.xmin = -10
        self.xmax = -1 

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot( 111 )

        self.linestyle = { 'LowCoverage': {'color':'#000055', 'width':5}, 'Gap': {'color':'#BA0000', 'width':4}, 'Primer': {'color':'#000000', 'width':2} }

    @property
    def ymin( self ):
        if self.primerfile is None:
            return -10
        else:
            return -30

    @property
    def ymax( self ):
        return self.numsamples

    @property
    def xstep( self ):
        return max( self.xmax / 100, 1)

    @property
    def ystep( self ):
        return max((self.ymax - self.ymin) / 10, 1)

    def parse_csv( self ):
        # What about incorrect lines
        with open( self.filename ) as fh:
            # For every line in the file 
            # - Store the name line[0] in self.names for later use
            # - Loop through each coverage region in the line and setup x and y axis regions(tuple( start, end ))
            # - Insert yaxis labels for the names with their ytick
            for line in fh:
                try:
                    line = line.strip().split( ',' )
                    self.names.append( line[0] )
                    lineno = len( self.names )
                    # Loop through each range section(start, stop, type)
                    for i in range( 1, len( line ), 3 ):
                        start, end, rtype = line[i:i+3]
                        start = int( start )
                        end = int( end )
                        # Get the maximum value for xaxis
                        self.xmax = max( self.xmax, end )

                        # Sets the start point
                        self.xaxis[rtype].append( (start, end) )
                        # Sets the endpoint
                        self.yaxis[rtype].append( (lineno, lineno) )
                    self.yaxislabels.append( (lineno, line[0]) )
                    # Keep track of how many samples are processed
                    self.numsamples += 1
                except ValueError:
                    continue

    def plotPrimers( self ):
        """ Plot Primer Sets """
        if self.primerfile is None:
            sys.stderr.write( "You cannot plot the primers without specifying the primer file location\n" )
            sys.exit( -1 )

        primer = Primer( self.primerfile )
        primergenes = primer.unique_genes()
        # Filter out the primer genes not in the names
        genelist = set()
        for gene in primergenes:
            for name in self.names:
                # If the gene is in any of the names then
                # add it to the genelist and break the loop
                if gene in name:
                    genelist.add( gene )
                    break

        if len( genelist ) == 0:
            raise ValueError( "None of the genes in %s were found to match any of the identifiers names in %s" % (str( primergenes ), str( self.names )) )

        genecount = 0
        labelf='Forward Primer'
        labelr='Reverse Primer'
        yaxislabels = []
        for gene in genelist:
            genecount -= 1
            primer_regions = primer.get_merged_primer_regions( gene )
            for region in primer_regions:
                if region.direction == 'R':
                    marker = '<'
                    self.ax.plot( (region.start, region.end), (genecount, genecount), color=self.linestyle['Primer']['color'], lw=self.linestyle['Primer']['width'], label=labelr, marker=marker )
                    labelr = ''
                else:
                    marker = '>'
                    self.ax.plot( (region.start, region.end), (genecount, genecount), color=self.linestyle['Primer']['color'], lw=self.linestyle['Primer']['width'], label=labelf, marker=marker )
                    labelf = ''
                self.xmax = max( self.xmax, region.end )
            yaxislabels.insert( 0, (genecount, gene) )
        self.yaxislabels = yaxislabels + [(0, '')] + self.yaxislabels

    def makeScatter( self ):
        """ Draw the low coverage and gap regions on an axes """
        for rtype in self.xaxis:
            label=rtype
            for x,y in zip( self.xaxis[rtype], self.yaxis[rtype]):
                self.ax.plot( x, y, color=self.linestyle[rtype]['color'], label=label, lw=self.linestyle[rtype]['width'] )
                label=''


    def save_png( self, outputfile, title = 'Gap / Low Coverage Plot' , dpi=300 ):
        yticks = [i for i,l in self.yaxislabels]

        width = 25
        height = max( 0.35 * len( yticks ), 4 )
        self.fig.set_size_inches( width, height )
        #plt.minorticks_on()

        # Create the x and y tick locations
        yticks = [i for i,l in self.yaxislabels]
        plt.yticks( yticks, yticks )
        yaxis = self.ax.get_yaxis()
        yaxis.set_ticklabels( [l for i,l in self.yaxislabels] )

        xticks = [i for i in range( 0, self.xmax + self.xstep, self.xstep )]
        plt.xticks( xticks, xticks, rotation='vertical' )

        self.ax.set_title( 'Gap/Low Coverage for %s' % title )
        self.ax.set_xlabel( 'Nucleotide Position' )
        self.ax.set_ylabel( 'Sample Index Number' )

        #self.ax.set_xlim( self.xmin, self.xmax + self.xstep )

        # Place legend in lower left(loc = 3)
        #plt.legend( loc=3 )
        plt.tight_layout( )
        plt.margins( 0.05, 0.5 )
        plt.savefig( outputfile, dpi=dpi )

def ops( ):
    parser = ArgumentParser()

    parser.add_argument( '--csv', dest='csvfile', required=True, help='CSV Gaps file to parse' )
    parser.add_argument( '-o', dest='outputfile', required=False, default='gaps.png', help='Filepath to put output image[Default: ./gaps.png]' )
    parser.add_argument( '-t', '--title', dest='title', required=False, default=None, help='Title for the scatterplot' )
    parser.add_argument( '-p', '--primer', dest='primerfile', required=False, default=None, help='Path to primer fasta file to plot. If not specified it won\'t be plotted' )
    
    return parser.parse_args()

def main( ops ):
    g = CSVGapFile( ops.csvfile, ops.primerfile )
    g.parse_csv()
    g.makeScatter()
    g.plotPrimers()
    g.save_png( ops.outputfile, ops.title )

if __name__ == '__main__':
    ops = ops()
    main( ops )
