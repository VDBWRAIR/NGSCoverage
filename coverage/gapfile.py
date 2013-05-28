import matplotlib.pyplot as plt
from wrairlib.primer import Primer

class CSVGapFile(object):
    def __init__( self, filename, primerfile=None ):
        self.filename = filename
        self.primerfile = primerfile

        # Will be a list of start, stop next to each other
        self.xaxis = {'LowCoverage': [], 'Gap': [], 'Reference': []}
        # Will just be a list of numbers which coorespond to the lines in the gapfile
        self.yaxis = {'LowCoverage': [], 'Gap': [], 'Reference': []}
        self.yaxislabels = []
        self._xaxislabels = { }
        self.xaxislabels = 0
        # Line names. Stores the line number with the first column value in the gapfile
        # These should be the References that were assembled against...probably
        self.names = []

        self.numsamples = 0
        #ymax is same as numsamples so use property instead
        self.xmin = -10
        #xmax is property

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot( 111 )

        self.linestyle = { 'LowCoverage': {'color':'#000055', 'width':5}, 'Gap': {'color':'#BA0000', 'width':4}, 'Primer': {'color':'#000000', 'width':2}, 'Reference': {'color': '#5F8700', 'width': 1} }

    @property
    def xaxislabels( self ):
        return self._xaxislabels.keys()

    @xaxislabels.setter
    def xaxislabels( self, value ):
        self._xaxislabels[value] = str(value)

    @property
    def xmax( self ):
        return self.__dict__.get( 'xmax', -1 )
    @xmax.setter
    def xmax( self, value ):
        value = int( value )
        if value > self.__dict__.get('xmax',-1):
            self.__dict__['xmax'] = value

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

    def split_namelen( self, in_namelen ):
        namelen = in_namelen.split('|')
        if len( namelen ) != 2 or namelen[1] == '':
            raise ValueError( "{} is incorrectly formatted. Does not contain |reflen".format(in_namelen) )
        return namelen

    def parse_csv( self ):
        # What about incorrect lines
        with open( self.filename ) as fh:
            # For every line in the file 
            # - Store the name line[0] in self.names for later use
            # - Loop through each coverage region in the line and setup x and y axis regions(tuple( start, end ))
            # - Insert yaxis labels for the names with their ytick
            for line in fh:
                try:
                    lineno = len( self.names ) + 1
                    line = line.strip().split( ',' )
                    refname, reflen = self.split_namelen( line[0] )
                    self.xmax = reflen
                    self.names.append( refname )
                    self.xaxis['Reference'].append( (0, int(reflen)) )
                    self.yaxis['Reference'].append( (lineno, lineno ) )
                    # Loop through each range section(start, stop, type)
                    for i in range( 1, len( line ), 3 ):
                        start, end, rtype = line[i:i+3]
                        start = int( start )
                        end = int( end )
                        # Sets the line so it starts and ends on the same line
                        self.xaxis[rtype].append( (start, end) )
                        self.yaxis[rtype].append( (lineno, lineno) )
                        # Uses descriptor above to add to dictionary
                        self.xaxislabels = start
                        self.xaxislabels = end
                    self.yaxislabels.append( (lineno, line[0]) )
                    # Keep track of how many samples are processed
                    self.numsamples += 1
                except ValueError as e:
                    print e
                    continue
        if len( self.xaxis['Gap'] ) == 0 and len( self.xaxis['LowCoverage'] ) == 0:
            raise EOFError( "Empty csv file" )

        # Now ensure that (0,0) is plotted on xaxis
        self.xaxis['Gap'].insert( 0, (0,0) )
        self.yaxis['Gap'].insert( 0, (0,0) )

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
                self.xmax = region.end
            yaxislabels.insert( 0, (genecount, gene) )
        self.yaxislabels = yaxislabels + [(0, '')] + self.yaxislabels

    def makeScatter( self ):
        """ Draw the low coverage and gap regions on an axes """
        for rtype in self.xaxis:
            label=rtype
            for x,y in zip( self.xaxis[rtype], self.yaxis[rtype]):
                xp = list( x )
                # If the gap is a single base, need to extend
                # the size by 1 so it actually shows up
                if x[1] - x[0] == 0 and x[0] != 0:
                    xp[1] = x[1] + 1
                self.ax.plot( xp, y, color=self.linestyle[rtype]['color'], label=label, lw=self.linestyle[rtype]['width'] )
                label=''

    def get_axis_margin( self, axisinstance, inchesdesired ):
        '''
            Given an axis intstance and inches desired change the margin for that axis
            so it matches the inches desired
        '''
        # Get the interval to work with
        s, e = axisinstance.get_data_interval()
        datainterval = float( e - s )
        if datainterval == 0:
            datainterval = 1
        # Return the ratio to be used with axes.set_?margin()
        try:
            return inchesdesired / datainterval
        except ZeroDivisionError as e:
            sys.stderr.write( "Axis data interval was: %s\n" % axisinstance.get_data_interval() )
            raise e

    def save_png( self, outputfile, title = 'Gap / Low Coverage Plot' , dpi=100 ):
        yaxis = self.ax.get_yaxis()
        xaxis = self.ax.get_xaxis()

        yticks = [i for i,l in self.yaxislabels]

        width = 25
        # Lower limit for height is 4 inches, upper limit is 14 inches
        height = min( max( 0.35 * len( yticks ), 4 ), 20 )
        self.fig.set_size_inches( width, height )

        # Create Y axis labels.
        # Crazy fun python here
        plt.yticks( *zip( *self.yaxislabels ) )
        yaxis.grid( color='grey' )

        xticks = [i for i in range( 0, self.xmax + self.xstep, self.xstep )]
        plt.xticks( xticks, xticks, rotation='vertical' )

        usetitle = 'Gap/Low Coverage'
        if title:
            usetitle += ' -- %s' % title
        self.ax.set_title( usetitle )
        self.ax.set_xlabel( 'Nucleotide Position' )
        #self.ax.set_ylabel( 'Sample Index Number' )
        plt.tight_layout( pad=3 )

        # Shink current axis's height by 10% on the bottom
        box = self.ax.get_position()
        self.ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])
        # Put a legend below current axis
        self.ax.legend(loc='upper center', bbox_to_anchor=(0.2, -0.125), fancybox=True, shadow=True, ncol=4)

        ymargin = self.get_axis_margin( yaxis, 1 ) 
        plt.margins( 0, ymargin )
        plt.savefig( outputfile, dpi=dpi )
