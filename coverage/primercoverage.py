from wrairlib.primer import Primer

from Bio import SeqIO
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

import re
from argparse import ArgumentParser, FileType
import logging

# Default pattern to match Reference files
REF_MATCH_PATTERN = '(?P<name>(?P<accession>.*?)_(?P<gene>.*?)_(?P<strain>.*))'

def logging_setup( level ):
    # Setup basic config for logging
    FORMAT = '%(asctime)-15s %(message)s'
    logging.basicConfig( format=FORMAT, level=getattr( logging, level ) )

def refs( reffile, pattern ):
    '''
        >>> from pprint import pprint
        >>> ref_match_pattern = '(?P<name>(?P<accession>.*?)_(?P<gene>.*?)_(?P<strain>.*))'
        >>> reffile = 'Examples/Ref/H1N1_boston.fasta'
        >>> primerfile = 'Examples/Primer/sH1N1.fasta' >>> r = refs( reffile, ref_match_pattern )
        >>> r == {'NS': 877, 'PB1': 2300, 'PB2': 2314, 'NA': 1446, 'PA': 2215, 'MP': 1002, 'NP': 1552, 'HA': 1750}
        True
    '''
    try:
        cp = re.compile( pattern )
    except:
        raise ValueError( "Incorrect pattern given: %s" % pattern )
    refs = {}
    for seq in SeqIO.parse( reffile, 'fasta' ):
        m = cp.match( seq.id )
        if m:
            refparts = m.groupdict()
            # Index the references by the gene name
            refs[refparts['gene']] = len( str( seq.seq ) )
        else:
            raise ValueError( "%s was in an unkown format" % reffile )

    return refs

def ref_coverage( reffile, primerfile, pattern ):
    '''
        >>> from pprint import pprint
        >>> ref_match_pattern = '(?P<name>(?P<accession>.*?)_(?P<gene>.*?)_(?P<strain>.*))'
        >>> reffile = 'Examples/Ref/H1N1_boston.fasta'
        >>> primerfile = 'Examples/Primer/sH1N1.fasta' 
        >>> rc = ref_coverage( reffile, primerfile, ref_match_pattern )
        >>> isinstance( rc, dict )
        True
    '''
    ref = refs( reffile, pattern )
    primer = Primer( primerfile )

    genes = {}
    count = 0

    # Labels for legend
    labelf = 'Forward Primer'
    labelr = 'Reverse Primer'
    # Ensure iteration in sorted order to ensure consistent iteration everywhere
    for gene in sorted( ref.keys() ):
        logging.debug( "Creating Line2D for %s" % gene )
        regions = primer.get_merged_primer_regions( gene )
        reflen = ref.get( gene, False )
        if not regions:
            raise ValueError( "%s had no regions for gene %s. You may need to refine the match pattern." % (primerfile, gene) )

        # Create a line for the reference segment
        refline = Line2D( [0, reflen], [count, count], linewidth=5, color='green', alpha=0.5 )
        logging.debug( "Reference Line: (0,%s), (%s,%s)" % (reflen, count, count) )
        primerlines = []
        for region in regions:
            primerlines.append( Line2D( [region.start, region.end], [count, count], linewidth=2, color='black' ) )
            if region.direction == 'R':
                primerlines[-1].set_marker( '<' )
                primerlines[-1].set_label( labelr )
                # Set reverse label to None so that only 1 line gets an actual label
                #  otherwise the legend tries to use all of them making it way to big
                labelr = ''
            else:
                primerlines[-1].set_marker( '>' )
                primerlines[-1].set_label( labelf )
                # Same principal as above
                labelf = ''
            logging.debug( "Primer Line: (%s,%s), (%s,%s)" % (region.start, region.end, count, count) )
        genes[gene] = {'Reference': refline, 'Primers': primerlines}
        count += 1

    # Now create just one labeled refrence line for the legend
    first = genes[genes.keys()[0]]
    first['Reference'].set_label( 'Reference' )

    return genes

def max_reflen( lines ):
    '''
        >>> ref_match_pattern = '(?P<name>(?P<accession>.*?)_(?P<gene>.*?)_(?P<strain>.*))'
        >>> reffile = 'Examples/Ref/H1N1_boston.fasta'
        >>> primerfile = 'Examples/Primer/sH1N1.fasta' 
        >>> rc = ref_coverage( reffile, primerfile, ref_match_pattern )
        >>> max_reflen( rc )
        2314
    '''
    return max( [pr['Reference'].get_xdata()[1] for pr in [gene for _, gene in lines.iteritems()]] )

def annotate_reference( line2d, axes ):
    '''
        Annotate the end of reference with its length
    '''
    start, end = line2d.get_xydata()
    logging.debug( "Annotate: Start(%s) End(%s)" % (start, end) )
    arrow = { 'width': 0.5, 'headwidth': 1 }
    axes.annotate( int( end[0] ), end, (end[0]-50, end[1]+0.2), arrowprops=arrow )

def make_image( reffile, primerfile, pattern, title='Primer coverage for Reference', outputname='refcov.png' ):
    '''
        >>> import os
        >>> ref_match_pattern = '(?P<name>(?P<accession>.*?)_(?P<gene>.*?)_(?P<strain>.*))'
        >>> reffile = 'Examples/Ref/H1N1_boston.fasta'
        >>> primerfile = 'Examples/Primer/sH1N1.fasta' 
        >>> make_image( reffile, primerfile, ref_match_pattern )
        >>> os.unlink( 'refcov.png' )
    '''
    # Create figure and add axes
    fig = plt.figure()
    fig.set_size_inches( 20, 8 )
    ax = fig.add_subplot( 111 )

    # Gather all the lines for each gene(Reference and Primers)
    lines = ref_coverage( reffile, primerfile, pattern )

    # Set the xaxis and yaxis limits to fit the values of the lines
    ax.set_ylim( -1, len( lines ) )
    ax.set_xlim( -50, max_reflen( lines ) + 50 )

    # Turn on minor tick marks
    # and then turn them off for y axis
    ax.minorticks_on()
    ax.tick_params( axis='y', which='minor', left='off', right='off' )

    # Label the axis
    ax.xaxis.set_label_text( 'Nucleotide Position' )
    ax.yaxis.set_label_text( 'Reference Segment' )

    # Set the title of the figure
    ax.set_title( title )

    # Turn on grid lines for xaxis
    ax.xaxis.grid( b=True )

    # ylabels will store the labels for each Gene on the yaxis
    # Need to include one below and one above the actual data
    # Has to be sorted to ensure iteration happens in the same order
    genes = sorted( lines.keys() )
    ylabels = [''] + genes + ['']
    ax.yaxis.set_ticks( range( -1, len( lines ) ) )
    ax.yaxis.set_ticklabels( ylabels )

    # Map every line onto the axes
    for gene in genes:
        rp = lines[gene]
        logging.debug( "Drawing %s" % gene )
        map( ax.add_line, rp['Primers'] )
        ax.add_line( rp['Reference'] )
        annotate_reference( rp['Reference'], ax )

    # Add the legend for all items that have labels
    ax.legend( loc='lower center', fancybox=True, shadow=True, ncol=3 )

    # Output the image to the filesystem
    fig.savefig( outputname, dpi=175 )

if __name__ == '__main__':
    import doctest
    doctest.testmod()
