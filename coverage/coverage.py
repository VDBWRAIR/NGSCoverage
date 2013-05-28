from roche.newbler import ProjectDirectory
from roche.newbler.util import reference_file_for_identifier
from roche.newbler.fileparsers.alignmentinfo import CoverageRegion

from Bio import SeqIO

import glob
import os.path
import sys

# Default low coverage threshold
MINREADDEPTH = 10

class AlignmentCoverage(object):
    """
    >>> basepath = os.path.dirname( os.path.abspath( __file__ ) )
    >>> pdir = glob.glob( '%s/Examples/05*' % basepath )[0]
    >>> 
    >>> ac = AlignmentCoverage( pdir )
    >>> ref_regions = ac.find_low_coverage()
    >>> known_refregions = {
    ... 'GQ377049_PB1_California04': [CoverageRegion(595, 632, 'LowCoverage'), CoverageRegion(1540, 1540, 'LowCoverage'), CoverageRegion(1541, 1575, 'Gap'), CoverageRegion(1576, 1714, 'LowCoverage')],
    ... 'GQ117044_HA_California04': [CoverageRegion(718, 731, 'LowCoverage')],
    ... 'FJ969515_PA_California04': [],
    ... 'FJ969512_NP_California04': [],
    ... 'FJ969516_PB2_California04': [CoverageRegion(1, 162, 'Gap'), CoverageRegion(1493, 1606, 'LowCoverage'), CoverageRegion(1607, 1684, 'Gap')],
    ... 'FJ969513_MP_California04': [],
    ... 'FJ969514_NS_California04': [CoverageRegion(842, 863, 'Gap')],
    ... 'FJ969517_NA_California04': []
    ... }
    >>> for ref, regions in ref_regions:
    ...   try:
    ...     assert len( regions ) == len( known_refregions[ref] )
    ...     for rregion, kregion in zip( regions, known_refregions[ref] ):
    ...       assert rregion == kregion
    ...   except AssertionError:
    ...       print ref
    ...       print regions
    ...       print known_refregions[ref]
    >>> pdir = glob.glob( '%s/Examples/08_31*' % basepath )[0]
    >>> ac = AlignmentCoverage( pdir )
    >>> ref_regions = ac.find_low_coverage()
    >>> known_refregions = {
    ... 'contig00001': [CoverageRegion(1,110,'LowCoverage'), CoverageRegion(1545,1571,'LowCoverage')],
    ... 'contig00002': [CoverageRegion(1,15,'LowCoverage'), CoverageRegion(1390,1536,'LowCoverage')],
    ... 'contig00003': [CoverageRegion(1,183,'LowCoverage'), CoverageRegion(1149,1150,'LowCoverage')]
    ... }
    >>> for ref, regions in ref_regions:
    ...   if ref not in known_refregions:
    ...     continue
    ...   try:
    ...     assert len( regions ) == len( known_refregions[ref] )
    ...     for rregion, kregion in zip( regions, known_refregions[ref] ):
    ...       assert rregion == kregion
    ...   except AssertionError:
    ...       print ref
    ...       print regions
    ...       print known_refregions[ref]
    """
    def __init__( self, projDir, readLowCoverage = None ):
        self.projDir = projDir
        self.project = ProjectDirectory( projDir )
        self.ai = self.project.AlignmentInfo

        if readLowCoverage:
            self.lct = readLowCoverage
        else:
            self.lct = MINREADDEPTH

        self.set_wanted_identifiers( )

    def find_low_coverage( self ):
        """ Return only low coverage or gap regions """
        # Get merged regions
        mregions = self.ai.merge_regions()
        if self.wanted_idents == "DENOVO":
            return self._low_coverage_assembly( mregions )
        else:
            return self._low_coverage_mapping( mregions )

    def _gaplc_regions( self, regions ):
        """ Return only the gap and lowcoverage regions from a list of regions """
        return [region for region in regions if region.rtype != 'Normal']

    def _low_coverage_mapping( self, mregions ):
        for ref, regions in mregions.iteritems():
            reflen = self.wanted_idents.get( ref, None )
            # If ref is not in wanted_idents then skip it
            if reflen is None:
                continue
            
            # If reference is not the same size as sequence bases
            # Need to make region for the end
            if regions[-1].end != reflen:
                reg = CoverageRegion( regions[-1].end + 1, reflen, 'Gap' )
                regions.append( reg )

            yield (ref, reflen, self._gaplc_regions( regions ))

    def _low_coverage_assembly( self, mregions ):
        for ref, regions in mregions.iteritems():
            yield (ref, None, self._gaplc_regions( regions ))

    def set_wanted_identifiers( self ):
        """
            Only interested in either references listed at the top
            of 454RefStatus.txt(or all contig names if assembly)

        """
        # Use try catch to distinguish mapping vs assembly
        # Also could just be a bad projDir fed in, but hey lets keep things interesting and
        # assume that isn't the case
        self.wanted_idents = {}
        if self.project.project_type == 'mapping':
            rs = self.project.RefStatus
            ref = rs.get_likely_reference()
        else:
            self.wanted_idents = 'DENOVO'
            return

        if ref is None:
            return

        try:
            ref_file = reference_file_for_identifier( ref, self.projDir )
        except IOError as e:
            sys.stderr.write( e.strerror + "\n" )
            sys.exit( -1 )

        # Check in case the reference name in 454RefStatus.txt is different than the refrence identifier
        # in the fasta file
        if ref_file is None:
            return

        # Store each reference sequence name with its length
        for seq in SeqIO.parse( ref_file, 'fasta' ):
            self.wanted_idents[seq.id] = len( seq.seq )

if __name__ == '__main__':
    import doctest
    doctest.testmod()
