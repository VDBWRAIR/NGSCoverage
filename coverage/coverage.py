from wrairlib import *
from wrairlib.parser.exceptions import *
from wrairlib.exceptions1 import *
from wrairlib.fff.fffprojectdir import *
from wrairlib.fff.mappingproject import MappingProject
from wrairlib.fff.qcxls import QCXLS
from wrairlib.util import *
from Bio import SeqIO

from wrairlib.fff.refstatus import RefStatus
from wrairlib.fff.alignmentinfo import AlignmentInfo, CoverageRegion

# Default low coverage threshold
MINREADDEPTH = 10

class AlignmentCoverage(object):
    def __init__( self, projDir, readLowCoverage = None ):
        self.projDir = projDir

        if readLowCoverage:
            self.lct = readLowCoverage
        else:
            self.lct = MINREADDEPTH

        self.set_wanted_identifiers( )

        if self.wanted_idents == 'DENOVO':
            aipath = os.path.join( projDir, 'assembly', '454AlignmentInfo.tsv' )
        else:
            aipath = os.path.join( projDir, 'mapping', '454AlignmentInfo.tsv' )
        self.ai = AlignmentInfo( aipath )

    def find_low_coverage( self ):
        if self.wanted_idents == "DENOVO":
            seqgen = self._low_coverage_assembly( )
        else:
            seqgen = self._low_coverage_mapping( )

        reflc = {}
        for seqalign in seqgen:
            if seqalign.name not in reflc:
                reflc[seqalign.name] = []
            [reflc[seqalign.name].append( region ) for region in self._gaplc_regions( seqalign )]

        return reflc

    def _gaplc_regions( self, seqalign ):
        """ Return only the gap and lowcoverage regions from a list of regions """
        for region in seqalign.regions:
            if region.rtype != 1:
                yield region

    def _low_coverage_mapping( self ):
        for seq in self.ai.seqs:
            reflen = self.wanted_idents.get( seq.name, None )
            # If ref is not in wanted_idents then skip it
            if reflen is None:
                continue
            
            # If reference is not the same size as sequence bases
            # Need to make region for the end
            if seq.regions[-1].end != reflen:
                seq.regions.append( CoverageRegion( seq.regions[-1].end + 1, reflen, 0 ) )

            yield seq

    def _low_coverage_assembly( self ):
        for seq in self.ai.seqs:
            yield seq

    def set_wanted_identifiers( self ):
        """
            Only interested in either references listed at the top
            of 454RefStatus.txt(or all contig names if assembly)

        """
        # Use try catch to distinguish mapping vs assembly
        # Also could just be a bad projDir fed in, but hey lets keep things interesting and
        # assume that isn't the case
        self.wanted_idents = {}
        try:
            rs = RefStatus( os.path.join( self.projDir, 'mapping', '454RefStatus.txt' ) )
            ref = rs.get_likely_reference()

            ref_file = reference_file_for_identifier( ref, self.projDir )
            # Store each reference sequence name with its length
            for seq in SeqIO.parse( ref_file, 'fasta' ):
                self.wanted_idents[seq.id] = len( seq.seq )
        except IOError as e:
            self.wanted_idents = 'DENOVO'

class LowCoverage:
    projDir = None
    ref = None
    crTable = None
    lct = None
    def __init__( self, projDir, reference, lct ):
        """
            >>> lc = LowCoverage( 'Examples/05_11_2012_1_TI-MID51_PR_2305_pH1N1', 'California', 100 )
            >>> print len( lc.crTable )
            154
        """
        self.projDir = projDir
        self.ref = reference
        self.lct = lct
        q = QCXLS( os.path.join( self.projDir, 'mapping', '454MappingQC.xls' ) )
        t = q.getCrossReferenceTable()
        self.crTable = self._filterTable( t[1:] )
        self._sortTable( self.crTable )

    def printTSV( self ):
        print "Reference\tRegion Center\tAverage Alignment Depth\tMinimum Alignment Depth\tMaximum Alignment Depth"
        for x in self.crTable:
            print "\t".join( [x[0],x[1],x[3],x[5],x[6]] )

    def _sort( self, x, y ):
        """
            >>> lc = LowCoverage( 'Examples/05_11_2012_1_TI-MID51_PR_2305_pH1N1', 'California', 100 )
            >>> lc._sort( ['a',0], ['b',1] )
            -1
            >>> lc._sort( ['a',1], ['b',0] )
            -1
            >>> lc._sort( ['b',1], ['a',0] )
            1
            >>> lc._sort( ['b',0], ['a',1] )
            1
            >>> lc._sort( ['a',0], ['a',1] )
            -1
            >>> lc._sort( ['a',1], ['a',0] )
            1
            >>> lc._sort( ['a',0], ['a',0] )
            0
        """
        if x[0] < y[0]:
            return -1
        elif x[0] > y[0]:
            return 1
        else:
            if int( x[1] ) < int( y[1] ):
                return -1
            elif int( x[1] ) > int( y[1] ):
                return 1
        return 0

    def _sortTable( self, table ):
        table.sort( cmp = self._sort )

    def _filterTable( self, table ):
        """
            Filter the given table(2d list) down to entries that have lower than the lct and also
            entries that contain text that matches the reference we are looking for

         Tests:
            >>> d = 'Examples/05_11_2012_1_TI-MID51_PR_2305_pH1N1'
            >>> r = 'California'
            >>> lct = 100
            >>> q = QCXLS( os.path.join( d, 'mapping', '454MappingQC.xls' ) )
            >>> lc = LowCoverage( d , r, lct )
            >>> t = lc._filterTable( q.getCrossReferenceTable()[1:] )
            >>> len( filter( lambda x: r not in x[0] and lct <= x[5], t ) )
            0
        """
        f = lambda x: int( x[5] ) < self.lct and self.ref in x[0]
        return filter( f, table )

if __name__ == '__main__':
    import doctest
    doctest.testmod()
