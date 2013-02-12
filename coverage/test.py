import glob
import coverage
from wrairlib.fff.alignmentinfo import AlignmentInfo, CoverageRegion
import time
import sys
import os
import os.path

# I'm so terrible at testing
basepath = os.path.dirname( os.path.abspath( __file__ ) )

# Mapping Run Tests
start = time.time()

pdir = glob.glob( '%s/Examples/05*' % basepath )[0]
ac = coverage.AlignmentCoverage( pdir )
ref_regions = ac.find_low_coverage()
known_refregions = {
    'GQ377049_PB1_California04': [CoverageRegion(595, 632, -1), CoverageRegion(1540, 1540, -1), CoverageRegion(1541, 2274, 0), CoverageRegion(1, 1575, 0), CoverageRegion(1576, 1714, -1)],
    'GQ117044_HA_California04': [CoverageRegion( 718, 731, -1 )],
    'FJ969515_PA_California04': [],
    'FJ969512_NP_California04': [],
    'FJ969516_PB2_California04': [CoverageRegion( 1, 162, 0 ), CoverageRegion( 1493, 1606, -1 ), CoverageRegion( 1607, 2280, 0 ), CoverageRegion( 1, 1684, 0 )],
    'FJ969513_MP_California04': [],
    'FJ969514_NS_California04': [],
    'FJ969517_NA_California04': [],
}

for ref in ref_regions:
    for region in zip( ref_regions[ref], known_refregions[ref] ):
        assert region[0] == region[1]

end = time.time()
print "Tested %s in %s seconds" % (pdir, end-start)


# Denovo Run Tests
start = time.time()

pdir = glob.glob( '%s/Examples/08_31*' % basepath )[0]
ac = coverage.AlignmentCoverage( pdir )
ref_regions = ac.find_low_coverage()
known_refregions = {
    'contig00001': [CoverageRegion(1,110,-1), CoverageRegion(1545,1571,-1)],
    'contig00002': [CoverageRegion(1,15,-1), CoverageRegion(1390,1536,-1)],
    'contig00003': [CoverageRegion(1,183,-1), CoverageRegion(1149,1150,-1)]
}

for ref in known_refregions:
    for region in zip( ref_regions[ref], known_refregions[ref] ):
        assert region[0] == region[1]

end = time.time()
print "Tested %s in %s seconds" % (pdir, end-start)
