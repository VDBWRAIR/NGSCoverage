# Default low coverage threshold
lct = 100

class LowCoverage:
    projDir = None
    ref = None
    crTable = None
    lct = None
    def __init__( self, projDir, reference, lct ):
        """
            >>> lc = LowCoverage( 'Examples/05_11_2012_1_TI-MID51_PR_2305_pH1N1', 'California', 100 )
            >>> lc.printTSV()
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

