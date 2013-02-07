from wrairlib.fff.refstatus import RefStatus
from wrairlib import *
from wrairlib.parser.exceptions import *
from wrairlib.exceptions1 import *
from wrairlib.fff.fffprojectdir import *
from wrairlib.fff.mappingproject import MappingProject
from wrairlib.util import *
import os
import os.path
import sys
from Bio import SeqIO
import re

class FindGaps:
    dir = None # Holds the main directory
    sample_dirs = [] # Holds all of the midkey directories
    reference_path = None # Holds the reference path
    reference_info = None # Holds the SeqIO.index for the reference
    parsed_reads = {} # Holds the parsed reads
    use_only_ref = "" # Holds the reference name
    single = False # Is this a single dir run

    def __init__( self, directory, use_only_ref ):
        if directory[-1] == '/':
            self.dir = directory[:-1]
        else:
            self.dir = directory
        self.use_only_ref = use_only_ref
        self._check_single()
        self._set_sample_dirs()
        self._set_reference( use_only_ref )
        self._set_reference_info( )

    def _check_single( self ):
        contents = list_dir( self.dir )
        for d in contents:
            if os.path.basename( d ) == '454Project.xml':
                self.single = True
                break

    def _set_sample_dirs( self ):
        """
            Sets the sample directories from the given main directory


            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/', 'california' )
            >>> len( fg.sample_dirs )
            48

            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2010_01_15/CoInfection_reAnalysis/', 'california' )
            >>> len( fg.sample_dirs )
            14
        """
        contents = list_dir( self.dir )
        if self.single:
            self.sample_dirs = [ os.path.basename( self.dir ) ]
        else:
            dirs = []
            for c in contents:
                if os.path.isdir( c ) and ('MID' in c or 'RL' in c):
                    dirs.append( c )
            self.sample_dirs = dirs

    def _set_reference( self, reference = None ):
        """
            Set the reference using the given reference or infer it from 454RefStatus.txt
        
            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/', 'manag' )
            >>> fg.reference_path
            '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/Ref/H3N2_Managua.fasta'
            
            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/', None )
            >>> fg.reference_path
            '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/Ref/H3N2_Managua.fasta'
        """
        # If the reference is blank then it indicates
        # that it should be automatically found using 454RefStatus file
        if reference is None:
            if self.single:
                refstatus = RefStatus( os.path.join( self.dir, 'mapping', '454RefStatus.txt' ) )
            else:
                refstatus = RefStatus( os.path.join( self.dir, self.sample_dirs[0], 'mapping', '454RefStatus.txt' ) )
                
            reference = refstatus.get_likely_reference()
            self._find_reference_by_sequencename( reference )
        else:
            # Return the searched out reference
            self._find_reference_by_filename( reference )

    def _find_reference_by_sequencename( self, seqname ):
        """ Given a sequence identifier name return the reference file that contains it in a 454project """
        refs = self._refs_for_proj()
        for ref in refs:
            for seq in SeqIO.parse( open( ref ), 'fasta' ):
                if seqname.lower() in seq.id.lower():
                    self.reference_path = ref

    def _refs_for_proj( self ):
        # Setup the mapping project parser
        if self.single:
            first_proj = os.path.join( self.dir, 'mapping', '454MappingProject.xml' )
        else:
            first_proj = os.path.join( self.dir, self.sample_dirs[0], 'mapping', '454MappingProject.xml' )
        m = MappingProject( first_proj )
        
        # Get all the references for the first sample directory(They should all be the same)
        refs = m.get_reference_files()
        return refs

    def _find_reference_by_filename( self, refFileName ):
        """
            Sets the reference to use by filtering all references used by this pipeline run
            down to only the one that contains the text in refFileName

        """
        refs = self._refs_for_proj()

        # Return only the reference that is specified by refFileName
        self.reference_path = filter( lambda x: refFileName.lower() in x.lower(), refs )
        if len( self.reference_path ) > 1:
            raise TooManyReferenceFilesException( refFileName, self.reference_path )
        elif len( self.reference_path ) < 1:
            raise TooFewReferenceFilesException( refFileName, self.reference_path )
        else:
            self.reference_path = self.reference_path[0]

    def _set_reference_info( self ):
        """
            Sets self.reference_info to a SeqIO.to_dict object for the reference
        """
        self.reference_info = SeqIO.to_dict( SeqIO.parse( self.reference_path, 'fasta' ) )

    def parse_description( self, description ):
        """
            Parses a description line for a sequence in a 454AllContigs.fna file

            Arguments:
                description -- String containing the description line of a sequence(Best fetched using Bio.SeqIO)
            Return:
                Dictionary keyed Contig,Ref,Read_Start,Read_End,Length,Num_Reads

            Tests:
            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/', 'california' )
            >>> a = fg.parse_description( "contig00001  FJ969516_PB2_California04, 163..1606  length=1443   numreads=1168" )
            >>> for k in sorted( a.keys() ):
            ...  print "%s => %s" % (k, a[k])
            Contig => contig00001
            Length => 1443
            Num_Reads => 1168
            Read_End => 1606
            Read_Start => 163
            Ref => FJ969516_PB2_California04

            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/', 'california' )
            >>> b = fg.parse_description( "contig00015  Human/1_(PB2)/H5N1/1/Thailand/2004, 2209..2333  length=123   numreads=26" )
            >>> for k in sorted( b.keys() ):
            ...  print "%s => %s" % (k, b[k])
            Contig => contig00015
            Length => 123
            Num_Reads => 26
            Read_End => 2333
            Read_Start => 2209
            Ref => Human/1_(PB2)/H5N1/1/Thailand/2004

            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/', 'california' )
            >>> b = fg.parse_description( "contig00001  CY018765_B/Yamagata/16/1988/4(HA), 523..795  length=273   numreads=2" )
            >>> for k in sorted( b.keys() ):
            ...  print "%s => %s" % (k, b[k])
            Contig => contig00001
            Length => 273
            Num_Reads => 2
            Read_End => 795
            Read_Start => 523
            Ref => CY018765_B/Yamagata/16/1988/4(HA)
       """
        parsed_description = None
        patterns = []
        patterns.append( "(?P<Contig>\w+)\s+(?P<Ref>\S+),\s(?P<Read_Start>\d+)..(?P<Read_End>\d+)\s+length=(?P<Length>\d+)\s+numreads=(?P<Num_Reads>\d+)" )
        for p in patterns:
            parsed_description = re.match( p, description )
            if parsed_description:
                return parsed_description.groupdict()
        raise UnknownIdentifierLineException( "%s. Expected format is \"%s\"" % (description, p) )

    def parse_sample_reads( self, fffAllContigs_records, sample_name ):
        """
            Iterates through the records and gathers the gaps in the file for each contig
            Returns {'Reference1':
                        {'PA': {
                                'Reads': [(Start1,End1),(Start2,End2)...],
                                'Accession': 'SomeAccession',
                                'Ref_Length': '1234'
                               }
                        }
                    }

        """
        for seq in fffAllContigs_records:
            try:
                description_dict = self.parse_description( seq.description )
            except UnknownIdentifierLineException:
                sys.stderr.write( "Skipping %s because it is in an unkown format\n" % seq.description )
                continue

            if self.use_only_ref.lower() not in description_dict['Ref'].lower():
                continue
            try:
                top_key = sample_name
                # If the reference name is not already in the dictionary create a new entry
                if top_key not in self.parsed_reads:
                    self.parsed_reads[top_key] = {}

                # Make sure that the bottom_key name is in the Reference dictionary
                bottom_key = description_dict['Ref']
                if bottom_key not in self.parsed_reads[top_key]:
                    self.parsed_reads[top_key][bottom_key] = {}
                    self.parsed_reads[top_key][bottom_key]['Reads'] = []
                    self.parsed_reads[top_key][bottom_key]['Ref_Length'] = len( self.reference_info[description_dict['Ref']] )

                # add new reads to the e dictionary's reads list
                self.parsed_reads[top_key][bottom_key]['Reads'].append( (int(description_dict['Read_Start']), int(description_dict['Read_End'])) )
            except (NoReferenceFileException,KeyError):
                #sys.stderr.write( "Skipping %s %s because the reference file cannot be found for it\n" % (description_dict['Accession'], top_key) )
                pass
            except TooManyReferenceFilesException:
                sys.stderr.write( "More than one reference file contains the identifier %s. Skipping...\n" % seq.description )

    def _parse_single( self ):
        path = os.path.join( self.dir, 'mapping', '454AllContigs.fna' )
        dir_info = parse_dir_path( self.dir )
        reads = SeqIO.parse( path, 'fasta' )
        self.parse_sample_reads( reads, dir_info['sample'] )

    def _parse_all( self ):
        for d in self.sample_dirs:
            path = os.path.join( self.dir, d, 'mapping', '454AllContigs.fna' )
            dir_info = parse_dir_path( d )
            reads = SeqIO.parse( path, 'fasta' )
            self.parse_sample_reads( reads, dir_info['sample'] )

    def get_uniq_keys( self ):
        """
            Returns a uniq set of Gene names gathered from all the reference files
            
            Tests:
            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/', 'california' )
            >>> fg.get_uniq_keys()
            ['FJ969512_NP_California04', 'FJ969513_MP_California04', 'FJ969514_NS_California04', 'FJ969515_PA_California04', 'FJ969516_PB2_California04', 'FJ969517_NA_California04', 'GQ117044_HA_California04', 'GQ377049_PB1_California04']
        """
        genes = get_idents_from_reference( self.reference_path )
        return sorted( set( genes ) )

    def get_gaps_from_reads( self, read_info, gene ):
        """
            Returns the gaps for a set of reads
            A set of reads will always start from 1 and go to ref_len

            Tests:
            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/', 'california' )
            >>> a = fg.get_gaps_from_reads( {'Reads': [(5,10)], 'Ref_Length': 10}, 'Test' )
            >>> zip( a.keys(), a.values() )
            [('Test_1_Start', 1), ('Test_1_End', 4)]
            >>> a = fg.get_gaps_from_reads( {'Reads': [(1,4),(7,10)], 'Ref_Length': 10}, 'Test' )
            >>> zip( a.keys(), a.values() )
            [('Test_1_Start', 5), ('Test_1_End', 6)]
            >>> a = fg.get_gaps_from_reads( {'Reads': [(1,5),(7,10)], 'Ref_Length': 10}, 'Test' )
            >>> zip( a.keys(), a.values() )
            [('Test_1_Start', 6), ('Test_1_End', 6)]
            >>> a = fg.get_gaps_from_reads( {'Reads': [(1,5)], 'Ref_Length': 10}, 'Test' )
            >>> zip( a.keys(), a.values() )
            [('Test_1_Start', 6), ('Test_1_End', 10)]
            >>> a = fg.get_gaps_from_reads( {'Reads': [(5,7)], 'Ref_Length': 10}, 'Test' )
            >>> zip( a.keys(), a.values() )
            [('Test_1_Start', 1), ('Test_2_Start', 8), ('Test_1_End', 4), ('Test_2_End', 10)]
            >>> a = fg.get_gaps_from_reads( {'Reads': [(1,2),(5,7),(10,12)], 'Ref_Length': 12}, 'Test' )
            >>> zip( a.keys(), a.values() )
            [('Test_1_Start', 3), ('Test_2_Start', 8), ('Test_1_End', 4), ('Test_2_End', 9)]
            >>> fg.get_gaps_from_reads( {'Reads': [(1,10)], 'Ref_Length': 10}, 'Test' )
            {}
        """
        # Extract the reads
        reads = read_info['Reads']
        # Extract the ref_length
        ref_length = read_info['Ref_Length']
        # Gap counter
        count = 1
        # Start the last end read at one since all reads start at 1
        last_end_read = 1

        # Easier if they are sorted.
        reads.sort()

        # The dictionary to fill with gaps
        gaps = {}

        # Get gaps in the beginning or middle
        for read in reads:
            start,end = read
            # Is the start of the read the same as the last end?
            # If not there is a gap before this read(last_end_read -> start-1)
            if last_end_read != start:
                # This part of the name is the same for start and end
                base_name = "%s_%s_" % (gene,count)
                # The start of the gap is where the last end read was set
                gaps[base_name + "Start"] = last_end_read
                # The end of the gap is one less than the beginning of this read
                gaps[base_name + "End"] = start - 1
                # Increment the counter
                count += 1

            # Now set the last_end_read one pas the end of this read
            last_end_read = end + 1

        # Finish up by making sure there are no gaps at the end
        if last_end_read < ref_length:
            base_name = "%s_%s_" % (gene,count)
            gaps[base_name + "Start"] = last_end_read
            gaps[base_name + "End"] = ref_length

        return gaps

    def get_all_gaps( self ):
        """
            Returns a dictionary filled with gap information for each Reference's genes
            {'Reference1Name': {
                                'Gene1': {'Gene1_1_Start': <gap_start>, 'Gene1_1_End': <gap_end>,...},
                                'Gene2': ...
                              }
            }
            Tests:
            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/', 'california' )
            >>> fg.get_all_gaps()['PR_2390']
            {'GQ377049_PB1_California04': {'GQ377049_PB1_California04_1_Start': 1, 'GQ377049_PB1_California04_2_Start': 2272, 'GQ377049_PB1_California04_1_End': 2038, 'GQ377049_PB1_California04_2_End': 2274}, 'GQ117044_HA_California04': {}, 'FJ969515_PA_California04': {}, 'FJ969512_NP_California04': {}, 'FJ969516_PB2_California04': {}, 'FJ969513_MP_California04': {}, 'FJ969514_NS_California04': {'FJ969514_NS_California04_2_Start': 633, 'FJ969514_NS_California04_1_Start': 1, 'FJ969514_NS_California04_2_End': 742, 'FJ969514_NS_California04_1_End': 211}, 'FJ969517_NA_California04': {}}
        """
        if self.single:
            self._parse_single()
        else:
            self._parse_all()
        # Start max_gaps at 1
        max_gaps = 1
        # Get all of the records that are parsed into Reference then gene containers
        records = self.parsed_reads
        # Get the global gene list to iterate over
        genes = self.get_uniq_keys()

        # Will hold gaps for every reference and gene
        gaps = {}

        # Loop through all of the references
        for ref, ref_stuff in records.iteritems():
            # Create the reference key if it does not exist
            if ref not in gaps:
                gaps[ref] = {}
            # Loop through all of the genes that are available and initialize the dictionary to nothing
            for g in genes:
                gaps[ref][g] = {}

            for rg in ref_stuff:
                # Extract the gaps for the current gene, or an empty list if none exist
                reads = ref_stuff.get( rg, {'Reads': [(1,1)], 'Ref_Length': 1} )
                # Get all the caps from the reads
                gaps[ref][rg] = self.get_gaps_from_reads( reads, rg )
                
        return gaps

    def output_csv( self, with_header = True ):
        """
            Outputs csv formatted gaps
            
            Tests:
            >>> fg = FindGaps( '/home/EIDRUdata/Data_seq454/2012_05_11/H3N2/', 'california' )
            >>> fg.output_csv()
        """
        # Filler string
        csv = ""

        # All the gaps
        all_gaps = self.get_all_gaps()
        starts = self.get_gap_headers( all_gaps )
        ends = self.get_gap_headers( all_gaps, 'End' )

        if with_header:
            gap_headers = self.get_csv_gap_headers( starts, ends )
            csv += "Sample,Reference," + gap_headers + "\n"

        # Build the csv file
        # Loop through each reference
        for sample, idents in all_gaps.iteritems():
            csv += "%s,%s," % (sample, self.use_only_ref)
            # Loop through each identifier and put start,end even if it doesn't exist
            for s,e in zip( starts, ends ):
                ref = s.rsplit( '_', 2 )[0]
                if s in idents[ref]:
                    csv += "%s,%s," % (idents[ref][s],idents[ref][e])
                else:
                    csv += "%s,%s," % ('-','-')
            csv += "\n"
        return csv

    def get_csv_gap_headers( self, starts, ends ):
        return ",".join( [",".join( z ) for z in zip( starts, ends )] )

    def get_gap_headers( self, gaps, SE = 'Start' ):
        # Gap headers
        gap_headers = []

        # Get the gap headers
        for ref, genes in gaps.iteritems():
            for gene, gaps in genes.iteritems():
                for g in gaps:
                    if SE in g:
                        gap_headers.append( g )
        gap_headers = list( set( gap_headers ) )
        return sorted( gap_headers )

if __name__ == '__main__':
    import doctest
    doctest.testmod()
