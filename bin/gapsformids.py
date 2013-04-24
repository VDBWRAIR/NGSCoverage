#!/usr/bin/env python

from wrairnaming import Formatter

from argparse import ArgumentParser
import multiprocessing
import subprocess
import sys
import os
import os.path
import glob

# I know you are supposed to avoid doing global shared
# but i want quick and dirty
mutex1 = multiprocessing.Lock()
mutex2 = multiprocessing.Lock()

# Yep, a global
outputdir = None

def findgap( indexline ):
    lineno, projdir = indexline

    if not projdir:
        return

    runaligncov( projdir, lineno )

def runaligncov( projdir, lineno ):
    global mutex1, mutex2, outputdir

    try:
        f = Formatter()
        sample = f.GsProject.directory_format.parse_input_name( projdir )['sample']
    except Exception:
        return

    mutex2.acquire()
    sys.stdout.write( "Processing directory #%s %s\n" % (lineno, projdir) )
    mutex2.release()
    out, err = subprocess.Popen( ['aligncoverage', '-d', projdir, '--csv'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE ).communicate()
    # Don't write empty gaps file
    if out:
        mutex1.acquire()
        fh = open( os.path.join( outputdir, sample.replace( ' ', '' ) + '.gap' ), 'a' )
        fh.write( projdir + "\n" )
        fh.write( out )
        fh.close()
        mutex1.release()

def readindex( indexfile ):
    try:
        if type( indexfile ) == str:
            fh = open( indexfile )
        elif type( indexfile ) == file:
            fh = indexfile
        else:
            print "I have no idea how to read given index input from %s" % indexfile
            sys.exit( -1 )
        linecount = 0
        for line in fh:
            yield (linecount, line.strip())
            linecount += 1
    except IOError:
        print "Failed to open %s" % indexfile
        sys.exit( -1 )

def cleanoutputdir( outputdir ):
    gapfiles = glob.glob( '%s/*.gap' % outputdir )
    for gapfile in gapfiles:
        try:
            os.unlink( gapfile )
        except IOError as e:
            sys.stderr.write( "Could not remove %s" % gapfile )
            raise e

def ops():
    global outputdir

    parser = ArgumentParser()

    cpus = multiprocessing.cpu_count() - 1
    parser.add_argument( '-i', '--index', dest='index', default=sys.stdin, help='GsMapper file list to use. Defaults to standard input' )
    parser.add_argument( '--cpus', dest='cpus', default=cpus, type=int, help='Number of cpus to use[Default: %s]' % cpus )
    parser.add_argument( '-o', '--output', dest='outputdir', default=None, help='Output directory[Default: Current directory]' )

    args = parser.parse_args()

    if args.outputdir is None:
        outputdir = './'
    else:
        outputdir = args.outputdir

    try:
        outputdir = os.path.abspath( outputdir )
    except Exception:
        print "I can't convert %s to an absolute path" % outputdir
        
    try:
        os.makedirs( outputdir )
    except OSError:
        cleanoutputdir( outputdir )

    return args

def main( ops ):
    global outputdir

    dirs = readindex( ops.index )
    pool = multiprocessing.Pool( ops.cpus )
    pool.map( findgap, dirs )

if __name__ == '__main__':
    main( ops() )
