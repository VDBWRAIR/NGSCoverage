FindGap
=======
Find all gaps for a Roche 454 assembly or mapping project

Prerequisites
=============
- pyWrairLib
https://github.com/VDBWRAIR/pyWrairLib

aligncoverage Usage
===================

aligncoverage utilizes the 454AlignmentInfo.tsv file to build a CoverageRegion for each contig that is located in the file.
The way it determines a CoverageRegion's type is base on the depth of each base and is not easily configurable at this time. #TODO#
If the project that is being analyzed is a Assembly then the AlignDepth is used otherwise the TotalDepth column is used. You can refer to the Roche documentation on all
of the available fields.

Since the 454AlignmentInfo.tsv file contains multiple contigs that coorespond to the same mapped reference(in a mapping project), aligncoverage merges all of the regions
produced from all SeqAlignments for every reference.

It then outputs those regions in human readable or csv format

Usage
=====
#> aligncoverage -h
usage: aligncoverage [-h] -d DIR [-c]

optional arguments:
  -h, --help         show this help message and exit
  -d DIR, --dir DIR  454 Project directory path to find gaps in
  -c, --csv          Output in csv format

Usage is fairly straight forward and simple
You have to provide it a Roche 454 project directory using the -d  or --dir option
If you want csv output then also specify --csv or -c


findgap Usage(deprecated)
=============
#>findgap.py -d <directory> -r <reference_name>

<reference_name> should be the reference that you want to find gaps for.
You don't need to give it the full name of the reference file but you can.

Reference is located using the following method:
1. Read <midkey sample directory>/mapping/454MappingProject.xml
2. Extract all Reference files located between <ReferenceFiles> and </ReferenceFiles>
3. Select only the reference that contains the reference you provide using the -r option

Example:
If you have the following references:
/some/project/dir/2012_05_11/H3N2/Ref/infB_Victoria.fasta
/some/project/dir/2012_05_11/H3N2/Ref/pdmH1N1_California.fasta
/some/project/dir/2012_05_11/H3N2/Ref/H3N2_Managua.fasta
/some/project/dir/2012_05_11/H3N2/Ref/H5N1_Thailand.fasta
/some/project/dir/2012_05_11/H3N2/Ref/H1N1_boston.fasta

The following would match the infB_Victoria.fasta reference
-r infB_Victoria.fasta
-r infb_victoria.fasta
-r victoria
-r infb

As you can see the -r option is very flexible and matches as long as the value you give it is uniq among the
reference files
!! Be Careful though !!
-r H1N1
This would match both pdmH1N1_California.fasta as well as H1N1_boston.fasta

Example Runs
============
#> findgap.py -d /some/project/dir/2012_05_11/H3N2 -r managua
#> findgap.py -d /some/project/dir/2012_06_22/FluB_Victoria -r bangladesh
#> findgap.py -d /some/project/dir/2012_07_19/R1R2_Den1_Ref_D1_IPCBIDV3786_2008_cambodia -r cambodia

Output
======
As with any Linux script you can redirect the output( using the > operator ) to a file which you will probably want to do
findgap.py -d /some/project/dir/2012_05_11/H3N2 -r managua > /some/project/dir/2012_05_11/H3N2/gaps.csv
findgap.py -d /some/project/dir/2012_06_22/FluB_Victoria -r bangladesh > /some/project/dir/2012_06_22/FluB_Victoria/gaps.csv
findgap.py -d /some/project/dir/2012_07_19/R1R2_Den1_Ref_D1_IPCBIDV3786_2008_cambodia -r cambodia > /some/project/dir/2012_07_19/R1R2_Den1_Ref_D1_IPCBIDV3786_2008_cambodia/gaps.csv

Opening the Output
==================
As this script generates comma separated value files you can open the file with Excel or OpenOffice.org and then select comma as the delimiter

Hint: You can open the file from the command line in Linux by using openoffice.org -calc gaps.csv
