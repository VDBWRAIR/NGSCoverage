FindGap
=======
Find all gaps for a 454 Pipeline run

How to read this file
=====================
Any line that begins with #> is something you run in the command line
Any line that begins with > is an continuation line from a #> line
All other lines are just text lines that are there as reference

Prerequisites
=============
- pyWrairLib
https://github.com/VDBWRAIR/pyWrairLib

Usage
=====
#>findgap.py -d <directory> -r <reference_name>
<directory> should be a pipeline run directory such as any of the below:
/some/project/dir/2012_07_19/R1R2_Den1_Ref_D1_IPCBIDV3786_2008_cambodia
/some/project/dir/2012_06_22/FluB_Victoria
/some/project/dir/2012_05_11/H3N2

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

Shortcuts
=========
If you don't want to type the full path to the 454 Pipeline run directory so much you can use the `pwd` shortcut
`pwd` will insert the current path wherever you put it

#> cd /some/project/dir/2012_05_11/H3N2 
#> findgap.py -d `pwd` -r managua > `pwd`/gaps.csv

#> cd /some/project/dir/2012_06_22/FluB_Victoria 
#> findgap.py -d `pwd` -r bangladesh > `pwd`/gaps.csv

#> cd /some/project/dir/2012_07_19/R1R2_Den1_Ref_D1_IPCBIDV3786_2008_cambodia 
#> findgap.py -d `pwd` -r cambodia > `pwd`/gaps.csv

Opening the Output
==================
As this script generates comma separated value files you can open the file with Excel or OpenOffice.org and then select comma as the delimiter

Hint: You can open the file from the command line in Linux by using openoffice.org -calc gaps.csv
