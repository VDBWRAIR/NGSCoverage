NGS Coverage
=======
Find all gaps for a Roche 454 assembly or mapping project
Eventually, will accept BAM/SAM file too

Example Images
==============

Many Samples of the Same type with the Primer graphed at the bottom
-------------------------------------------------------------------

![Many Samples With Primer](https://github.com/VDBWRAIR/NGSCoverage/raw/master/coverage/Examples/testgraphs/allh3n2pb1withprimer.png)

Single Sample With primer
-------------------------
![Single Sample With Primer](https://github.com/VDBWRAIR/NGSCoverage/raw/master/coverage/Examples/testgraphs/mid51withprimer.png)

Prerequisites
=============
Please check under setup.py inside of install_requires for details on what python packages are required.
https://github.com/VDBWRAIR/NGSCoverage/blob/master/setup.py#L33

aligncoverage
===================

aligncoverage utilizes the 454AlignmentInfo.tsv file to build a CoverageRegion for each contig that is located in the file.
The way it determines a CoverageRegion's type is base on the depth of each base and is not easily configurable at this time.

See: [Configurable LowCoverage](https://github.com/VDBWRAIR/NGSCoverage/issues/1)

If the project that is being analyzed is a Assembly then the AlignDepth is used otherwise the TotalDepth column is used. You can refer to the Roche documentation on all
of the available fields.

Since the 454AlignmentInfo.tsv file contains multiple contigs that coorespond to the same mapped reference(in a mapping project), aligncoverage merges all of the regions
produced from all SeqAlignments for every reference.

It then outputs those regions in human readable or csv format

Usage
-----
```
#> aligncoverage -h
usage: aligncoverage [-h] -d DIR [-c]

optional arguments:
  -h, --help         show this help message and exit
  -d DIR, --dir DIR  454 Project directory path to find gaps in
  -c, --csv          Output in csv format
```

Usage is fairly straight forward and simple
You have to provide it a Roche 454 project directory using the -d  or --dir option
If you want csv output then also specify --csv or -c
--csv is useful for when you use gapstoscatter

Note: All output is written to the terminal(Standard Out) so you will likely want to redirect it using the Unix redirect operators(>)

gapstoscatter
=============

gapstoscatter will read in any csv file generated by aligncoverage and produce a graphic that displays the gaps and low coverage visually
If you give it a primer file it will map the primers to visually compare with.

Usage
-----
```
#> gapstoscatter -h
usage: gapstoscatter [-h] --csv CSVFILE [-o OUTPUTFILE] [-t TITLE]

optional arguments:
  -h, --help            show this help message and exit
  --csv CSVFILE         CSV Gaps file to parse
  -o OUTPUTFILE         Filepath to put output image[Default: ./gaps.png]
  -t TITLE, --title TITLE
                        Title for the scatterplot
```

Full examples
------------

Generate Graphic without Primer
```
#> aligncoverage -d /some/454/project/dir --csv > gaps.csv
#> gapstoscatter --csv gaps.csv
```

Generate Graphic with Primer
```
#> aligncoverage -d /some/454/project/dir --csv > gaps.csv
#> gapstoscatter --csv gaps.csv -p /path/to/primer
```

Notes: 
- You can also change the title of the graphic by utilizing the -t option. If you don't specify the -t it will use the name of the file
  from the --csv option
- If you do not specify the -o option it will use gaps.png as the output

Primer File
-----------
Some quick notes about the structure expected in the primer file
To know where the beginning and end of a primer are and what direction each primer identifier
needs to be specifically formatted.

Basically you need to ensure that the sequence id field looks similar to this:
R followed by a number
-- or --
F followed by a number

Where R represents reverse and F represents forward primer
The number represents the start position of the primer if it is forward or the end if it is reverse

The code that parses this is located here:
https://github.com/VDBWRAIR/pyWrairLib/blob/master/wrairlib/primer.py#L358

gapsformids.py
==============
This script is useful if you want to run a bunch of existing newbler project folders through the aligncoverage script.
Allows you to generate graphics easily for many projects so you can see how well your coverage has been.

The index file is simply a text file containing newbler project directories listed one per line
This script has multiprocessor support to speed it up.

Usage
-----
```

usage: gapsformids.py [-h] [-i INDEX] [--cpus CPUS] [-o OUTPUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        Index to use. Defaults to standard
                        input
  --cpus CPUS           Number of cpus to use[Default: 1]
  -o OUTPUTDIR, --output OUTPUTDIR
                        Output directory[Default: Current directory]
```


findgap(deprecated)
=============
```
#>findgap.py -d <directory> -r <reference_name>
```

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
```
-r infB_Victoria.fasta
-r infb_victoria.fasta
-r victoria
-r infb
```

As you can see the -r option is very flexible and matches as long as the value you give it is uniq among the
reference files
!! Be Careful though !!
-r H1N1
This would match both pdmH1N1_California.fasta as well as H1N1_boston.fasta

Example Runs
------------
```
#> findgap.py -d /some/project/dir/2012_05_11/H3N2 -r managua
#> findgap.py -d /some/project/dir/2012_06_22/FluB_Victoria -r bangladesh
#> findgap.py -d /some/project/dir/2012_07_19/R1R2_Den1_Ref_D1_IPCBIDV3786_2008_cambodia -r cambodia
```

Output
------
As with any Linux script you can redirect the output( using the > operator ) to a file which you will probably want to do
```
findgap.py -d /some/project/dir/2012_05_11/H3N2 -r managua > /some/project/dir/2012_05_11/H3N2/gaps.csv
findgap.py -d /some/project/dir/2012_06_22/FluB_Victoria -r bangladesh > /some/project/dir/2012_06_22/FluB_Victoria/gaps.csv
findgap.py -d /some/project/dir/2012_07_19/R1R2_Den1_Ref_D1_IPCBIDV3786_2008_cambodia -r cambodia > /some/project/dir/2012_07_19/R1R2_Den1_Ref_D1_IPCBIDV3786_2008_cambodia/gaps.csv
```

Opening the Output
------------------
As this script generates comma separated value files you can open the file with Excel or OpenOffice.org and then select comma as the delimiter

Hint: You can open the file from the command line in Linux by using openoffice.org -calc gaps.csv
