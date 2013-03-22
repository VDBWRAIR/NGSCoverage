thisdir=$(cd $(dirname $0) && pwd)
primerdir=$(dirname ${thisdir})/Primer
bindir=$(dirname $(dirname $(dirname ${thisdir})))/bin

pushd ${thisdir}

pushd ../05_11_2012_1_TI-MID51_PR_2305_pH1N1
echo "Generating mid51withprimer.png"
${bindir}/gapstoscatter --csv gaps.csv -p ${primerdir}/sH1N1.fasta -o ${thisdir}/mid51withprimer.png
echo "Generating mid51withoutprimer.png"
${bindir}/gapstoscatter --csv gaps.csv -o ${thisdir}/mid51withoutprimer.png
popd

pushd ../08_06_2012_1_Ti-MID30_D84_140_Dengue3
echo "Generating mid30withprimer.png"
${bindir}/gapstoscatter --csv gaps.csv -p ${primerdir}/Den3.fasta -o ${thisdir}/mid30withprimer.png
echo "Generating mid30withoutprimer.png"
${bindir}/gapstoscatter --csv gaps.csv -o ${thisdir}/mid30withoutprimer.png
popd

echo pb1_allh3n22012.csv
echo "Generating allh3n2pb1withprimer.png"
${bindir}/gapstoscatter --csv ../pb1_allh3n22012.csv -p ${primerdir}/H3N2.fasta -o ${thisdir}/allh3n2pb1withprimer.png
echo "Generating allh3n2pb1withoutprimer.png"
${bindir}/gapstoscatter --csv ../pb1_allh3n22012.csv -o ${thisdir}/allh3n2pb1withoutprimer.png
        
echo H3N2__Indiana__2011__NS.png
${bindir}/gapstoscatter --csv ../H3N2__Indiana__2011__NS.gaps -o ${thisdir}/H3N2__Indiana__2011__NS.png

echo SingleGapsAtEnds.png
${bindir}/gapstoscatter --csv ../singlegaplcatends.gaps -o ${thisdir}/SingleGapsAtEnds.png

popd
