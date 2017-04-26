#!/usr/bin/env bash

#
#	Set parameters below to fit your system
#



BED=
GENOME=
OUTDIR=

function usage()
{
cat <<EOF
usage: $0 -b <bedfile> -g <genome> -o <outdir> 
EOF
exit 1;
}

while getopts b:g:o: opt
do
case ${opt} in
b) BED=${OPTARG};;
g) GENOME=${OPTARG};;
o) OUTDIR=${OPTARG};;
*) usage;;
esac
done

if [ "${BED}" = "" ]; then usage; fi
if [ "${GENOME}" = "" ]; then usage; fi

foo=$BED
tmpfile=${foo##*/}

# echo $tmpfile

#array=(0.01 0.02 0.03 0.04 0.05 0.06 0.07)
array=(0.02)

len=${#array[*]}

for (( i=18; i <= 24; i+=1 )); do
    j=0
    
    while [ $j -lt $len ]
    do
	outname=$OUTDIR"/"$tmpfile"_LEN"$i"_ERROR"$j".fa"
	printf "%s\n" $outname;
	#  echo $i $j ${array[$j]} $BED $GENOME   $BED.simldshfksj ;
	bedsim sim  $GENOME $BED -n 100000 --read-len $i --start-error ${array[$j]}   --stop-error ${array[$j]} -o $outname
	
	let j++
    done	
done


