#!/usr/bin/env bash
#
#	Set parameters below to fit your system
#

BWA_DIR="bwa"

function run_aln () {
    printf "Processing %s (bwa)" $file;
    

    $BWA_DIR aln -n 2 -t 4  $GENOME_DIR   $file  | bwa samse $GENOME_DIR   - $file  -n 50   | samtools sort -@ 4 -O bam -l 9 -T ~/tmp - >  $file.bwa.bam

    if [ "$?" -eq "0" ]
    then
	printf "%30s\n" "Success";
    else
	printf "\nERROR: bwa FAILED!\n\n";
	exit 1;
    fi
}

GENOME=
INDIR=

function usage()
{
    cat <<EOF
usage: $0  -g <genome file> -i <input dir>
EOF
    exit 1;
}

while getopts g:i: opt
do
    case ${opt} in
	g) GENOME_DIR=${OPTARG};;
	i) INDIR=${OPTARG};;
	*) usage;;
    esac
done

if [ "${GENOME_DIR}" = "" ]; then usage; fi
if [ "${INDIR}" = "" ]; then usage; fi

for file in $INDIR/*
do
    if [ -f $file ]; then 
	if [[ $file =~ .f[q,a]$ ]]; then	    
	    if [ -f $file.bwa.bam ]; then
		echo "$file already processed"
	    else
		run_aln $file
	    fi
	fi
    fi
done
