#!/usr/bin/env bash
# Download BAMs from the Genexus machine
# Version 1.0.1 (2022-05-04) (Niema Moshiri)

# important constants
GENEXUS_USERNAME='ionservice'
GENEXUS_IP='132.239.149.188'
GENEXUS_REPORTS='/serverdata/results/analysis/output/reports'
GENEXUS_BAM='merged.bam.ptrim.bam'
OUT_BAM_SUF='trimmed.bam'

# check usage
if [ "$#" -ne 2 ] ; then
    echo "USAGE: $0 <list_of_EXC_names> <output_folder>"; exit 1
fi

# check list of EXC names
if [ -f "$1" ] ; then
    echo "Using list of EXC names: $1"
else
    echo "List of EXC names not found: $1" ; exit 1
fi
for exc in $(cat "$1") ; do
    if [[ ("$exc" != EXC_*_* || "$exc" == EXC_*_*_*) && ("$exc" != CALM_*_*_*) ]] ; then
        echo "Invalid EXC/CALM name: $exc"
        echo "Must follow EXC_???_?????? (e.g. EXC_MW5_503759) or CALM_???_??????_?? (e.g. CALM_SEP_008217_16)"
        exit 1
    fi
done

# check output folder
if [ -d "$2" ] ; then
    echo "Using existing output folder: $2"
elif [ -f "$2" ] ; then
    echo "Specified output folder already exists as a file: $2" ; exit 1
else
    echo "Creating output folder: $2"
    mkdir "$2"
fi

# download BAMs
for exc in $(cat "$1") ; do
    echo "Downloading: $exc"
    scp "$GENEXUS_USERNAME@$GENEXUS_IP:$GENEXUS_REPORTS/*$exc*/$GENEXUS_BAM" "$2/$exc.$OUT_BAM_SUF"
done
