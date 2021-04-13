#!/bin/bash

while read catalog
do
    echo $(date) Starting $catalog
    snakemake -pr -j21 --config vcfname=$catalog basepath=$PWD -d output/$catalog
done
