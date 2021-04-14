#!/bin/bash

while read sid
do
  echo Downloading sample $sid
  rm -f $sid.fastq
  fasterq-dump -p --skip-technical --split-files $sid
  cat ${sid}_1.fastq ${sid}_2.fastq >> $sid.fastq
  rm ${sid}_1.fastq ${sid}_2.fastq
  gzip -9v $sid.fastq &
done < selected_SRA_IDs.txt
