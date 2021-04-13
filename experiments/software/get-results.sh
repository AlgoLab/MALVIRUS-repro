#!/bin/bash

for catalog in */
do
    pushd $catalog
    for sample in */
    do
        csvjoin -I \
            <(echo -e "catalog,sample\n${catalog%/},${sample%/}") \
            <(csvcut -c 2,3,5,6 $sample/malva.lineage.csv) \
            <(python ../../../software/parse-malva-csv.py malva_ < $sample/malva.csv) \
            <(python ../../../software/parse-time-output.py malva_ < $sample/malva.time) \
            > $sample/summary.csv
    done
    popd
done
