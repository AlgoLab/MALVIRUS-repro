#!/bin/bash

test -f env-prepared && exit 0

echo Preparing the environment. It is a long process.

conda env create -f ./environment.yml

pushd ./experiments/input/samples

conda run -n malvirus-repro-exp ./download-samples.sh

popd

pushd ./experiments/software/malva

conda run -n malvirus-repro-exp make

popd

pushd ./catalogs/software

conda run -n malvirus-repro-exp make

popd

echo "Please activate environment malvirus-repro-exp with 'conda activate malvirus-repro-exp' to proceed with the experiments"

echo env-prepared
