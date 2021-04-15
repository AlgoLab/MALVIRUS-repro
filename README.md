# MALVIRUS - Reproducible experiments

This repository contains all the workflows and the scripts needed to reproduce the experiments about MALVIRUS accuracy.

The only software prerequisite is the open source package manager [conda](https://conda.io/en/), needed to retrieve and install all the other software used in the experiments.
Notice that the experiments have been run on Ubuntu 20.04 and we did not test them on different OSs (although MALVIRUS is cross-platform).

After you cloned (or downloaded this repository, in order to install all the software and retrieve the read sample data from SRA you should use the following command:

``` shell
bash ./prepare-environment.sh
```

Since the terms of use of GISAID do not allow us to redistribute their data, you should also download them from GISAID and copy them into the repository.
In particular, you have to download the package containing all the sequences in FASTA format (the filename should be `sequences_fasta_<date>.tar.xz`), you should decompress it and extract the file `sequences.fasta`. Then you have to recompress in XZ format and place it in the directory `catalogs/input/gisaid/<date>/` with the name `sequences.fasta.xz`.
Similarly, download the metadata package (the filename should be `metadata_tsv_<date>.tar.xz`), extract the file `metadata.tsv`, recompress it in GZIP format and place the resulting file in the same directory.
Then, update the file `catalogs/Snakefile` at line 5 by replacing `"210407"` with the date indicated in the files you previously downloaded.

Finally, download the sequences whose accession IDs are listed in file [`gisaid_accession_ids.txt`](./gisaid_accession_ids.txt) and place each one of them in a file `<GISAID_ACCESSION_ID>.fasta` of directory `experiments/input/gisaid`. Replace all the FASTA headers with `>sample`.

After all the software and the data have been retrieved, you can run the workflow contained in directory `catalogs` with the command (replace `4` with the number of CPU cores you want to use for the execution):

``` shell
snakemake --use-conda -j4
```

When all the catalogs have been built, you can test MALVIRUS accuracy on the read samples using the workflow contained in directory `experiments`. You should use the following command to start the workflow on every catalog (from directory `experiments`):

``` shell
ls vcf | ./multi-exec.sh
```

At the end, you will have a set of CSV files `experiments/output/<catalog name>/summary.k31.csv` summarizing the results for each catalog.
