# CASP related tools

## Motivations
Improve protein structure prediction accuracy by MSA searching strategies

## FAQ
1. How to build the docker image for [Alphafold2](http://github.com/bigict/alphafold)
   ```shell
   $ git clone git@github.com:bigict/CASP.git
   $ cd CASP
   $ git submodule update
   $ cd alphafold
   $ docker build -t alphafold -f docker/Dockerfile .
   ```
2. How to run [Alphafold2](http://github.com/bigict/alphafold)
   ```
   sh run_alphafold.sh -d <db_dir> -o <output_dir> <fasta_file_1> <fasta_file_2> ...
   ```
