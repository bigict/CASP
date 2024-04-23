# CASP related tools

## Motivations
Improve protein structure prediction accuracy by MSA searching strategies

## FAQ
1. How to build the docker image for [Alphafold2](http://github.com/bigict/alphafold)
   ```shell
   $ git clone git@github.com:bigict/CASP.git
   $ cd CASP
   $ git submodule update --init --recursive --remote
   $ cd alphafold
   $ docker build -t alphafold -f docker/Dockerfile .
   ```
2. How to run [Alphafold2](http://github.com/bigict/alphafold)
   ```
   $ pip install -r requirements.txt
   $ sh run_alphafold.sh -d <db_dir> -o <output_dir> <fasta_file_1> <fasta_file_2> ...
   ```
3. How to run msa searching tool
   ```
   $ cd CASP
   $ sh run_msa_tool.sh -d <db_dir> -o . examples/T1024-D1.fasta  # db_dir - BFD database dir
   ```
4. How to run the whole pipeline
   ```
   $ sh run_casp.sh -d <db_dir> -o . examples/T1024-D1.fasta
   ```
