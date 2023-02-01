# Kintun-GIP
:warning: This project is currently under development :warning:

## Dependencies
- tRNAscan-SE 2.0.11 
- barrnap 0.9
- biopython1.80
- glob2 0.7
- mmseqs 14.7e284
- sibeliaz 1.2.5

## Recommended installation 
We have been successful installing the dependencies in a conda environment and running the code in that environment.
```
conda create -n gmi-kintun -c bioconda -c conda-forge "python>=3.10" aragorn barrnap biopython glob2 mmseqs2 sibeliaz pandas
conda activate gmi-kintun
git clone https://github.com/GMI-Lab/Kintun-GIP.git
cd Kintun-GIP/
pip install .
```

## Kintun-GIP clust
This tool will annotate and cluster all the t(m)DNAs that are coded in a collection of strains. It is required 
that all chromosomes of the strains are fully assembled, i.e. each chromosome is in only one contig. It is mandatory 
to choose a prefix for a new nomenclature scheme to be resolved. This avoids confusing different schemes.

```
kintunGIP clust [-h] -i FASTA_DIR -x FILE_EXTENSION -o RESULTS_DIR -p PREFIX [-t THREADS] [-log LOG_FILE] [-v]

options:
  -h, --help         show this help message and exit
  -i FASTA_DIR       Directory with fasta files (REQUIRED)
  -x FILE_EXTENSION  Extension for fasta files (REQUIRED;default: fasta)
  -o RESULTS_DIR     Name for results directory (REQUIRED)
  -p PREFIX          Prefix for the t(m)DNA clustering scheme (REQUIRED)
  -t THREADS         Number of threads for nucmer and SibeliaZ
  -log LOG_FILE      Name of log file (Optional)
  -v                 show program's version number and exit
```

