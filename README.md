### Setup
```
git clone https://github.com/julianzaugg/viral_annotate
cd viral_annotate
conda create -n viral_annotate -c bioconda snakemake conda=4.8.4 -y
conda activate viral_annotate
```

### Usage
Example:
```commandline
FASTA_FILE=my_viral_sequences.fasta
OUT_DIR=my/output/dir
THREADS=4

# (optional) Location to install conda envs
CONDA_PREFIX=/srv/home/$(whoami)/.conda/envs 

mkdir -p $OUT_DIR
cp config.yaml $OUT_DIR

snakemake --keep-going \
--latency-wait 20 \
--cores $THREADS \
--use-conda \
--conda-prefix $CONDA_PREFIX \
run_pipeline \
--directory $OUT_DIR \
--printshellcmds \
--configfile $OUT_DIR/config.yaml \
--snakefile $PWD/Snakefile \
--config \
VIRAL_INPUT=$FASTA_FILE \
MAX_THREADS=$THREADS 
```
