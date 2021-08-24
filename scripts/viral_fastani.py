import os
import subprocess
from Bio import SeqIO
from pathlib import Path
import glob


split_input_dir = "data/viral_clustering/fastani/fasta_inputs"
fastani_out_dir = "data/viral_clustering/fastani"


try:
    os.makedirs(split_input_dir)
except FileExistsError:
    pass


def split_fasta(input_fasta_filename, output_directory):
    with open(input_fasta_filename, 'r') as fh:
        sequences = SeqIO.parse(fh, "fasta")
        for s in sequences:
            outname = "{}.fasta".format(s.name)
            with open(Path(output_directory) / outname, 'w') as fh2:
                SeqIO.write(s, fh2, format="fasta-2line")

# Split viral FASTA
# split_fasta(snakemake.input.checkv_selected, split_input_dir)
split_fasta(snakemake.input.all_viral_sequences, split_input_dir)

# Get list of fasta files and write to file
input_viral_paths = [os.path.abspath(filename) for filename in glob.iglob(f"{split_input_dir}/*.fasta")]
with open(f"{fastani_out_dir}/viral_genomes.txt", 'w') as fh:
    fh.write("\n".join(input_viral_paths))

# Run FastANI
subprocess.Popen(
    f"""
    fastANI \
    --ql {fastani_out_dir}/viral_genomes.txt \
    --rl {fastani_out_dir}/viral_genomes.txt \
    --fragLen {snakemake.params.frag_len} \
    --minFraction {snakemake.params.min_fraction} \
    -t {snakemake.threads} \
    -o {fastani_out_dir}/fastani_viral.tsv
    """,
    shell=True).wait()
