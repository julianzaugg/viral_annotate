from snakemake.utils import logger, min_version
min_version("6.0")

import sys
import os
import re
import tempfile
import glob
import subprocess
import pandas as pd

sys.path.append(
    os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), "scripts")
)

configfile: "config.yaml"

INPUT_FASTA = config["VIRAL_INPUT"]

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

onstart:
    if INPUT_FASTA == "none":
        sys.exit("Need to specify input FASTA file")

ABSOLUTE_DATA_PATH = os.getcwd()
SNAKE_PATH= workflow.basedir

include: "rules/viral_annotate.smk"

rule run_pipeline:
    input:
        "data/viral_summary/viral_summary.tsv",
        "finished_viral_annotation",
        "finished_viral_lineage",
    output:
        touch("finished_pipeline")
