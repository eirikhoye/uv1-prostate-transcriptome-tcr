#import findspark
#findspark.init()
import pyspark
import os
import sys
import imnet
import numpy as np
import pandas as pd
import networkx as nx
import pickle
from pyspark.sql import SparkSession


"""
OLD DEPRECATED
# Usage: python scripts/run_imnet.py {path/to/input.tsv} {name_outfile} OLD DEPRECATED


Usage: opt/spark/bin/spark-submit --master local[*] {path/to/input.tsv} {name_outfile}
i.e.
Usage: opt/spark/bin/spark-submit --master local[*] workflow/scripts/run_imnet.py tcrseqs_for_imnet/expanded_clones_raw_seqs.txt tcrseqs_graphs/expanded_clones
"""

os.environ['SPARK_HOME'] = '/storage/mathelierarea/processed/eirikhoy/tcr_uv1/opt/spark/'
#sc = pyspark.SparkContext('local[*]')

spark = SparkSession.builder.getOrCreate()
sc = spark.sparkContext


def run_imnet(file_path, out_path, min_ld=1, max_ld=1):
    # Read rearrangement.tsv and pull amino acid sequences
    df = pd.read_csv(file_path, sep=',')
    strings = list(df.amino_acid)
    
    # generate graph from amino acid strings using pyspark
    g = imnet.process_strings.generate_graph(strings, sc=sc, min_ld=1, max_ld=1)
    g.remove_edges_from(nx.selfloop_edges(g))   
    param_glob = {}
    
    # Write graph in a graphml file format
    with open(out_path, 'wb') as f:
        pickle.dump(g, f)

# Input needs to be Adaptive rearrangement style .tsv file, most important is that column containing clonal amino acid strings is named amino_acid
run_imnet(sys.argv[1], sys.argv[2], min_ld=1, max_ld=1)

sc.stop()