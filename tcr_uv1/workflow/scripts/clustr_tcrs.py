import pandas as pd
from clustcr import Clustering
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--clusters", required=True)
parser.add_argument("--motifs", required=True)
parser.add_argument("--aa_colname", required=True)
parser.add_argument("--threads", type=int, default=1)
args = parser.parse_args()

data = pd.read_csv(args.input)[args.aa_colname]
clustering = Clustering(n_cpus=args.threads)
result = clustering.fit(data)
result.write_to_csv(path=args.clusters)
result.summary().to_csv(args.motifs)
