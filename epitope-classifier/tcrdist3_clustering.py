import pandas as pd
from os.path import join as opj
from tcrdist.mappers import vdjdb_to_tcrdist2, vdjdb_to_tcrdist2_mapping_TRA, vdjdb_to_tcrdist2_mapping_TRB
import json
from tcrdist.repertoire import TCRrep
from tcrdist.adpt_funcs import get_basic_centroids
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, help='Path to input .tsv file')
parser.add_argument('-o', '--output', type=str, help='Path to output file')
parser.add_argument('--chain', type=str, help='Chain to be used in the analysis. Can be one of "alpha" or "beta".')

args = parser.parse_args()
path = args.input
output_path = args.output
chain = args.chain

raw = pd.read_table(path, sep = '\t')
trb = TCRrep(cell_df = raw, organism = 'human', chains = [chain]
trb = get_basic_centroids(trb, max_dist = 200)

# write output
trb.centroids_df.to_csv(output_path, sep = '\t')
