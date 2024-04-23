"""Clustering MSA based on amino acid mutation site
"""
import os
import re
import pickle

import numpy as np

from absl import app
from absl import flags
from absl import logging

from alphafold.common import residue_constants

flags.DEFINE_string(
    'pkl_file', None, 'Paths to FASTA files, each containing a prediction '
    'target that will be folded one after another. If a FASTA file contains '
    'multiple sequences, then it will be folded as a multimer. Paths should be '
    'separated by commas. All FASTA paths must have a unique basename as the '
    'basename is used to name the output directories for each prediction.')
flags.DEFINE_string('output_dir', None, 'Path to a directory that will '
                    'store the results.')

FLAGS = flags.FLAGS

def to_str_msa(int_msa):
  return ''.join(residue_constants.ID_TO_HHBLITS_AA[res] for res in int_msa)

def main(argv):
  if len(argv) < 2:
    raise app.UsageError('Too few command-line arguments.')

  with open(FLAGS.pkl_file, 'rb') as f:
    feature_dict = pickle.load(f)
  assert 'msa' in feature_dict

  str_msa = list(map(to_str_msa, feature_dict['msa']))

  filter_rules = [re.compile(arg, re.I) for arg in argv[1:]]

  new_order = [0]  # Always contain MSA_0
  for i, msa in enumerate(str_msa[1:], start=1):
    for j, filter_rule in enumerate(filter_rules, start=1):
      if filter_rule.match(msa):
        new_order.append(i)
        logging.info('Select MSA [%.5d] %s %s', i, str_msa[i], argv[j])
        break

  # filter msa based on new_order
  new_order = np.array(new_order)
  # 'deletion_matrix_int', 'msa', 'num_alignments', 'msa_species_identifiers'
  for key, value in feature_dict.items():
    if isinstance(value, np.ndarray):
      print(key, value.shape)
  for key in ('deletion_matrix_int', 'msa', 'msa_species_identifiers'):
    feature_dict[key] = np.take(feature_dict[key], new_order, axis=0)
  print(feature_dict['num_alignments'])
  feature_dict['num_alignments'].fill(len(new_order))
  print(feature_dict['num_alignments'])
  for key, value in feature_dict.items():
    if isinstance(value, np.ndarray):
      print(key, value.shape)

  # balabala
  os.makedirs(FLAGS.output_dir, exist_ok=True)
  with open(os.path.join(FLAGS.output_dir, 'features.pkl'), 'wb') as f:
    pickle.dump(feature_dict, f, protocol=4)

if __name__ == '__main__':
  flags.mark_flags_as_required([
      'output_dir',
      'pkl_file',
  ])
  app.run(main)
