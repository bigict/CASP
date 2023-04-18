"""Tools for split fasta to single sequence file, run
     ```bash
     $python fasta_split.py -h
     ```
     for further help.
"""
import os
import functools
import glob
import multiprocessing as mp
import re
import logging

from profold2.data.parsers import parse_fasta


logger = logging.getLogger(__file__)

def domain_index_split(text):
    seq_index_pattern = '(\d+)-(\d+)'
    for s in text.split(','):
      r = re.match(seq_index_pattern, s)
      assert r
      yield tuple(map(int, r.group(1, 2)))

def process(item, args=None):  # pylint: disable=redefined-outer-name
  input_file, domain_index = item
  logger.info('process: %s', input_file)

  positions = list(domain_index_split(domain_index))

  filename, _ = os.path.splitext(os.path.basename(input_file))
  with open(input_file, 'r') as f:
    fasta_string = f.read()
  sequences, descriptions = parse_fasta(fasta_string)
  assert len(sequences) == 1
  for seq, pid in zip(sequences, descriptions):
    pid = pid.split()[0]
    if args.filename_as_pid:
      pid = filename
    for k, (i, j) in enumerate(positions):
      with open(os.path.join(args.prefix, f'{pid}_{k}.fasta'), 'w') as f:
        f.write(f'>{pid}@{k} {pid}:{i}-{j}\n')
        f.write(seq[i-1:j])
      logger.info('pid: %s@%s->%s', pid, k, input_file)
    break
  return input_file

def main(args):  # pylint: disable=redefined-outer-name
  logger.info('args - %s', args)

  os.makedirs(args.prefix, exist_ok=True)
  inputs = [(args.inputs[i], args.inputs[i+1]) for i in range(0, len(args.inputs), 2)]

  with mp.Pool() as p:
    for _ in p.imap(functools.partial(process, args=args), inputs):
      pass


if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument('inputs', metavar='file', type=str, nargs='+',
      help='inputs')
  parser.add_argument('-o', '--prefix')
  parser.add_argument('--filename_as_pid', action='store_true', help='')
  args = parser.parse_args()

  main(args)
