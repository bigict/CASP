"""
msa_assemble Module

without using multimer pair and merge, directly horizontal splice

Author: TanQian
"""
#import subprocess
import os
import argparse
import re
#import collections
#import itertools
#import string
import sys

from alphafold.data import parsers
#from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Set


def write_spliced_msa(output_path, spliced_msa):
  with open(output_path, 'w', encoding='utf-8') as f:
    f.write(spliced_msa)

  print(f'saved to {output_path}')


def generate_bfd_paths(input_folder, protein_name):
  #查找前缀是protein_name+“Domain"的目录eg:T1125Domain1
  domain_dirs = sorted([
      dir for dir in os.listdir(input_folder)
      if dir.startswith(protein_name + 'Domain')
  ],
                       key=lambda x: int(x[len(protein_name + 'Domain'):]))

  if not domain_dirs:
    print(f'No matching directories found for {protein_name} in bfd norms\n')
    sys.exit(1)

  print(f'Found {len(domain_dirs)} matching directories: {domain_dirs}')
  bfd_paths = {}
  for domain in domain_dirs:
    domain_path = os.path.join(input_folder, domain, 'msas',
                               'bfd_uniref_hits.a3m')
    bfd_paths[domain] = domain_path

  return bfd_paths


def generate_uniref90_paths(input_folder, protein_name):
  #查找前缀是protein_name+“Domain"的目录eg:T1125Domain1
  domain_dirs = sorted([
      dir for dir in os.listdir(input_folder)
      if dir.startswith(protein_name + 'Domain')
  ],
                       key=lambda x: int(x[len(protein_name + 'Domain'):]))
  if not domain_dirs:
    print(f'No matching directories found for {protein_name} in uniref90\n')
    sys.exit(1)

  uniref90_paths = {}
  for domain in domain_dirs:
    uniref90_path = os.path.join(input_folder, domain, 'msas',
                                 'uniref90_hits.sto')
    uniref90_paths[domain] = uniref90_path
  return uniref90_paths


def generate_mgnify_paths(input_folder, protein_name):
  #查找前缀是protein_name+“Domain"的目录eg:T1125Domain1
  domain_dirs = sorted([
      dir for dir in os.listdir(input_folder)
      if dir.startswith(protein_name + 'Domain')
  ],
                       key=lambda x: int(x[len(protein_name + 'Domain'):]))
  if not domain_dirs:
    print(f'No matching directories found for {protein_name} in mgnify\n')
    sys.exit(1)

  mgnify_paths = {}
  for domain in domain_dirs:
    mgnify_path = os.path.join(input_folder, domain, 'msas', 'mgnify_hits.sto')
    mgnify_paths[domain] = mgnify_path
  return mgnify_paths


def read_msa_from_bfd(
    bfd_paths):  #bfd_paths是一个字典，key是domain，value是bfd_uniref_hits.a3m的路径
  bfd_msas = {}
  for domain, path in bfd_paths.items():
    with open(path, 'r', encoding='utf-8') as file:
      bfd_msas[domain] = file.read()
  return bfd_msas


def read_msa_from_uniref90(uniref90_paths):  #同时转换成a3m格式
  uniref90_msas = {}
  for domain, path in uniref90_paths.items():
    with open(path, 'r', encoding='utf-8') as file:
      temp = file.read()
      uniref90_msas[domain] = parsers.convert_stockholm_to_a3m(temp)
  return uniref90_msas


def read_msa_from_mgnify(mgnify_paths):
  mgnify_msas = {}
  for domain, path in mgnify_paths.items():
    with open(path, 'r', encoding='utf-8') as file:
      temp = file.read()
      mgnify_msas[domain] = parsers.convert_stockholm_to_a3m(temp)
  return mgnify_msas


def combine_db_msa(bfd_msas, uniref90_msas, mgnify_msas):
  #将bfd_msas,uniref90_msas,mgnify_msas合并成一个字典，key是domain，value是一个字符串，是拼接后的MSA
  db_msas = {}
  for domain in bfd_msas.keys():
    print(f'Combining {domain} MSAs')
    db_msas[
        domain] = bfd_msas[domain] + uniref90_msas[domain] + mgnify_msas[domain]
  return db_msas


# def combine_domain(db_msas):
#   # 将 db_msas 根据 domain 进行横向拼接
#   domain_names = sorted(list(db_msas.keys()),
#                         key=lambda x: int(re.search(r'\d+', x).group()))

#   # domain_msas_lines 是个列表，其中的每个元素是一个 domain 的 MSA，该 MSA 被按照换行符拆分成了行的列表。
#   domain_msas_lines = [msa.split('\n') for msa in db_msas.values()]

#   # 统一domain内部长度
#   for i in range(len(domain_msas_lines)):
#     #处理各个domain：
#     max_length = max(len(line) for line in domain_msas_lines[i])
#     for j in range(len(domain_msas_lines[i])):
#       if len(domain_msas_lines[i][j]) < max_length:
#         while len(domain_msas_lines[i][j]) < max_length:
#           domain_msas_lines[i][j] += '-'

#   # 统一各个domaini的深度
#   max_depth = max(len(msa) for msa in domain_msas_lines)
#   for i in range(len(domain_msas_lines)):
#     if len(domain_msas_lines[i]) < max_depth:
#       while len(domain_msas_lines[i]) < max_depth:
#         domain_msas_lines[i].append('-' * len(domain_msas_lines[i][0]))

#   # 横向拼接
#   spliced_msa = ''
#   for row in range(len(domain_msas_lines[0])):
#     for i in range(len(domain_names)):
#       spliced_msa += domain_msas_lines[i][row]
#     spliced_msa += '\n'

#   return spliced_msa


def combine_domain(db_msas, protein_name, input_folder):
  spliced_msa = ''
  with open(os.path.join('seq_parse/FUpredOutput',
                         f'{protein_name}-FUpred.txt'),
            'r',
            encoding='utf-8') as fu_file:
    domain_info = fu_file.read()

  domain_info = domain_info.split(';')
  domain_info_list = [domain.split(',') for domain in domain_info]
  #domain_info_list示例：['1-78','111-138','228-321','79-110','139-196'],现在将这些列表元素按照首个数字排序
  domain_info_list = sorted(domain_info_list,
                            key=lambda x: int(x[0].split('-')[0]))
  print(f'domain_info_list: {domain_info_list}')

  #将head_line和msa_line视为一组，所以max_depth/2
  max_depth = max(len(msa) in db_msas.values()) / 2

  for row_index in range(max_depth):
    head_line = ''
    msa_line = ''
    #开始对这一行按照domain_info_list的顺序进行横向拼接
    for i, region in enumerate(domain_info_list):
      #region:eg:'1-78',是个字符串
      start = int(region.split('-')[0])
      end = int(region.split('-')[1])
      msa_len = end - start + 1

      #如果domain_info_list不能覆盖原有序列，则以'-'填充
      if i == 0:
        if start > 1:
          head_line += '>' + ' ' * (start - 1)
          msa_line += '-' * (start - 1)
      else:
        if start > 1 + int(domain_info_list[i - 1].split('-')[1]):
          head_line += ' ' * (start - 1 -
                              int(domain_info_list[i - 1].split('-')[1]))
          msa_line += '-' * (start - 1 -
                             int(domain_info_list[i - 1].split('-')[1]))

      found_index = None
      for j, domain in enumerate(domain_info, start=1):
        if region in domain:
          print(f'{region}  is in Domain{j}: {domain}')
          found_index = j
          break

      if found_index is None:
        print(f'{region} not found in domain_info_list')
        sys.exit(1)

      #找到了含有region的domain后，找到这个domain的msa
      #input_folder= msa_tool_data_dir/msa
      region_msa = db_msas[
          f'{input_folder}/{protein_name}Domain{found_index}'].split('>')
      head_line += region_msa[row_index].split('\n')[0]
      msa_line += region_msa[row_index].split('\n')[1]
    spliced_msa += head_line + '\n' + msa_line + '\n'

  return spliced_msa


def main():
  parser = argparse.ArgumentParser(description='Assemble MSA')
  parser.add_argument('-i',
                      '--input_folder',
                      type=str,
                      help='${CWD}/data/${timestamp}/msa')
  parser.add_argument('-n',
                      '--protein_name',
                      type=str,
                      help='Name of the protein')
  parser.add_argument('-o',
                      '--output_folder',
                      type=str,
                      help='Output folder for assembled MSA file')
  args = parser.parse_args()

  if not args.input_folder or not args.protein_name or not args.output_folder:
    print(
        'Error:input_folder,protein_name,output_folder arguments are required.')
    sys.exit(1)
  try:
    bfd_paths = generate_bfd_paths(args.input_folder, args.protein_name)
  except SystemExit as e1:
    print(f'find no a3m file in {args.input_folder},exit code: {e1.code}\n')
    sys.exit(1)
  try:
    uniref90_paths = generate_uniref90_paths(args.input_folder,
                                             args.protein_name)
  except SystemExit as e2:
    print(f'find no sto file in {args.input_folder},exit code: {e2.code}\n')
    sys.exit(1)
  try:
    mgnify_paths = generate_mgnify_paths(args.input_folder, args.protein_name)
  except SystemExit as e3:
    print(f'find no sto file in {args.input_folder},exit code: {e3.code}\n')
    sys.exit(1)

  #bfd_msas是一个字典，key是domain，value是bfd_uniref_hits.a3m的内容,是一个字符串
  #此时uniref90&mgnify都转成了a3m格式
  bfd_msas = read_msa_from_bfd(bfd_paths)
  uniref90_msas = read_msa_from_uniref90(uniref90_paths)
  mgnify_msas = read_msa_from_mgnify(mgnify_paths)
  db_msas = combine_db_msa(bfd_msas, uniref90_msas, mgnify_msas)

  spliced_msa = combine_domain(db_msas, args.protein_name, args.input_folder)
  output_path = os.path.join(args.output_folder,
                             f'{args.protein_name}-assembled.a3m')
  write_spliced_msa(output_path, spliced_msa)


if __name__ == '__main__':
  try:
    main()
  except SystemExit as e:
    print('msa_assemble.py has no domains to assemble,thus Canceled\n')
    sys.exit(1)
  else:
    print('msa_assemble.py run successfully\n')
    sys.exit(0)
