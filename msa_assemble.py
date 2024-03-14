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


def parse_domain(domain):
  n = 0
  #n：前面包含多少个align的字符，就是Start_at
  for region in map(lambda x: x.split('-'), domain.split(',')):
    s, e = map(int, region)
    assert s <= e
    yield s, e, n
    n += e - s + 1


def split_domains(domain_info):
  i = 0
  #i就是found_index
  for domain in map(lambda x: x.strip(), domain_info.split(';')):
    if domain:
      for region in parse_domain(domain):
        yield i, region
      i += 1


def combine_domain(db_msas, protein_name):
  spliced_msa = ''
  with open(os.path.join('seq_parse/FUpredOutput',
                         f'{protein_name}-FUpred.txt'),
            'r',
            encoding='utf-8') as fu_file:
    domain_info = fu_file.read()

  #domain_info删掉文章末尾可能出现的换行符和';'，并分割
  #domain_info = domain_info.strip().strip(';')
  #domains = domain_info.split(';')
  # domain_info = [
  #     item for sub_string in domain_info.split(';')
  #     for item in sub_string.split(',')
  # ]
  #domain_info示例：['1-78','111-138','228-321','79-110','139-196'],现在将这些列表元素按照首个数字排序
  # domain_info = sorted(domain_info, key=lambda x: int(x.split('-')[0]))
  # print(f'domain_info: {domain_info}')
  # domain_info_len = len(domain_info)
  # print(f'domain_info_len: {domain_info_len}')

  domains = list(split_domains(domain_info))
  print(f'domains: {domains}')
  #domains: [(0, (1, 78, 0)), (0, (111, 138, 78)), (0, (228, 321, 106)), (0, (354, 559, 200)), (0, (603, 687, 406)), (0, (801, 1089, 491)), (0, (1128, 1280, 780)), (0, (1311, 1424, 933)), (1, (79, 110, 0)), (2, (139, 196, 0)), (3, (197, 227, 0)), (4, (322, 353, 0)), (5, (560, 602, 0)), (6, (688, 743, 0)), (7, (744, 800, 0)), (8, (1090, 1127, 0)), (9, (1281, 1310, 0))]
  #将domains按照(0, (1, 78, 0))中(1, 78, 0)的第一个元素排序
  domains = sorted(domains, key=lambda x: x[1][0])

  print(f'domains: {domains}')
  #将head_line和msa_line视为一组，所以max_depth/2.值为db_msas.values()中最大的msa的行数
  max_depth = max(len(msa.split('\n')) for msa in db_msas.values()) // 2
  print(f'max_depth: {max_depth}')

  for row_index in range(max_depth):  #注意row_index从0开始！ 只有domains的名字从1开始
    print(f'assemble the {row_index}th line')
    head_line = ''
    msa_line = ''
    #开始对这一行按照domain_info的顺序进行横向拼接
    for i, region in enumerate(domains):
      #i:region的序号，region:（start,end,start_at）
      print(f'it is the {i} th region:{region}')
      start = region[1][0]
      end = region[1][1]
      start_at = region[1][2]
      msa_len = end - start + 1
      found_index = int(region[0])
      print(f'msa_len:{msa_len}')
      #如果region不能覆盖原有序列，则以'-'填充
      if i == 0:
        if start > 1:
          head_line += '>' + ' ' * (start - 1)
          msa_line += '-' * (start - 1)
      else:
        if start > 1 + int(domains[i - 1][1][1]):
          head_line += ' ' * (start - 1 - int(domains[i - 1][1][1]))
          msa_line += '-' * (start - 1 - int(domains[i - 1][1][1]))

      #found_index是domains中region所在的domain的索引
      # found_index = None
      # region_index = None

      # found_domain_info = None
      # for j, domain in enumerate(domains, start=1):
      #   if region in domain:
      #     found_index = j
      #     found_domain_info = domain
      #     print(f'{region}  is in Domain{j}: {found_domain_info}')
      #     has_region_cnt = len(domain.split(','))
      #     print(f'has_region_cnt:{has_region_cnt}')
      #     #查找该region是domain中的第几个
      #     found_domain_info = found_domain_info.split(',')
      #     for n, sub_region in enumerate(found_domain_info, start=0):
      #       if sub_region == region:
      #         print(f'{region} is the {n}th region in Domain{j}')
      #         region_index = n
      #     break

      # if found_index is None:
      #   print(f'{region} not found in domains')
      #   sys.exit(1)

      region_msa = db_msas[f'{protein_name}Domain{found_index+1}'].split('\n')
      #if 此行对应的domain的msa行数不够，则以'-'填充
      if row_index > ((len(region_msa) // 2) - 1):
        print(
            f'row_index: {row_index} is out of range,max_row_index: {len(region_msa) // 2 - 1}'
        )
        head_line += ' ' * msa_len
        msa_line += '-' * msa_len
      else:
        #msa中只有大写字母或者‘-’计数要达到msa_len个，小写字母保留但不计数.msa_line根据p_start和p_end截取
        #has_region_cnt =domains列表中元素[0]是region[0]的个数
        has_region_cnt = len([x for x in domains if x[0] == region[0]])
        print(f'has_region_cnt:{has_region_cnt}')
        if has_region_cnt == 1:
          print(
              f'the Domain only has 1 region:{region},so we assemble msa_line:{region_msa[row_index * 2 + 1]}'
          )
          head_line += region_msa[row_index * 2]
          msa_line += region_msa[row_index * 2 + 1]
          continue

        #wrong:version:不适用于有小写字母的情况
        # if start == 1:
        #   print(
        #       f'start=1,is the beginning,so we choose Domain{found_index+1}\'s {region} from [{start-1}: {end}],msa_line:{region_msa[row_index * 2 + 1][start - 1:end]}'
        #   )
        #   head_line += region_msa[row_index * 2][start - 1:end]
        #   msa_line += region_msa[row_index * 2 + 1][start - 1:end]
        #   continue

        print(f'start assemble msa at {start_at}')
        p, q = 0, 0
        k = 0
        while k < start_at and start != 1:
          if not region_msa[row_index * 2 + 1][p].islower():
            k += 1
          p += 1
        #untile here: p是该region之前实际有效的氨基酸应该截取到的下标，k=start-1个大写字母
        '''
        eg： region:'1-3','4-5',所以看msa_line[0,1,2,3,4,5,6]：AaBbCDE,结果应该是AaBbC'DE'：msa_line[5:7](左闭右开)
        start_at=3,msa_len=2,p=0,k=0
        |index|char| after_p|after_k|p|k|
        |  0  |  A |   1    |   1   |0|0|
        |  1  |  a |   2    |   1   |1|1|
        |  2  |  B |   3    |   2   |2|1|
        |  3  |  b |   4    |   2   |3|2|
        |  4  |  C |   5    |   3   |4|2|--loop ended
        |  5  |  D |   -    |   -   |0|0|
        |  6  |  E |   -    |   -   |0|0|
        所以现在的region:'4-5'就应该从p=5开始截取（因为index从0开始）
        p=5,q=5
        '''

        q = p
        k = 0
        while k < msa_len:
          if not region_msa[row_index * 2 + 1][q].islower():
            k += 1
          q += 1
        #untile here: q 是msa_line实际要截取到的下标(开区间），k=msa_len
        '''
        eg： region:'1-3','4-5',所以看msa_line[0,1,2,3,4,5,6]：AaBbCDE,结果应该是AaBbC'DE'：msa_line[5:7](左闭右开)
        start-1=3,q=5,k=0,msa_len=2
        |index|char| after_q|after_k|q|k|
        |  5  |  D |   6    |   1   |5|0|
        |  6  |  E |   7    |   2   |6|1|--loop end
        |  7  |  - |   -    |   -   |
        q=7,而我们要截取的是'DE',所以是msa_line[5:7]=[p:q]
        '''

        print(
            f'now we choose Domain{found_index+1}\'s {region} from [{p}:{q}],msa_line:{region_msa[row_index * 2 + 1][p:q]}'
        )
        head_line += region_msa[row_index * 2][p:q]
        msa_line += region_msa[row_index * 2 + 1][p:q]

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

  spliced_msa = combine_domain(db_msas, args.protein_name)
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
