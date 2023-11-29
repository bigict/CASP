# brutal process using reformat.pl from HH-suite3
import subprocess
import os
import argparse
import re

def read_msa_from_bfd(file_path):
    """
    读取给定路径的MSA文件，并返回序列列表。
    """
    sequences = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith(">"):
                # 忽略注释行
                continue
            sequences.append(line.strip())
    
    #对序列列表进行预处理
    # 删除原氨基酸序列即sequence[0]
    sequences.pop(0)
    return sequences

def read_msa_from_sto(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("#") or line.startswith("//"):
                # 忽略注释行
                continue
            match = re.match(r"\S+\s+([A-Za-z.-]+)", line)
            if match:
                sequences.append(match.group(1))#group(1)返回第一个括号匹配的内容

    #对序列列表进行预处理
    # 删除原氨基酸序列即sequence[0]
    sequences.pop(0)
    return sequences

def find_max_depth(bfd_domain_msas):
    """
    找到给定结构域MSA的最大深度。
    """
    max_depth = max(len(msa) for msa in bfd_domain_msas.values())
    return max_depth

def pad_msas(msa, max_depth):
    """
    对msa添加或者不添加行，使得其列表长度等于max_depth。
    """
    depth_diff = max_depth - len(msa)
    if depth_diff > 0:
        # 添加行
        msa += ['-' * len(msa[0])] * depth_diff
    return msa

def pad_sequences(sequences):

    # 右侧填充'-'
    max_length = max(len(seq) for seq in sequences)
    padded_sequences = [seq.ljust(max_length, '-') for seq in sequences]
    return padded_sequences

def sto_splice_msas(uniref90_paths, mgnify_paths):
    """
    将给定结构域的UNIREF90和MGNIFY的MSA进行纵向拼接。
    """
    # 读取每个结构域的UNIREF90和MGNIFY的MSA
    uniref90_domain_msas = {}
    for domain, path in uniref90_paths.items():
        msa = read_msa_from_sto(path)
        uniref90_domain_msas[domain] = msa

    mgnify_domain_msas = {}
    for domain, path in mgnify_paths.items():
        msa = read_msa_from_sto(path)
        mgnify_domain_msas[domain] = msa
    
    #直接将mgnify_domain_msas按照domain的顺序拼接到uniref90_domain_msas后面
    sto_spliced_msas_list = []
    for domain in uniref90_paths.keys():
        sto_domain_msas_len=len(uniref90_domain_msas[domain])+len(mgnify_domain_msas[domain])
        for i in range(sto_domain_msas_len):
            if i<len(uniref90_domain_msas[domain]):
                sto_spliced_msas_list.append(uniref90_domain_msas[domain][i])
            else:
                sto_spliced_msas_list.append(mgnify_domain_msas[domain][i-len(uniref90_domain_msas[domain])])
    
    sto_spliced_msas_list = pad_sequences(sto_spliced_msas_list)

    sto_spliced_msa = ""
    for i in range(len(sto_spliced_msas_list)):
        sto_spliced_msa += sto_spliced_msas_list[i]
        sto_spliced_msa += '\n'

    return sto_spliced_msa
  

def bfd_splice_msas(bfd_paths):
    """
    将给定结构域的MSA进行横向拼接。
    """
    # 读取每个结构域的MSA
    bfd_domain_msas = {}
    for domain, path in bfd_paths.items():
        msa = read_msa_from_bfd(path)
        bfd_domain_msas[domain] = msa

    # 对每个结构域的序列进行预处理，充为等长序列
    for domain, msa in bfd_domain_msas.items():
        bfd_domain_msas[domain] = pad_sequences(msa)

    # 找到所有 bfd_domain_msas 中最大的深度（序列数量）
    max_depth = find_max_depth(bfd_domain_msas)

    # 对每个结构域的序列进行预处理，使深度一致
    for domain, msa in bfd_domain_msas.items():
        bfd_domain_msas[domain] = pad_msas(msa, max_depth)


    # 开始横向拼接
    bfd_spliced_msa = ""
    for i in range(max_depth):
        for domain in bfd_paths.keys():
            bfd_spliced_msa += bfd_domain_msas[domain][i]
        bfd_spliced_msa += '\n'  # 每行以'\n'结尾

    return bfd_spliced_msa



def write_spliced_msa(output_path, spliced_msa):
    """
    将拼接后的MSA写入文件。
    """
    with open(output_path, 'w') as file:
        file.write(spliced_msa)

def generate_bfd_paths(input_folder,protein_name):
    #查找前缀是protein_name+“Domain"的目录eg:T1125Domain1
    domain_dirs = [dir for dir in os.listdir(input_folder) if dir.startswith(protein_name+'Domain') ]
    bfd_paths = {}
    for domain in domain_dirs:
        domain_path = os.path.join(input_folder, domain, 'msas', 'bfd_uniref_hits.a3m')
        bfd_paths[domain] = domain_path

    return bfd_paths

def generate_uniref90_paths(input_folder,protein_name):
    #查找前缀是protein_name+“Domain"的目录eg:T1125Domain1
    domain_dirs = [dir for dir in os.listdir(input_folder) if dir.startswith(protein_name+'Domain') ]
    uniref90_paths = {}
    for domain in domain_dirs:
        uniref90_path = os.path.join(input_folder, domain, 'msas', 'uniref90_hits.sto')
        uniref90_paths[domain] = uniref90_path

    return uniref90_paths

def generate_mgnify_paths(input_folder,protein_name):
    #查找前缀是protein_name+“Domain"的目录eg:T1125Domain1
    domain_dirs = [dir for dir in os.listdir(input_folder) if dir.startswith(protein_name+'Domain') ]
    mgnify_paths = {}
    for domain in domain_dirs:
        mgnify_path = os.path.join(input_folder, domain, 'msas', 'mgnify_hits.sto')
        mgnify_paths[domain] = mgnify_path

    return mgnify_paths


def main():
    parser=argparse.ArgumentParser(description="Splice MSA")
    parser.add_argument('-i','--input_folder',type=str,help='${CWD}/data/${timestamp}/msa')
    parser.add_argument('-n', '--protein_name', type=str, help='Name of the protein')
    parser.add_argument('-o', '--output_folder', type=str, help='Output folder for spliced MSA file')
    args = parser.parse_args()

    if not args.input_folder or not args.protein_name or not args.output_folder:
        print("Error: All input_folder, protein_name, and output_folder arguments are required.")
        exit(1)

    bfd_paths= generate_bfd_paths(args.input_folder,args.protein_name)
    uniref90_paths= generate_uniref90_paths(args.input_folder,args.protein_name)
    mgnify_paths= generate_mgnify_paths(args.input_folder,args.protein_name)

    bfd_spliced_msa = bfd_splice_msas(bfd_paths)
    sto_spliced_msa = sto_splice_msas(uniref90_paths, mgnify_paths)
    spliced_msa = bfd_spliced_msa + sto_spliced_msa


    output_path = os.path.join(args.output_folder,f'{args.protein_name}-SpliceDone.a3m')

    write_spliced_msa(output_path, spliced_msa)

    print(f"Spliced MSA with bdf&uniref90&mgnify saved to {output_path}")

if __name__ == "__main__":
    main()
    