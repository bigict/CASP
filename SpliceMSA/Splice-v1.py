import os
import argparse

def read_msa(file_path):
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
    return sequences

def find_max_depth(domain_msas):
    """
    找到给定结构域MSA的最大深度。
    """
    max_depth = max(len(msa) for msa in domain_msas.values())
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

def del_pad_sequences(sequences):
    """
    对序列列表进行预处理
    """
    # 删除原序列行，即sequence[0]
    sequences.pop(0)

    # 右侧填充'-'
    max_length = max(len(seq) for seq in sequences)
    padded_sequences = [seq.ljust(max_length, '-') for seq in sequences]
    return padded_sequences

def splice_msas(domain_paths):
    """
    将给定结构域的MSA进行横向拼接。
    """
    # 读取每个结构域的MSA
    domain_msas = {}
    for domain, path in domain_paths.items():
        msa = read_msa(path)
        domain_msas[domain] = msa

    # 对每个结构域的序列进行预处理，删除原序列行并填充为等长序列
    for domain, msa in domain_msas.items():
        domain_msas[domain] = del_pad_sequences(msa)

    # 找到所有 domain_msas 中最大的深度（序列数量）
    max_depth = find_max_depth(domain_msas)

    # 对每个结构域的序列进行预处理，使深度一致
    for domain, msa in domain_msas.items():
        domain_msas[domain] = pad_msas(msa, max_depth)


    # 开始横向拼接
    spliced_msa = ""
    for i in range(max_depth):
        for domain in domain_paths.keys():
            spliced_msa += domain_msas[domain][i]
        spliced_msa += '\n'  # 每行以'\n'结尾

    return spliced_msa

def write_spliced_msa(output_path, spliced_msa):
    """
    将拼接后的MSA写入文件。
    """
    with open(output_path, 'w') as file:
        file.write(spliced_msa)

def generate_domain_paths(input_folder,protein_name):
    #查找前缀是protein_name的目录eg:T1125Domain1
    domain_dirs = [dir for dir in os.listdir(input_folder) if dir.startswith(protein_name) and os.path.isdir(os.path.join(input_folder,dir))]
    
    domain_paths = {}
    for domain in domain_dirs:
        domain_path = os.path.join(input_folder, domain, 'msas', 'bfd_uniref_hits.a3m')
        domain_paths[domain] = domain_path

    return domain_paths

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Splice MSA")
    parser.add_argument('-i','--input_folder',type=str,help='${CWD}/data/${timestamp}/msa')
    parser.add_argument('-n', '--protein_name', type=str, help='Name of the protein')
    parser.add_argument('-o', '--output_folder', type=str, help='Output folder for spliced MSA file')
    args = parser.parse_args()

    if not args.input_folder or not args.protein_name or not args.output_folder:
        print("Error: All input_folder, protein_name, and output_folder arguments are required.")
        exit(1)

    domain_paths= generate_domain_paths(args.input_folder,args.protein_name)
    
    spliced_msa = splice_msas(domain_paths)

    output_path = os.path.join(args.output_folder,f'{args.protein_name}-SpliceDone.a3m')

    write_spliced_msa(output_path, spliced_msa)

    print(f"Spliced MSA saved to {output_path}")