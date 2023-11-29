#use parsers.convert_stockholm_to_a3m, without using multimer pair and merge
import subprocess
import os
import argparse
import re
import collections
import itertools
import re
import string
import sys
sys.path.append('/public/dbu/CASP/alphafold/')

from alphafold.data import parsers
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Set

def generate_bfd_paths(input_folder,protein_name):
    #查找前缀是protein_name+“Domain"的目录eg:T1125Domain1
    domain_dirs = [dir for dir in os.listdir(input_folder) if dir.startswith(protein_name+'Domain') ]
    print(f"Found {len(domain_dirs)} matching directories: {domain_dirs}")
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

def read_msa_from_bfd(bfd_paths):#bfd_paths是一个字典，key是domain，value是bfd_uniref_hits.a3m的路径
    bfd_msas={}
    for domain,path in bfd_paths.items():
        with open(path,'r') as file:
            bfd_msas[domain]=file.read()
    return bfd_msas 

def read_msa_from_uniref90(uniref90_paths):#同时转换成a3m格式
    uniref90_msas={}
    for domain,path in uniref90_paths.items():
        with open(path,'r') as file:
            temp=file.read()
            uniref90_msas[domain]=parsers.convert_stockholm_to_a3m(temp)
    return uniref90_msas

def read_msa_from_mgnify(mgnify_paths):
    mgnify_msas={}
    for domain,path in mgnify_paths.items():
        with open(path,'r') as file:
            temp=file.read()
            mgnify_msas[domain]=parsers.convert_stockholm_to_a3m(temp)
    return mgnify_msas

def combine_db_msa(bfd_msas,uniref90_msas,mgnify_msas):
    #将bfd_msas,uniref90_msas,mgnify_msas合并成一个字典，key是domain，value是一个字符串，是拼接后的MSA
    db_msas={}
    for domain in bfd_msas.keys():
        print(domain)
        db_msas[domain]=bfd_msas[domain]+uniref90_msas[domain]+mgnify_msas[domain]
    return db_msas



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
    
    #bfd_msas是一个字典，key是domain，value是bfd_uniref_hits.a3m的内容,是一个字符串
    #此时uniref90&mgnify都转成了a3m格式
    bfd_msas=read_msa_from_bfd(bfd_paths)
    uniref90_msas=read_msa_from_uniref90(uniref90_paths)
    mgnify_msas=read_msa_from_mgnify(mgnify_paths)
    db_msas=combine_db_msa(bfd_msas,uniref90_msas,mgnify_msas)
    
    spliced_msa = combine_domain(db_msas)
    #output_path = os.path.join(args.output_folder,f'{args.protein_name}-SpliceDone.a3m')
    #write_spliced_msa(output_path, spliced_msa)
    #print(f"Spliced MSA saved to {output_path}")

if __name__ == "__main__":
    main()
    