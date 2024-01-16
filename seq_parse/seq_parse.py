"""
seq_parse Module

This module provides functions for parsing sequences using fupred.

Author: TanQian
"""
import argparse
import os
import subprocess
import sys
from Bio import SeqIO


def parse_arguments():
    parser = argparse.ArgumentParser(description="Sequence Parser")
    # 添加输入文件参数
    parser.add_argument("-i",
                        "--input",
                        help="Input sequence file (e.g., seq.fasta)",
                        required=True)
    # 添加输出文件夹参数
    parser.add_argument("-o",
                        "--output",
                        help="Output folder (e.g., resfolder)",
                        required=True)
    return parser.parse_args()


def create_output_folder(output_folder):
    #创建输出文件夹,如果没有则新建
    folder_path = os.path.abspath(output_folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    return folder_path


def run_fupred(input_file, output_folder):
    # 运行FUpred程序 ./run-FUpred.pl seq.fasta -o fu.txt
    try:

        command = f"./FUpred/run-FUpred.pl {input_file}"
        subprocess.run(command, shell=True, check=True)
        # 获取输入文件的基本名称（去除路径和扩展名）
        input_base_name = os.path.splitext(os.path.basename(input_file))[0]

        with open(os.path.join("FUpredOutput",
                               f"{input_base_name}-FUpred.txt"),
                  "r",
                  encoding="utf-8") as fu_file:
            domain_info = fu_file.read()  #文件中内容全部存在domain_info中
            #print domain_info in terminal:
            print(domain_info)
        #按照结构域信息裁剪fasta文件
        domain_list = domain_info.split(";")  #将domain_info按照;分割，得到一个列表

        if len(domain_list) == 1:
            print("FUpred predicts only one domain in this protein! Stop:)\n")
            #if there is only one domain. exit the function
            return

        for i, domain in enumerate(domain_list, start=1):
            if not domain.strip():
                continue

            output_file = os.path.join(output_folder,
                                       f"{input_base_name}Domain{i}.fasta")
            with open(output_file, "w", encoding="utf-8") as output_fasta:
                #分割好的domain也有1-456,969-1147;的格式,所以将domain1的两个region的序列拼接在一起
                output_fasta.write(
                    f">{input_base_name}_Domain{i}:{domain_list[i-1].strip()}\n"
                )
                for region in domain.split(","):
                    start, end = map(
                        int,
                        region.strip().split("-"))  # 使用strip()去除换行符，再进行分割
                    output_fasta.write(
                        str(
                            SeqIO.read(input_file,
                                       "fasta").seq[start - 1:end]))
                output_fasta.write("\n")

    except subprocess.CalledProcessError as error_message:
        print("an error occurred while running FUpred: ", error_message)
        sys.exit(1)


def main():
    # 解析命令行参数
    args = parse_arguments()
    # 创建输出文件夹
    output_folder = create_output_folder(args.output)
    # 打印一些信息，以便确认程序是否正确执行
    print(f"Input file: {args.input}")  #f代表format，可以将变量的值直接放入字符串中
    print(f"Output folder: {output_folder}")
    # 运行FUpred程序
    run_fupred(args.input, output_folder)


if __name__ == "__main__":
    try:
        main()
    except argparse.ArgumentError:
        print("Please input the correct parameters!")
        sys.exit(1)
