import subprocess
import os
import sys

def sto_to_a3m_with_perl(sto_file):
    command = f"perl reformat.pl {sto_file} uniref90_hits.a3m"
    # 执行命令行
    subprocess.run(command, shell=True, check=True)

def main():
    if sys.platform.startswith('win'):
        print("当前是 Windows 环境")
    else:
        print("无法确定当前操作系统")
        return
    
    sto_file_path = r"D:\CASP\CASP\SpliceMSA\uniref90_hits.sto"
    
    print(os.path.exists(sto_file_path))
    
    sto_to_a3m_with_perl(sto_file_path)

if __name__ == "__main__":
    main()