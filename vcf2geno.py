#! /usr/bin/python3
#! coding: utf-8

'''
このスクリプトはPython 3.6以上でサポートされています。


このスクリプトは入力VCFから、選択したモードに応じて以下のファイルを書き出します。

SCANモード:入力VCFの10行目以降、各系統のジェノタイプデータのうちのGTのみ残したものを出力します。
例) 1/1:0,25:25:75:1115,75,0 -> 1/1
また、フォーマットに沿わないジェノタイプデータが含まれる場合、それらは取り除かれることなく出力されます。
これらに関してはログにその情報が記載されるので、そちらを参考にして下さい。

これは一部のVariant CallerもしくはVariant Call時のBAMファイルが原因で
VCFの定義に沿わないジェノタイプデータを出力してしまうという不具合に対応させたものです。


PCAモード:入力VCFの10行目以降、各系統のジェノタイプデータを以下の規則のもとGTのみ残し数値データに変換したものを出力します。
0/0(0|0) -> -1
1/1(1|1) -> 1
0/1(0|1), 1/0(1|0) -> 0
例) 1/1:0,25:25:75:1115,75,0 -> 1
なお欠損値(./.)が含まれる場合はプログラムが中断されます。
このためPCAモードを利用する際はBeagle等で欠損値の穴埋め(imputation)を行ったものを入力VCFとして下さい。


GWASモード:入力VCFの10行目以降、各系統のジェノタイプデータを以下の規則のもとGTのみ残し数値データに変換したものを出力します。
0/0(0|0) -> -1
1/1(1|1) -> 1
else -> NA
例) 1/1:0,25:25:75:1115,75,0 -> 1
ジェノタイプが野生型ホモの場合-1、変異型ホモの場合1に変換されます。
ヘテロもしくは欠損値はNAに変換されます。
なお、この形式はRのrrBLUPに対応した形式です。

必須パラメータは-i(--input_file_path),-o(--output_file_path),-m(--mode)です。
-mに関してはSCAN, PCA, GWASの中から出力したい形式に応じて選んで下さい。
##で始まるheader部分を残したい場合は-k(--keep-header)のオプションをTrueにして下さい。デフォルトはFalseです。

入出力ファイルは.vcf(.vcf.gz)もしくは.txtをサポートしています。
入出力に圧縮形式を利用する場合は、処理時間が大幅に伸びる可能性があります。ディスクに余裕がある場合は、非圧縮形式を推奨します。
multi allelic site(変異が2種類以上ある箇所)は自動的に取り除かれます。
'''


import argparse
import datetime
import gzip
from logging import basicConfig, getLogger, StreamHandler, FileHandler, DEBUG, INFO, Formatter
import os
import re
import sys
import time
from typing import IO, List, Pattern, Union


## Setting of logger
logger = getLogger(__name__)
logger.setLevel(DEBUG)

sh = StreamHandler()
sh.setLevel(INFO)
sh.setFormatter(Formatter("%(asctime)s %(levelname)8s %(message)s"))

fh = FileHandler(filename=__file__ + datetime.datetime.now().isoformat() +".log")
fh.setLevel(DEBUG)
fh.setFormatter(Formatter("%(asctime)s %(levelname)8s %(message)s"))

logger.addHandler(sh)
logger.addHandler(fh)


def runtime_counter(start: float, end: float) -> str: 
    """
    Arguments:
    ===
        start: float
            Start time
        end: float
            Finish time
        You can get each argumemts using time.time()

    Returns:
    ===
        str
            Run time which calculated end - start
    """
    runtime: str = str(round((end - start), 2))
    return runtime


def open_file(file_path: str, mode: str) -> Union[IO, bool]:
    """
    Arguments:
    ===
        file_path: str
            Path to input/output VCF
            .vcf, .vcf.gz, .txt format is supported
        mode: str
            
    Returns:
    ===
        IO
            If input/output file is supported, return opened file
        bool
            If input/output file is not supported, return False
    """
    if file_path.endswith(".vcf.gz"):
        return gzip.open(file_path, mode=mode+"t")
    elif file_path.endswith(".vcf") or file_path.endswith(".txt"):
        return open(file_path, mode=mode)
    else:
        return False


def check_alt(alt: str) -> bool:
    """
    Arguments:
    ===
        alt: str
            Each ALT field of input VCF
    Returns:
    ===
        bool
            If input VCF contains multi allelic site, return True
    """
    #alt: str = line.split("\t", 9)[4]
    if alt.count(","):
        return True


def remain_only_GT(geno: str) -> str:
    """
    Arguments:
    ===
        geno: str
            Each genotype field of input VCF
    Returns:
    ===
        str
           Return only "GT" in genotype field
    """
    #geno: str = line.split("\t", 9)[9]
    ##VCFのgenotype fieldは定義上、下記の正規表現で全てカバーできるはず。
    ##Variant callerによっては当てはまらない場合があるかも。
    geno_pattern: Pattern = re.compile(r"([.|\d][/|][.|\d])[0-9,|_\-:\.ATGC]+")
    return re.sub(geno_pattern, r"\1", geno)


def scan_GT(GT: str) -> Union[set, bool]:
    """
    Arguments:
    ===
        GT: str
            GT of genotype field of input VCF
    Returns:
    ===
        set
            If input VCF has unexpected GT, return itself
        bool
            If input VCF dose not have unexpected GT, return False
    """
    expected_GT: set = {"0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1", "./."}
    observed_GT: set = set(GT.split("\t"))
    unexpected_GT: set = observed_GT - expected_GT
    if unexpected_GT:
        return unexpected_GT
    else:
        return False


def GT2numeric(GT: str, mode: str) -> Union[str, bool]:
    """
    Arguments:
    ===
        GT: str
            GT of genotype field of input VCF
        mode:  str
            PCA or GWAS
    Returns:
    ===
        str
            GT will be changed by the following rules to numeric data

            If PCA is used as mode, change GT field as follow:
                0/0(0|0) -> -1
                1/1(1|1) -> 1
                0/1(0|1), 1/0(1|0) -> 0
                if GT contains ./. , return False
            if GWAS is use as mode, change GT field as follow:
                0/0(0|0) -> -1
                1/1(1|1) -> 1
                else -> NA
    """
    if mode == "PCA":
        if GT.count("./."):
            return False
        else:
            GT: str = GT.replace("0/0","-1").replace("0|0","-1")
            GT: str = GT.replace("1/1","1").replace("1|1","1")
            GT: str = GT.replace("0/1","0").replace("0|1","0")
            GT: str = GT.replace("1/0","0").replace("1|0","0")
    if mode == "GWAS":
        GT: str = GT.replace("0/0","-1").replace("0|0","-1")
        GT: str = GT.replace("1/1","1").replace("1|1","1")
        GT: str = GT.replace("0/1","NA").replace("0|1","NA")
        GT: str = GT.replace("1/0","NA").replace("1|0","NA")
        GT: str = GT.replace("./.", "NA")
    return GT


def main():
    parser=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    parser.add_argument("-i", "--input_file_path", type=str, action="store",
        dest="inputFilePath", required=True, help="Path to input file")
    parser.add_argument("-o", "--output_file_path", type=str, action="store",
        dest="outputFilePath", required=True, help="Path to output file")
    parser.add_argument("-m", "--mode", type=str, action="store",
        dest="mode", required=True, choices=["SCAN", "PCA", "GWAS"])
    parser.add_argument("-k", "--keep-header", action="store_true",
        dest="keep_header", default=False,
        help="If use this argument, leave vcf header.")

    args=parser.parse_args()

    input_file_path: str = args.inputFilePath
    output_file_path: str = args.outputFilePath
    mode: str = args.mode
    keep_header: bool = args.keep_header

    start: float = time.time()

    logger.debug(__file__ + f"\n\
        \t\t\t\t--input_file_path {input_file_path}\n\
        \t\t\t\t--output_file_path {output_file_path}\n\
        \t\t\t\t--mode {mode}\n\
        \t\t\t\t--keep-header {keep_header}\n")


    logger.info("=======================================================")
    logger.info("Start program...")

    
    if not open_file(input_file_path, "r"):
        logger.info("Invalid format!")
        logger.info("Only VCF format (.vcf, .vcf.gz, .txt) is supported for input files.")
        sys.exit()
    else:
        input_file: IO = open_file(input_file_path, "r")

    if not open_file(output_file_path, "w"):
        logger.info("Invalid format!")
        logger.info("Only VCF format (.vcf, .vcf.gz, .txt) is supported for output files.")
        sys.exit()
    else:
        output_file: IO = open_file(output_file_path, "w")
    
    count_line: int = 1
    count_SNPs: int = 0
    unexpected_line: int = 0
    multi_alt_site: int = 0
    
    for line in input_file:
        line: str = line.rstrip("\n|\r|\r\n")
        if line[0:2] == "##": #headerのうち##から始まるもの
            if keep_header:
                output_file.write(line + "\n")
            else:
                pass
        elif line[0:6] == "#CHROM": #headerのうち#CHROMから始まるもの
            if keep_header:
                output_file.write(line + "\n")
            else:
                output_file.write(line[1:] + "\n") #"#CHROMの"#を取り除いて出力
        else: #bodyの部分
            line: List[str] = line.split("\t", 9)
            alt: str = line[4]
            geno: str = line[9]
            if check_alt(alt):
                multi_alt_site += 1 #multi allelic siteは書き出さない
            else:
                GT: str = remain_only_GT(geno)
                line[9] = GT
                if mode == "SCAN":
                    scan_res: Union[str, bool] = scan_GT(GT)
                    if scan_res: # unexpected_GTがある場合logに書き出す。
                        unexpected_line += 1
                        logger.debug(f"Line {count_line} has unexpected GT:{scan_res}.")
                    line: str = "\t".join(line)
                    output_file.write(line + "\n")
                    count_SNPs += 1
                else:
                    new_GT: Union[str, bool] = GT2numeric(GT, mode)
                    if not new_GT:
                        logger.info("Invalid format in input VCF!")
                        logger.info("PCA mode cannot allow ./. genotype in VCF.")
                        logger.info("If you want to use your VCF for PCA, you need to imputate your VCF by using Beagle.")
                        sys.exit()
                    else:
                        line[9] = new_GT
                        line: str = "\t".join(line)
                        output_file.write(line + "\n")
                        count_SNPs += 1    
        count_line += 1

    input_file.close()
    output_file.close()
    end: float = time.time()
    logger.info("Success processing!")
    logger.info(f"Run Time = {runtime_counter(start, end)} seconds")
    logger.info(f"{count_SNPs} SNPs were written in your {output_file_path} .")
    if unexpected_line:
        logger.info(f"However your {input_file_path} has {unexpected_line} invalid format lines in genotype field and keep remain.")
        logger.info(f"Check log file")
    if multi_alt_site:
        logger.info(f"{multi_alt_site} SNPs were multi allelic site, and they were removed.")

    logger.info("=======================================================")




if __name__=="__main__":
    main()