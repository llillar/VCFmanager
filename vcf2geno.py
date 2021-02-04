#! /usr/bin/python3
#! coding: utf-8

'''
このスクリプトはPython 3.6以上でサポートされています。

このスクリプトはVCFを入力ファイルとし、下流の解析フローへ適した形へと変換します。
特に何も指定しない場合は次のようなファイルになります。
・ヘッダーが取り除かれます。
・多型が複数含まれるSNP(Multi allelic site)が取り除かれます。
・各系統のジェノタイプデータが下記の様にして数値データに変換されます。
    野生型ホモ -> 1
    ヘテロ -> 0
    変異型ホモ -> -1

以下はオプションの説明です。
・VCFのヘッダーを残したい場合は -kh(--keep-header)のフラグを立て下さい。
・ジェノタイプデータを数値データに変換したくない場合は -kg(--keep-GT)のフラグを立て下さい。
  この場合ジェノタイプデータは下記のようにGT部分のみを残した形になります。
  1/1:0,25:25:75:1115,75,0 -> 1/1
・ジェノタイプを数値データに変換する際の規則を変更したい場合は -cr(--convert-rule)で指定して下さい。
  例のようにコロンで区切り 野生型ホモ:ヘテロ:変異型ホモ の順番で指定して下さい。
  例) -cr -1:NA:1
・マイナーアレル頻度によるフィルタリングを行いたい場合は -mM(--min-MAF)で指定して下さい。
  0 ~ 0.5の間で指定し、指定された頻度以下のアレルを持つSNP(行)は取り除かれます。
・欠損値の割合によるフィルタリングを行いたい場合は -mN(--max-NA)で指定して下さい。
  指定した割合以上の欠損値を含むSNP(行)は取り除かれます。
・入力VCF中で不要なフィールド(列)があれば -rf(--remove-fields)で指定して下さい。
  例のようにコロンで区切り指定して下さい。
  仕様上、系統も指定して取り除くことができます。
  例) -rf POS:ID:FILTER:SAMPLE2

入出力ファイルは.vcf(.vcf.gz)もしくは.txtをサポートしています。
入出力に圧縮形式を利用する場合は、処理時間が大幅に伸びる可能性があります。ディスクに余裕がある場合は、非圧縮形式を推奨します。
'''


import argparse
from collections import Counter
import datetime
import gzip
from logging import basicConfig, getLogger, StreamHandler, FileHandler, info, INFO, Formatter
import os
import re
import sys
import time
from typing import IO, List, Pattern, Union


################ Setting of logger ################
logger = getLogger(__name__)
logger.setLevel(INFO)

sh = StreamHandler()
sh.setLevel(INFO)
sh.setFormatter(Formatter("%(asctime)s %(levelname)8s %(message)s"))

fh = FileHandler(filename=__file__ + datetime.datetime.now().isoformat() +".log")
fh.setLevel(INFO)
fh.setFormatter(Formatter("%(asctime)s %(levelname)8s %(message)s"))

logger.addHandler(sh)
logger.addHandler(fh)
################ End of setting of logger ################

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


# def GT2numeric(GT: str, mode: str) -> Union[str, bool]:
#     """
#     Arguments:
#     ===
#         GT: str
#             GT of genotype field of input VCF
#         mode:  str
#             PCA or GWAS
#     Returns:
#     ===
#         str
#             GT will be changed by the following rules to numeric data

#             If PCA is used as mode, change GT field as follow:
#                 0/0(0|0) -> -1
#                 1/1(1|1) -> 1
#                 0/1(0|1), 1/0(1|0) -> 0
#                 if GT contains ./. , return False
#             if GWAS is use as mode, change GT field as follow:
#                 0/0(0|0) -> -1
#                 1/1(1|1) -> 1
#                 else -> NA
#     """
#     if mode == "PCA":
#         if GT.count("./."):
#             return False
#         else:
#             GT: str = GT.replace("0/0","-1").replace("0|0","-1")
#             GT: str = GT.replace("1/1","1").replace("1|1","1")
#             GT: str = GT.replace("0/1","0").replace("0|1","0")
#             GT: str = GT.replace("1/0","0").replace("1|0","0")
#     if mode == "GWAS":
#         GT: str = GT.replace("0/0","-1").replace("0|0","-1")
#         GT: str = GT.replace("1/1","1").replace("1|1","1")
#         GT: str = GT.replace("0/1","NA").replace("0|1","NA")
#         GT: str = GT.replace("1/0","NA").replace("1|0","NA")
#         GT: str = GT.replace("./.", "NA")
#     return GT


def main():
    ################ Setting command line arguments ################
    parser=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    ##入力ファイルのパス(必須)    
    parser.add_argument("-i", "--input-file-path", type=str, action="store",
        dest="inputFilePath", required=True, help="Path to input file.")
    
    ##出力ファイルのパス(必須)
    parser.add_argument("-o", "--output-file-path", type=str, action="store",
        dest="outputFilePath", required=True, help="Path to output file.")
    
    ##ヘッダーを残すか否か(デフォルトは残さない)
    parser.add_argument("-kh", "--keep-header", action="store_true",
        dest="keep_header", default=False,
        help="If use this argument, leave vcf header.(default=False)")

    ##ジェノタイプフィールドを数値データに変換するか否か(デフォルト(False)は変換する)
    parser.add_argument("-kg", "--keep-GT", action="store_true",
        dest="keep_GT", default=False,
        help="If use this argument, genotype will not be converted to numeric.\
            (default=False)")
    
    ##各ジェノタイプの変換ルール 野生型ホモ:ヘテロ:変異型ホモ(keep_GTがTrueのときは使用されない)
    parser.add_argument("-cr", "--convert-rule", type=str, action="store", 
        dest="convert_rule", default="1:0:-1",
        help="Conversion rule for converting genotype to numeric data.\
        Specify in the following order. \
        REF:HETERO:ALT (default=1:0:-1)")
    
    ##マイナーアレル頻度によるフィルタリング(デフォルトはFalse、フィルタリングしない)
    parser.add_argument("-mM", "--min-MAF", type=Union[float, bool], action="store",
        dest="min_MAF", default=0.0, help="SNP below min-MAF will be removed.\
        0 ~ 0.5 (default=False, nothing will be removed)")
    
    ##欠損値の割合によるフィルタリング(デフォルトはFalse、フィルタリングしない)
    parser.add_argument("-mN", "--max-NA", type=Union[float, bool], action="store",
        dest="max_NA", default=1.0, help="SNP above max-NA will be removed.\
        0 ~ 1 (default=False, nothing will be removed)")
    
    ##VCFの不要なフィールドを指定(デフォルトは何も指定していない、False)
    ##仕様上、系統を指定して除くことも出来る。
    parser.add_argument("-rf", "--remove-fields", type=Union[str,bool], action="store",
        dest="remove_fields", default=False,
        help="Fileds of VCF to remove. Specify as the following, separated by colons(:)\
        ID:QUAL:FORMAT:SAMPLE2 (default=False)")


    args = parser.parse_args()

    input_file_path: str = args.inputFilePath
    output_file_path: str = args.outputFilePath
    keep_header: bool = args.keep_header
    keep_GT: bool = args.keep_GT
    convert_rule: list = list(arg.convert_rule.split(":")) ##リストに変換する
    min_MAF: Union[float, bool] = args.min_MAF
    max_NA: Union[float, bool] = args.min_NA
    if min_MAF and (min_MAF < 0.0 or min_MAF > 0.5):
        logger.info("min_MAF must be 0 ~ 1")
        sys.exit()
    if max_NA and (max_NA < 0.0 or max_NA > 1.0):
        logger.info("max_NA must be 0 ~ 1")
        sys.exit()
    remove_fields: Union[list, bool] = \
        args.remove_field.split(":") if args.remove_field else False ##指定された場合リストに変換する
    
    ################ End of setting command line arguments ################

    start: float = time.time()

    logger.info(__file__ + f"\n\
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
            # else:
            #     GT: str = remain_only_GT(geno)
            #     line[9] = GT
            #     if mode == "SCAN":
            #         scan_res: Union[str, bool] = scan_GT(GT)
            #         if scan_res: # unexpected_GTがある場合logに書き出す。
            #             unexpected_line += 1
            #             logger.info(f"Line {count_line} has unexpected GT:{scan_res}.")
            #         line: str = "\t".join(line)
            #         output_file.write(line + "\n")
            #         count_SNPs += 1
            #     else:
            #         new_GT: Union[str, bool] = GT2numeric(GT, mode)
            #         if not new_GT:
            #             logger.info("Invalid format in input VCF!")
            #             logger.info("PCA mode cannot allow ./. genotype in VCF.")
            #             logger.info("If you want to use your VCF for PCA, you need to imputate your VCF by using Beagle.")
            #             sys.exit()
            #         else:
            #             line[9] = new_GT
            #             line: str = "\t".join(line)
            #             output_file.write(line + "\n")
            #             count_SNPs += 1    
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