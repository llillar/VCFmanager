#! /usr/bin/python3
#! coding: utf-8
'''
Python >= 3.7
VCF version4.2
(https://samtools.github.io/hts-specs/VCFv4.2.pdf)

10_after_imputation.py
    -i (--input-file-path)
    -o (--output-file-path)
    -cr (--convert-rule)
    -mM (--min-MAF)
    -mN (--max-NA)
    -rf (--remove-fields)

BeagleによるImputationを行った後、RやPythonで解析を進めるための前処理用スクリプト。
ジェノタイプを数値データに変換し、不要な行、列を除く。
'''

import argparse
from collections import Counter
import datetime
from logging import basicConfig, getLogger, StreamHandler, FileHandler, info, INFO, Formatter
import re
import sys
import time
from typing import Counter, List, Pattern

from src.my_utils import Runtime_counter, Multi_pop
from src.my_vcf import Check_alt, GT2numeric, Remain_only_GT, Calc_MAF, Calc_NA_rate, Change_chrom


def main():
    ################ Setting command line arguments ################
    parser=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # 入力ファイルのパス(必須)
    parser.add_argument("-i", "--input-file-path", type=str, action="store",
        dest="inputFilePath", required=True, help="Path to input file.")
    
    # 出力ファイルのパス(必須)
    parser.add_argument("-o", "--output-file-path", type=str, action="store",
        dest="outputFilePath", required=True, help="Path to output file.")
    
    # 各ジェノタイプの変換ルール
    # [野生型ホモ:ヘテロ:変異型ホモ]
    parser.add_argument("-cr", "--convert-rule", type=str, action="store", 
        dest="convert_rule", default="[1:0:-1]",
        help="Conversion rule for converting genotype to numeric data.\
        Specify in the following order. \
        REF:HETERO:ALT (default=[1:0:-1])")
    
    # マイナーアレル頻度によるフィルタリング
    # (デフォルトはNA、フィルタリングしない)
    parser.add_argument("-mM", "--min-MAF", action="store",dest="min_MAF", \
        default="NA", help="SNP below min-MAF will be removed.0 ~ 0.5 \
            (default=\"NA\", nothing will be removed)\
            ")
    
    # 欠損値の割合によるフィルタリング
    # (デフォルトはNA、フィルタリングしない)
    parser.add_argument("-mN", "--max-NA", action="store", dest="max_NA", \
        default="NA", help="SNP above max-NA will be removed. 0 ~ 1 \
            (default=\"NA\", nothing will be removed)\
            ")

    
    # VCFの不要なフィールドを指定
    # (デフォルトは何も指定していない、False)
    parser.add_argument("-rf", "--remove-fields", action="store",
        dest="remove_fields", default=False,
        help="Fileds of VCF to remove. Specify any or all of, \
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT \
            separated by colons(:) (default=False)")
    
    args = parser.parse_args()
    input_file_path: str = args.inputFilePath
    output_file_path: str = args.outputFilePath

    # []つきで受け取る
    # -1などが先頭に来ると他の引数と認識されるため
    if not args.convert_rule.startswith("[") or not args.convert_rule.endswith("]"):
        print(f"Argument --convert-rule must be enclosed in []")
        sys.exit()
    # []を取り除いてリストに変換する
    convert_rule: List[str] = list(args.convert_rule[1:-1].split(":"))

    min_MAF: Union[str, NA] = args.min_MAF
    if min_MAF != "NA":
        min_MAF = float(min_MAF)
        if min_MAF < 0.0 or min_MAF > 0.5:
            print("min_MAF must be 0 ~ 0.5")
            sys.exit()
    
    max_NA: Union[str, NA] = args.max_NA
    if max_NA != "NA":
        max_NA = float(max_NA)
        if max_NA < 0.0 or max_NA > 1.0:
            print("max_NA must be 0 ~ 1")
            sys.exit()

    # []つきで受け取る
    if args.remove_fields:
        if not args.remove_fields.startswith("[") or not args.remove_fields.endswith("]"):
            print("Argument --remove-fields must be enclosed in []")
            sys.exit()
    # []を取り除いてリストに変換する、何も指定されていない場合空のリストを作る。
    remove_fields: List[str] = \
        args.remove_fields[1:-1].split(":") if args.remove_fields else []
    
    All_fields: dict = {"CHROM":0, "POS":1, "ID":2, "REF":3, "ALT":4, \
        "QUAL":5, "FILTER":6, "INFO":7, "FORMAT":8}
    remove_fields_index: List[int] = []
    for field in remove_fields:
        try:
            remove_fields_index.append(All_fields[field])
        except KeyError:
            print(f"{field}は入力ファイルに含まれていません。")
            sys.exit()
    ################ End of setting command line arguments ################


    ################ Setting of logger ################
    logger = getLogger(__name__)
    logger.setLevel(INFO)
    sh = StreamHandler()
    sh.setLevel(INFO)
    sh.setFormatter(Formatter("%(asctime)s %(message)s"))
    fh = FileHandler(filename=__file__ \
        + datetime.datetime.now().isoformat() +".log")
    fh.setLevel(INFO)
    fh.setFormatter(Formatter("%(asctime)s %(message)s"))
    logger.addHandler(sh)
    logger.addHandler(fh)
    ################ End of setting of logger ################


    ################ Main process ################
    start: float = time.time()
    
    logger.info(__file__ + f"\n\
        \t\t\t\t--input_file_path {input_file_path}\n\
        \t\t\t\t--output_file_path {output_file_path}\n\
        \t\t\t\t--convert-rule {convert_rule}\n\
        \t\t\t\t--min-MAF {min_MAF}\n\
        \t\t\t\t--max-NA {max_NA}\n\
        \t\t\t\t--remove-fields {remove_fields}\n")
    logger.info("=======================================================")
    logger.info("Start program...")

    #count_line: int = 0
    count_SNPs: int = 0
    multi_alt_site: int = 0
    under_MAF_site: int = 0
    above_NA_site: int = 0
    try:
        with open(input_file_path, "r") as input_file, \
            open(output_file_path, "w") as output_file:
            for line in input_file:
                line: str = line.rstrip("\n|\r|\r\n")
                if line.startswith("##"): # Meta-information line
                    pass # Meta-information lineは除く
                elif line.startswith("#CHROM"): # Header line
                    splited_line: List[str] = line.split("\t")
                    # #CHROMの#の部分は要らない。
                    # Rで読み込めなくなるから。
                    splited_line[0] = "CHROM"
                    splited_line = Multi_pop(splited_line, remove_fields_index)
                    new_line = "\t".join(splited_line)
                    output_file.write(new_line + "\n")
                else: # Data line
                    splited_line: List[str] = line.split("\t")

                    # 縦棒が残っているとMAFの計算に影響が出るので変換する
                    # Remain_only_GTを使うのは、
                    # 生のVCFから直接このスクリプトを動かす時に必要なため
                    splited_line[9:] = \
                        list(map(Remain_only_GT, splited_line[9:]))
                    if Check_alt(splited_line[4]):
                        multi_alt_site += 1 # multi allelic siteの場合は書き出さない
                    elif min_MAF != "NA" and Calc_MAF(splited_line[9:]) <= min_MAF:
                        under_MAF_site += 1 # min_MAF以下のSNPは書き出さない
                    elif max_NA != "NA" and Calc_NA_rate(splited_line[9:]) >= max_NA:
                        above_NA_site += 1 # max_NA以上のNAの割合のSNPは書き出さない
                    else:
                        # #CHROM fieldを染色体番号だけに変える。
                        splited_line[0] = Change_chrom(splited_line[0])

                        # ID fieldになにも記述がなければ("."ならば)
                        # "染色体番号"-"物理位置"の形式に書き換える。
                        if splited_line[2] == ".":
                            splited_line[2] = \
                                splited_line[0] + "-" +splited_line[1]
                        
                        # GTを数値データに変換する
                        splited_line[9:] = GT2numeric(splited_line[9:], convert_rule)

                        # 不要な列を除く
                        splited_line = \
                            Multi_pop(splited_line, remove_fields_index)
                        
                        new_line: str = "\t".join(splited_line)
                        output_file.write(new_line + "\n")
                        count_SNPs += 1
                #count_line += 1 
    except FileNotFoundError as fene:
        logger.info("Error!")
        logger.info(f"File: {fene.filename} does not exisit.")
        logger.info("Suspend the process.")
        logger.info("=======================================================")
        sys.exit()
    except UnicodeDecodeError:
        logger.info("Error!")
        logger.info("Maybe your file is compressed.")
        logger.info("Check it out.")
        logger.info("Suspend the process.")
        logger.info("=======================================================")
        sys.exit()
    
    end: float = time.time()
    logger.info("Success processing!")
    logger.info(f"Run Time = {Runtime_counter(start, end)} seconds")
    logger.info(f"{count_SNPs} SNPs were written in your {output_file_path} .")
    if multi_alt_site:
        logger.info(f"{multi_alt_site} SNPs were multi allelic site, \
            and they were removed.")
    if under_MAF_site:
        logger.info(f"{under_MAF_site} SNPs were under {min_MAF}, and they were removed.")
    if above_NA_site:
        logger.info(f"{above_NA_site} SNPs were above {max_NA}, and they were removed.")
    if remove_fields:
        logger.info(f"Field: {remove_fields} were removed.")
    logger.info("=======================================================")
    ################ End of main process ################


if __name__=="__main__":
    main()