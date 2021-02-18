#! /usr/local/bin/python3
#! coding: utf-8
'''
Python >= 3.7
VCF version4.2
(https://samtools.github.io/hts-specs/VCFv4.2.pdf)

00_before_imputation.py
    -i (--input-file-path)
    -o (--output-file-path)

BeagleによるImputationを行う際の前処理用スクリプト。
VCFのData lineを対象とし、そのうち genotype fieldからGT(genotype)だけを取り出す。
これが2倍体のフォーマットに沿わない場合、欠損値(./.)に変換して出力する。
GTAKで2倍体にも関わらず半数体のジェノタイプが出たことがあり、
下流の解析に詰まったことがあるため。
'''

import argparse
import datetime
from logging import getLogger, StreamHandler, FileHandler, INFO, Formatter
import sys
import time
from typing import List

sys.path.append("./src")
from my_utils import Runtime_counter
from my_vcf import Remain_only_GT


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
    
    args = parser.parse_args()
    input_file_path: str = args.inputFilePath
    output_file_path: str = args.outputFilePath
    ################ End of setting command line arguments ################


    ################ Setting of logger ################
    logger = getLogger(__name__)
    logger.setLevel(INFO)
    sh = StreamHandler()
    sh.setLevel(INFO)
    sh.setFormatter(Formatter("%(asctime)s %(message)s"))
    fh = FileHandler(filename=__file__ + datetime.datetime.now().isoformat() +".log")
    fh.setLevel(INFO)
    fh.setFormatter(Formatter("%(asctime)s %(message)s"))
    logger.addHandler(sh)
    logger.addHandler(fh)
    ################ End of setting of logger ################


    ################ Main process ################
    start: float = time.time()

    logger.info(__file__ + f"\n\
        \t\t\t\t--input_file_path {input_file_path}\n\
        \t\t\t\t--output_file_path {output_file_path}\n")
    logger.info("=======================================================")
    logger.info("Start program...")

    try:
        with open(input_file_path, "r") as input_file, \
            open(output_file_path, "w") as output_file:
            for line in input_file:
                line: str = line.rstrip("\n|\r|\r\n")
                if line.startswith("#"): # Meta-information or header line
                    output_file.write(line + "\n")
                else: #Data line
                    splited_line: List[str] = line.split("\t")
                    splited_line[8] = "GT"
                    splited_line[9:] = \
                        list(map(Remain_only_GT, splited_line[9:]))
                    new_line: str = "\t".join(splited_line)
                    output_file.write(new_line + "\n")
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
    logger.info("Next step is Imputation!")
    logger.info("=======================================================")
    ################ End of main process ################

if __name__=="__main__":
    main()