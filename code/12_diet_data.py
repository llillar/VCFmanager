#! /usr/local/bin/python3
#! coding: utf-8
'''
Python >= 3.7

12_diet_data.py
    -i (--input-file-path)
    -o (--output-file-path)


データ量が多くメモリに乗り切らない計算を行う場合において
データを削減するスクリプト。デフォルトでは1/10に削減する。
PCAなどデータを要約する場合向け。
'''


import argparse
import datetime
from logging import getLogger, StreamHandler, FileHandler, INFO, Formatter
import os
import random
import subprocess
import sys
import time


sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/src")
from my_utils import Runtime_counter


def main():
    random.seed(0)
    ################ Setting command line arguments ################
    parser=argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # 入力ファイルのパス(必須)
    parser.add_argument(
        "-i", "--input-file-path", type=str, action="store",
        dest="inputFilePath", required=True, help="Path to input file.")
    
    # 出力ファイルのパス(必須)
    parser.add_argument(
        "-o", "--output-file-path", type=str, action="store",
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
    fh = FileHandler(
        filename=__file__ + datetime.datetime.now().isoformat() +".log")
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

    # シェルコマンドで行数を数える
    proc_res: subprocess.CompletedProcess = subprocess.run(
        args=["wc", "-l", input_file_path], shell=False, 
        stdout=subprocess.PIPE)
    # 入力ファイルが存在しない場合
    if not proc_res.stdout.decode():
        logger.info("Error!")
        logger.info("Maybe input file does not exist.")
        logger.info("Suspend the process.")
        logger.info("=======================================================")
        sys.exit()
    num_lines: int = int(proc_res.stdout.decode().split(" ")[0])

    # 全体のdiet_rate分の1をランダムに選び出力する。
    # ただし1行目はヘッダーとして必ず残す。
    # in演算子を使う場合listよりsetの方が高速なのでsetに変換する。
    diet_rate: int = 10
    outlines: set = \
        set([1] + random.sample(range(2, num_lines+1), k=int(num_lines/diet_rate)))
    i: int = 1
    try:
        with open(input_file_path, "r") as input_file, \
            open(output_file_path, "w") as output_file:
            for line in input_file:
                if i in outlines:
                    output_file.write(line)
                i += 1
    except FileNotFoundError as fene:
        logger.info("Error!")
        logger.info(f"File: {fene.filename} does not exisit.")
        logger.info("Suspend the process.")
        logger.info("=======================================================")
        sys.exit()

    end: float = time.time()
    logger.info("Success processing!")
    logger.info(f"Run Time = {Runtime_counter(start, end)} seconds")
    logger.info("=======================================================")
    ################ End of main process ################


if __name__=="__main__":
    main()