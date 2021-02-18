#! /usr/local/bin/python3
#! coding: utf-8
'''
Python >= 3.7
numpy==1.20.1
pandas==1.2.2
scikit-learn==0.24.1

10_after_imputation.py
    -i (--input-file-path)
    -od (--output-dir)

入力ファイルの想定
ID    sample1    sample2    sample3    sample4
SNP1     1          1          0          1
SNP2    -1          0          1         -1
SNP3    -1          1          1         -1
SNP4     1          1          0          0
SNP5     0         -1         -1          0
'''


import argparse
import datetime
from logging import getLogger, StreamHandler, FileHandler, INFO, Formatter
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/src")
from my_utils import Runtime_counter


def main():
    ################ Setting command line arguments ################
    parser=argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # 入力ファイルのパス(必須)
    parser.add_argument(
        "-i", "--input-file-path", type=str, action="store",
        dest="inputFilePath", required=True, help="Path to input file.")
    
    # ファイルの出力先(必須)
    parser.add_argument(
        "-od", "--output-dir", type=str, action="store",
        dest="output_dir", required=True, help="Directory to output files.")
    
    # 正規化するか否か
    # parser.add_argument(
    #     "-s", "--standardized", type=bool, action="store", \
    #     dest="standardized", default=True, help="If True, standardize data.\
    #     (default=True)")
    
    args = parser.parse_args()
    input_file_path: str = args.inputFilePath
    out_dir: str = args.output_dir
    # Make directory if does not exist.
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    # standardized: bool = args.standardized
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
        \t\t\t\t--output_dir {out_dir}\n")
    logger.info("=======================================================")
    logger.info("Start program...")

    # Reading data as pandas dataframe
    try:
        logger.info("Reading data...")
        df: pd.DataFrame = pd.read_table(input_file_path, index_col=0)
    except FileNotFoundError as fene:
        logger.info("Error!")
        logger.info(f"File: {fene.filename} does not exisit.")
        logger.info("Suspend the process.")
        logger.info("=======================================================")
        sys.exit()

    # Standardizing by each line(SNP).
    logger.info("Standardizing data...")
    df = df.apply(lambda x: (x-x.mean())/x.std(ddof=True), axis=1)
    df = df.T

    # Performing PCA
    logger.info("Performing Principal Component Analysis ...")
    pca: PCA = PCA()
    try:
        pca.fit(df)
    except ValueError as ve:
        logger.info("Error!")
        logger.info("Maybe your input file contains NA.")
        logger.info("Please imputate your file before PCA.")
        logger.info("=======================================================")
        sys.exit()

    res: np.ndarray = pca.transform(df)

    # 主成分スコア
    pca_score: pd.DataFrame = pd.DataFrame(
        data=res,
        columns=[f"PC{x}" for x in range(1, len(df.index)+1)],
        index=df.index)
    pca_score.to_csv(f"{out_dir}/PCA_Score.txt", sep="\t")

    # 寄与率
    evr:pd.DataFrame = pd.DataFrame(
        data=pca.explained_variance_ratio_,
        columns=["explained_variance_ratio"],
        index=[f"PC{x}" for x in range(1, len(df.index)+1)])
    evr.to_csv(f"{out_dir}/Expl_Var_Ratio.txt", sep="\t")

    # # 固有値
    # ev:pd.DataFrame = pd.DataFrame(
    #     data=pca.explained_variance_,
    #     column=["explained_variance"],
    #     index=[f"PC{x}" for x in range(1, len(df.index)+1)])
    # ev.to_csv(f"{out_dir}/Expl_Var.txt", sep="\t")
    end: float = time.time()
    logger.info("Success processing!")
    logger.info(f"Run Time = {Runtime_counter(start, end)} seconds")
    logger.info("=======================================================")
    ################ Main process ################


if __name__=="__main__":
    main()