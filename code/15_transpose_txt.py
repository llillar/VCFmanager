#! /usr/local/bin/python3
#! coding: utf-8
'''
Python >= 3.7
pandas==1.2.2

テーブル形式のファイルをテキスト処理で転置する(行と列を入れ替える)スクリプト.
本体は同名のシェルスクリプト.

PCA(主成分分析)などの前処理用.
メモリに乗り切らない大きめのファイルを扱う時に使う.
目安としてメモリ32GBのメモリで255系統・6,000,000SNPsを
PythonやRのデータフレーム形式で読み込んで
転置しようとすると、そこそこ工夫しないと途中でメモリ不足になる.
'''

import argparse

import pandas as pd

def main():
    ################ Setting command line arguments ################
    parser=argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # 入力ファイルのパス(必須)
    parser.add_argument(
        "-i", "--input-file-path", type=str, action="store",
        dest="inputFilePath", required=True, help="Path to input file.")
    
    # チャンクサイズ(一度に読み込む行)
    parser.add_argument(
        "-c", "--chunk-size", type=int, action="store",
        dest="chunk_size", default=500000, required=True, 
        help="Chunk size(lines) to read at one time. (default=500000)")
    
    args = parser.parse_args()
    input_file_path: str = args.inputFilePath
    chunk_size: int = args.chunk_size
    ################ End of setting command line arguments ################


    ################ Main process ################
    df: pd.DataFrame = pd.read_table(input_file_path, chunksize=chunk_size)
    i: int =  1
    for chunk in df:
        filename: str = f"xxxTMPDIRxxx/chunk{i}.txt"
        if i == 1:
            chunk.T.to_csv(filename, header=None, sep="\t")
        else:
            chunk.T.to_csv(filename, index=None, header=None, sep="\t")
        i += 1
    ################ End of main process ################

if __name__=="__main__":
    main()