#! /bin/bash

# テーブル形式のファイルをテキスト処理で転置する(行と列を入れ替える)スクリプト.
# こちらが本体.
# メモリに乗り切らない大きめのファイルを扱う時に使う.
# 目安としてメモリ32GBのメモリで255系統・6,000,000SNPsを
# PythonやRのデータフレーム形式で読み込んで
# 転置しようとすると、そこそこ工夫しないと途中でメモリ不足になる.

# I/Oが多くかなりのボトルネックになるので、時間はかかる。

# 引数
INPUT=$1;
OUTPUT=$2;
CHUNK=${3:-500000}
#echo $CHUNK

TMPDIR="xxxTMPDIRxxx"
rm -rf $TMPDIR;
mkdir $TMPDIR;

# 対象のファイルをチャンクごとに転置して出力する。
python 15_transpose_txt.py -i $INPUT -c $CHUNK

# 各チャンクの改行コードを変換する
for chunk in `ls -v $TMPDIR`;
  do tr -d "\r" < $TMPDIR/$chunk > $TMPDIR/s_$chunk;
  rm -rf $TMPDIR/$chunk;
done;

file_list=(`ls -v $TMPDIR`);
cp $TMPDIR/${file_list[0]} $OUTPUT;

# 結合
for FILE in ${file_list[@]:1};
  do paste $OUTPUT $TMPDIR/$FILE > tmp_out.txt;
  mv tmp_out.txt $OUTPUT;
done;

rm -rf $TMPDIR