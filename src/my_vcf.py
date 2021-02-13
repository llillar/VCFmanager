#! /usr/bin/python3
#! coding: utf-8
'''
このモジュールはVCFハンドリング関連の関数をまとめたものです。
'''

from collections import Counter
from typing import Counter, List

def Check_alt(alt: str) -> bool:
    """
    This function check ALT field of VCF. 
    If ALT has more than two alleles, this function will return True. 
    
    Arguments:
    ----------
    alt: str
        ALT field of input VCF.
    
    Returns:
    ----------
    return: bool
        If ALT has more than two alleles, return True.
        If not, return False.
    """
    #alt: str = line.split("\t")[4]
    #各Data lineの5列目
    if alt.count(","):
        return True
    else:
        return False


def Remain_only_GT(geno: str) -> str:
    """
    This function will return GT(genotype) from Genotype field.
    
    Arguments:
    ----------
    geno: str
        Genotype field of input VCF.
    
    Returns:
    ----------
    GT: str
        GT(genotype) from Genotype field.
    """
    expected_GT: List[str] = ["0/0", "0/1", "1/0", "1/1", "./."]

    #geno_list: List[str] = line.split("\t")[9:]
    #各Data lineの9列目以降、各列がgenotype fieldに相当
    GT: str = geno.split(":")[0].replace("|", "/")

    #GATKを使った際にフォーマットに則らないジェノタイプが出てきたことがある。(仕様？)
    #そうしたジェノタイプは欠損値に変換する。
    if GT not in expected_GT:
        GT = "./."
    return GT


def GT2numeric(GT_list: List[str], convert_rule: List[str]) -> List[str]:
    """
    This function change GT to numeric data.
    
    Arguments:
    ----------
    GT_list: List[str]
        GT list generated by remain_only_GT function.
    convert_rule: List[str]
        [REF, HETERO, ALT]
    
    Returns:
    ----------
    num_list: List[str]
        GT will be changed by the following rules to numeric data or NA.
            0/0 -> REF
            0/1, 1/0 -> HETERO
            1/1 -> ALT
            ./. -> NA
    """
    REF, HETERO, ALT = convert_rule
    # 一度文字列にして置換するのが一番速い??
    tmp = "\t"/join(GT)
    tmp.replace("0/0", REF).replace("0/1", HETERO).replace("1/0", HETERO).replace("1/1", ALT).replace("./.", "NA")
    num_list: List[str] = tmp.split("\t")
    return num_list


def calc_MAF(GT_list: List[str]) -> float:
    """
    Calculate Minor Allele Frequency.
    
    Arguments:
    ----------
    GT_list: List[str]
        GT list generated by remain_only_GT function.
    
    Returns:
    ----------
    MAF: float
        Minor Allele Frequency.
    """
    Count_Summary: Counter = Counter(GT_list)
    Genotyped: int = (len(GT_list) - Count_Summary("./.")) * 2
    ALT_num: int = Count_Summary("0/1") + Count_Summary("1/0") + Count_Summary("1/1") * 2
    AAF: float = ALT_num / Genotyped
    MAF: float = min(AAF, 1 - AAF)
    return MAF


def calc_NA_rate(GT_list: List[str]) -> float:
    """
    Description of this function.
    
    Arguments:
    ----------
    GT_list: List[str]
        GT list generated by remain_only_GT function.
    
    Returns:
    ----------
    NA_rate: float
        Percentage of NA.
    """
    Count_Summary: Counter = Counter(GT_list)
    Num: int = len(GT_list)
    NA_Num: int = sum(Count_Summary["./."])
    NA_rate: float = NA_Num / Num
    return NA_rate


if __name__=="__main__":
    main()