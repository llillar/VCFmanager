#! /usr/local/bin/python3
#! coding: utf-8
'''
このモジュールはよく使う汎用的な??関数をまとめたものです。
'''

from typing import Any, List, Tuple

def Runtime_counter(start: float, end: float) -> str: 
    """
    This function returns runtime from start to end.
    
    Arguments:
    ----------
    start: float
        Start time
    end: float
        End time
    Can be got from time module time function.
    
    Returns:
    ----------
    runtime: str
        Runtime which calculated end - start.
    """
    runtime: str = str(round((end - start), 2))
    return runtime


def Get_index(target_list: List[Any], element: Any) -> Tuple[int, None]:
    """
    Get index of element from target_list
    
    Arguments:
    ----------
    target_list: List[Any]
        description
    element: Any

    Returns:
    ----------
    return1: type
        description
    """
    try:
        index: int = target_list.index(element)
        return index
    except ValueError:
        # リストに無いelementが指定された時。
        # boolだとintと区別が難しいので、Noneを返す。
        return None


def Multi_pop(target_list: List[Any], index_list: List[int]) -> list:
    """
    This function removes element(s) 
    from a target_list using index_list.
    
    Arguments:
    ----------
    target_list: list
        Target list which you want to remove some element(s).
    index_list: List[int]
        Index number(s) list of the element(s) 
        to be removed from the target_list.
    
    Returns:
    ----------
    return: list
        List with some element(s) removed.
    """
    index_list.sort(reverse=True)
    # IndexErrorを避ける処理
    while index_list[0] >= len(target_list):
        index_list.pop(0)
    # target_listから除く
    for i in index_list:
        target_list.pop(i)
    return target_list

def main():
    print("Hello, this is my_utils.py")

if __name__=="__main__":
    main()