#!usr/bin/env python
# coding=utf-8

from sklearn.externals import joblib
import numpy as np

# 复杂indel标签是1，其他为0

class com_indel:
    def __init__(self,indel,readinfo):
        self.read_id = indel.read_id
        # self.read_str = indel.read_string
        # self.pattern = indel.pattern
        # self.ref_str = indel.ref_string
        count_s = 0
        for word in readinfo.tstring[:indel.start_pos]:
            if(word=="-"):
                count_s += 1
        # if(indel.start_pos==8969):
        #     a = readinfo.tstring[:indel.start_pos]
        #     b = indel.ref_string
        #     print("error")
        self.ref_start_pos = indel.start_pos - count_s + indel.ref_start  #复杂indel在ref上真正的开始位置（去掉-之后的位置）

        count_e = 0
        for word in readinfo.tstring[:indel.end_pos]:
            if(word=="-"):
                count_e += 1
        self.ref_end_pos = indel.end_pos - count_e + indel.ref_start #复杂indel在ref的结束位置
        self.length = self.ref_end_pos - self.ref_start_pos
        self.new_read = ""
    def add_read(self,read):
        self.new_read = read


# com_indel_set中每一个元素是一个复杂indel,即com_indel类
def get_com_indel(total_indel_set,readInfo):

    clf = joblib.load("svm.model")
    com_indel_set = []
    for indel in total_indel_set:
        com_indel_set.append(com_indel(indel, readInfo[indel.read_id - 1]))
    return com_indel_set



