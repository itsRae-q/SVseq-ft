#!usr/bin/env python
# coding=utf-8

import os


# 存放一个文件的数据信息
class indi_data:

    def __init__(self,data_info):
        self.read = data_info[0]
        self.pattern = data_info[1]
        self.ref = data_info[2]


# 提取特征  match_ratio,mismatch_ratio,del_ratio,ins_ratio,sub_ratio,max_continue_mis_num (最大连续未匹配碱基数),ave_continue_mis_num,max_contine_ins_num,max_contine_del_num,max_contine_sub_num
def extract_feature(indi):
    read,pattern,ref = list(indi.read),list(indi.pattern),list(indi.ref)
    indi_len = len(pattern)

    # initialize
    inita = []
    for i in range(7):
        inita.append(0)
    ins_num, del_num, sub_num, contine_mis_num,contine_ins_num, contine_del_num, contine_sub_num = [x for x in inita]

    contine_ins_set, contine_del_set, contine_sub_set, contine_mis_set = [1],[1],[1],[1]

    for num in range(indi_len):
        if(pattern[num] == "*"):
            if (read[num] != "-" and  ref[num] == "-"):  # 插入
                ins_num += 1
                contine_mis_num += 1
                contine_ins_num += 1

                contine_del_set.append(contine_del_num)
                contine_del_num = 0
                contine_sub_set.append(contine_sub_num)
                contine_sub_num = 0
            if (read[num] == "-" and  ref[num] != "-"):  # 删除
                del_num += 1
                contine_mis_num += 1
                contine_del_num += 1

                contine_ins_set.append(contine_ins_num)
                contine_ins_num = 0
                contine_sub_set.append(contine_sub_num)
                contine_sub_num = 0
            if (read[num] != "-" and  ref[num] != "-"):  # 替换
                sub_num += 1
                contine_mis_num += 1
                contine_sub_num += 1

                contine_ins_set.append(contine_ins_num)
                contine_ins_num = 0
                contine_del_set.append(contine_del_num)
                contine_del_num = 0

        else:
            contine_mis_set.append(contine_mis_num)
            contine_mis_num = 0
            contine_ins_set.append(contine_ins_num)
            contine_ins_num = 0
            contine_del_set.append(contine_del_num)
            contine_del_num = 0
            contine_sub_set.append(contine_sub_num)
            contine_sub_num = 0
    # 去掉0
    contine_mis_set = [x for x in contine_mis_set if x!= 0]
    mis_num = ins_num + del_num + sub_num
    mis_ratio = round(mis_num/indi_len,5)
    match_ratio = round(1-mis_ratio,5)
    ins_ratio = round(ins_num/indi_len,5)
    del_ratio = round(del_num/indi_len,5)
    sub_ratio = round(sub_num/indi_len,5)
    max_continue_mis, max_continue_ins, max_continue_del, max_continue_sub = [x for x in list(map(max,[contine_mis_set,contine_ins_set,contine_del_set,contine_sub_set]))]

    max_continue_mis_ratio, max_continue_ins_ratio, max_continue_del_ratio, max_continue_sub_ratio = max_continue_mis/indi_len, max_continue_ins/indi_len, max_continue_del/indi_len, max_continue_sub/indi_len

    ave_continue_mis = sum(contine_mis_set)/len(contine_mis_set)
    ave_continue_mis_ratio = ave_continue_mis/indi_len

    feature_list = [mis_ratio, match_ratio, ins_ratio, del_ratio, sub_ratio, max_continue_mis_ratio, max_continue_ins_ratio, max_continue_del_ratio, max_continue_sub_ratio, ave_continue_mis_ratio]

    return feature_list


# 对indel_set里的read进行提取特征
def get_indel_feature(total_indel_set):
    feature_list_set = []  # 二维列表，每一个元素是一个read的特征，一个read有10个特征
    for elem in total_indel_set:
        indi = indi_data([elem.read_string,elem.pattern,elem.ref_string])
        feature_list_set.append(extract_feature(indi))
    return feature_list_set
