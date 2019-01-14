#!usr/bin/env python
# coding=utf-8

import numpy as np
from collections import Counter

class csv:

    def __init__(self,left_bp,right_bp,correct_str):
        self.left_bp = left_bp
        self.right_bp = right_bp
        self.correct_str = correct_str

# 计算距离矩阵和操作矩阵
def get_operate_matrix(long_string, short_string):

    row_num, col_num = len(short_string) + 1, len(long_string) + 1

    #  initialize matrix
    dis_matrix = np.zeros((row_num, col_num))
    for i in range(1,row_num):
        dis_matrix[i][0] = i
    for j in range(1,col_num):
        dis_matrix[0][j] = j

    operate_matrix = np.zeros((row_num -1, col_num-1))

    for i in range(1,row_num):
        for j in range(1,col_num):
            if(long_string[j-1] == short_string[i-1]):
                cost = 0
            else:
                cost = 1

            # edit_matrix
            dis_matrix[i][j] = min(dis_matrix[i][j-1]+1, dis_matrix[i-1][j]+1, dis_matrix[i-1][j-1]+cost)

            # operate_matrix
            if (dis_matrix[i][j] == dis_matrix[i][j-1]+1):
                operate_matrix[i-1][j-1] = 1
            else:
                if(dis_matrix[i][j] == dis_matrix[i-1][j-1] + cost):
                    operate_matrix[i - 1][j - 1] = 0
                else:
                    operate_matrix[i-1][j-1] = 2




    # print(dis_matrix)
    # print(operate_matrix)
    return operate_matrix

# 比对两个字符串  input:ATCAG   ATTCAA  output: AT-CAG ATTCAA
def compare_str(long_string, short_string):

    operate_matrix = get_operate_matrix(long_string,short_string)
    row_num, col_num = operate_matrix.shape
    long_string_list, short_string_list = list(long_string), list(short_string)

    i,j = row_num, col_num
    # if(i!=0 and j!=0):
    # 0表示对应位置替换，下次看左上方
    # 1表示左方shor_string下个位置插入-，下次看左边
    # 2表示上方long_string下个位置插入-，下次看上方
    while (True):
        if (operate_matrix[i - 1][j - 1] == 1):
            short_string_list.insert(i, "-")
            j -= 1
        if (operate_matrix[i - 1][j - 1] == 2):
            long_string_list.insert(j, "-")
            i -= 1
        if (operate_matrix[i - 1][j - 1] == 0):
            i -= 1
            j -= 1
        if (i == 0 and j == 0):
            break
        if (j == 0 and i != 0):
            for num in range(i):
                long_string_list.insert(0, "-")
            break
        if (i == 0 and j != 0):
            for num in range(j):
                short_string_list.insert(0, "-")
            break

    # long_string, short_string = "".join(long_string_list), "".join(short_string_list)

    return long_string_list, short_string_list

# 修正对照的string(1.此处的输入为string对应的list,不是string )   input A-TTC-  ATT-C   output:A-TT-C- and 插入-的位置
def updata_source_str(source, sample):
    i,j = 0, 0
    add_pos = [] # 记录sample中插入-的位置
    source_copy,sample_copy = source.copy(), sample.copy()

    while(True):
        if (source_copy[i] == sample_copy[j]):
            i += 1
            j += 1
        else:
            if(source_copy[i] == "-" and sample_copy[j] != "-"):
                sample_copy.insert(j,"-")
                #add_pos.append(j-len(add_pos))  # 减去len(add_pos)是为了记录插入的-在sample中的原始位置
                add_pos.append(j)
                j += 1
                i += 1
            if(source_copy[i] != "-" and sample_copy[j] == "-"):
                source_copy.insert(i,"-")
                i += 1
                j += 1
        if(i>len(source_copy)-1 ):
            for num in range(i,len(sample_copy)):
                source_copy.insert(num,"-")
            break
        if(j>len(sample_copy)-1):
            for num in range(j,len(source_copy)):
                sample_copy.insert(num,"-")
                add_pos.append(j)
                j += 1
            break

    return source_copy,sample_copy,add_pos

#求出一个比对的标准字符串  input:{A-TTC-  ATT-C -A-TTC} output:-A-TT-C-
def get_target_str(source_set):


    target_str = source_set[0]
    # 计算比对标准字符串
    for num in range(1,len(source_set)):
        target_str = updata_source_str(target_str,source_set[num])[0]
    # 求添加-的位置


    add_pos_list = []
    for num in range(len(source_set)):
        add_pos_list.append(updata_source_str(target_str,source_set[num])[2])

    return target_str,add_pos_list

# 将所有的indel对齐  string_set为比过一次后的str，不是原始的str
def get_compare_result(string_set,add_pos_list):


    string_set_copy = string_set.copy()
    for (string,pos) in zip(string_set_copy,add_pos_list):
        if(pos != []):
            for num in pos:
                string.insert(num,"-")

    return string_set_copy

# input: 同一类的indel集合，输出矫正后的indel_string
def self_correct(indel_set):

    first_indel = indel_set[0]
    source_set = []
    string_set = []
    source_set.append(list(first_indel))
    string_set.append(list(first_indel))

    for indel in indel_set[1:]:
        comp_str = compare_str(first_indel,indel)
        source_set.append(comp_str[0])    # first_indel比对后结果
        string_set.append(comp_str[1])    # 其他indel比对后结果

    pos = get_target_str(source_set)[1]    # 需要插入-的位置
    comp_result = get_compare_result(string_set,pos)
    # for string in comp_result:
    #     print("".join(string))

    # 按照多数原则决定一个位置上的元素
    correct_string = []
    col_word = []
    for col in range(len(comp_result[0])):
        for line in range(len(comp_result)):
            col_word.append(comp_result[line][col])
        dic = Counter(col_word).most_common(1)
        for word in dic:
            if(word[0] != "-"):
                correct_string.append(word[0])
            break
        col_word = []
    correct_string_str = "".join(correct_string)
    return correct_string_str

def delestr(string):
    del_str_list = [str for str in list(string) if(str!="-" and str!=" " and str!="\n")]
    del_str = "".join(del_str_list)
    return del_str

def get_csv(com_clus_set):

    csv_set = []

    for clus_elem in com_clus_set:
        wait_correct = []
        for com_elem in clus_elem.clus_list:
            if(com_elem.new_read != ""):    # 去除待比对字符串列表中为空的元素，因为new_read里有的为""
                wait_correct.append(delestr(com_elem.new_read))
        correct_str = self_correct(wait_correct)
        csv_set.append(csv(clus_elem.left_bp,clus_elem.right_bp,correct_str))

    return csv_set





