#!usr/bin/env python
# coding=utf-8

from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import numpy as np
import read

# 参数eps,  min_samples

class posinfo:

    def __init__(self,start_pos,end_pos):
        self.start_pos = start_pos
        self.end_pos = end_pos

# 需要参数：1.read_id 2.该行中一个变异区域的变异位置信息 3.该行的read信息
class indelinfo:

    def __init__(self,read_id,pos_info,line_read):
        self.read_id = read_id    # read_id是从1开始的
        self.start_pos = pos_info.start_pos   # 第一个碱基位置为0，以此类推
        self.end_pos = pos_info.end_pos
        self.read_string = line_read.qstring[self.start_pos:self.end_pos+1]
        self.pattern = line_read.pattern[self.start_pos:self.end_pos + 1]
        self.ref_string = line_read.tstring[self.start_pos:self.end_pos + 1]
        self.ref_start = line_read.start
        self.ref_end = line_read.end
        # self.is_com_indel = None





# 计算一个read的聚类结果,返回的是一个posinfo类
def get_line_cluster(line_var_pos):

    line_cluster_mid = DBSCAN(eps=20, min_samples=12).fit_predict(np.array(line_var_pos).reshape(-1, 1))   # 0.85  eps=20, min_samples=12    0.9可能是30,12      21   14  148  83
    num_cluster = 1
    for num in line_cluster_mid[::-1]:
        if(num != -1):
            num_cluster = num + 1   #聚类个数
            break

    line_pos = []
    if(num_cluster == 1):  #只有一类
        line_pos.append(posinfo(0,len(line_var_pos)))
    else:
        for cluster in range(num_cluster):
            cluster_index = np.argwhere(line_cluster_mid == cluster)
            start_index, end_index = int(cluster_index[0]),int(cluster_index[-1])
            start_pos, end_pos = line_var_pos[start_index], line_var_pos[end_index]
            line_pos.append(posinfo(start_pos,end_pos))
            # print(start_pos,end_pos)
    return line_pos

# 将sam文件信息存放到total_var_pos列表中并返回，列表中的每个元素是一个indelinfo类
def get_indel_info(readInfo):
    # filename = "result.sam"
    # readInfo, reflength = read.readsam(filename)

    total_var_pos = []  # sam文件中所有未匹配点的集合  二维列表
    read_id = 0
    for read_elem in readInfo:
        # print(read_elem.pattern)
        read_id += 1
        line_var_pos = []    # 每一个read中未匹配点的位置
        for num in range(len(read_elem.pattern)):
            if(read_elem.pattern[num] == "*"):
                line_var_pos.append(num)
        line_pos = get_line_cluster(line_var_pos)
        # total_var_pos.append(line_pos)
        for num in range(len(line_pos)):
            total_var_pos.append(indelinfo(read_id,line_pos[num],read_elem))
    return total_var_pos



