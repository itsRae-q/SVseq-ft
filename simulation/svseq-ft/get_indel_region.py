#!usr/bin/env python
# coding=utf-8

import read
import math
import numpy as np
import matplotlib.pyplot as plt
# import kalman
# import kalman1
import scipy.signal as signal
from scipy.ndimage import filters


class posinfo:

    def __init__(self,start_pos,end_pos):
        self.start_pos = start_pos  # read中变异开始的位置（没有去掉-）
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

# def kalman_filter(line_score):
#
#     x = [x for x in range(len(line_score))]
#     y = line_score
#     y_smoth = signal.medfilt(y,101)
#     # print(x)
#     # print(y_smoth)
#     return y_smoth

# 期望为exp，标准差为std的高斯函数
def gd(x,std=1,exp=0):
    left=1/(math.sqrt(2*math.pi)*std)
    right=math.exp(-math.pow(x-exp,2)/(2*math.pow(std,2)))
    result = round(left*right,5)
    return result

# 初始化list
def get_init_list(len):
    init_list = []
    for i in range(len):
        init_list.append(0)
    return init_list

# left为中间左边个数，right为右边个数
def comp_gd(left,right):
    gd_result = []
    std = 1
    step = std*3/win_size
    for num in range(left+right+1):
        gd_result.append(gd((num-left)*step,std)*step)  #step = sd*2/k
    return gd_result



def get_total_score(readInfo):
    total_score = []
    gd_list = comp_gd(win_size,win_size)
    cir_num = 0
    for read_elem in readInfo:
        cir_num += 1
        point_pos = -1  # 字符在read上的位置
        read_len = len(read_elem.pattern)
        line_score_list = [0]*read_len
        for point in read_elem.pattern:
            point_pos += 1
            if (point == "*" and point_pos>win_size+1 and point_pos<read_len-win_size-1):
                line_score_list[point_pos-win_size:point_pos+win_size+1] = [x+y for x,y in zip(line_score_list[point_pos-win_size:point_pos+win_size+1],gd_list)]
        total_score.append(line_score_list)
        # for elem in line_score_list:
        #     fout.write(str(elem)+" ")
        # fout.write("\n")
        # total_score.append(line_score_list)
        if(cir_num % 1000 == 0):
            print(cir_num)

    return total_score

def get_filter_score(total_score,filter_var):
    total_filter_score = []
    for elem in total_score:
        filter_score = filters.gaussian_filter(elem, filter_var)
        total_filter_score.append(filter_score)

    return total_filter_score

def threshold_choose(total_filter_score, readInfo,threshold):
    read_id = 0
    total_var_pos = []

    for (elem, read_elem) in zip(total_filter_score, readInfo):
        read_id += 1
        # ave = elem.mean()
        # var = elem.var()
        # threshold = ave + thre_ratio*var
        # threshold = 0.25
        site_pos = -1     # site位置从0开始计数
        start_pos, end_pos = -1, -1
        line_pos = []
        for site_score in elem:
            site_pos += 1
            if(site_score>threshold and start_pos==-1):
                start_pos = site_pos
            if(start_pos!=-1 and site_score<threshold and end_pos==-1):
                end_pos = site_pos
            if(start_pos!=-1 and end_pos!= -1):
                line_pos.append(posinfo(start_pos,end_pos))
                start_pos, end_pos = -1, -1
        for num in range(len(line_pos)):
            total_var_pos.append(indelinfo(read_id,line_pos[num],read_elem))
    return total_var_pos

def test(total_indel_set):
    count = 0
    for elem in total_indel_set:
        if(abs(elem.start_pos-elem.end_pos)>400 and abs(elem.start_pos-elem.end_pos)<600):
            count += 1
    print(len(total_indel_set))
    print(count)

def get_indel_region(readInfo,win,threshold,filter_var):

    global win_size
    win_size = win
    # global thre_ratio
    # thre_ratio = 7



    total_score = get_total_score(readInfo)   #所有read的变异分数，二维数组
    total_filter_score = get_filter_score(total_score,filter_var)
    total_indel_set = threshold_choose(total_filter_score, readInfo,threshold)

    return total_indel_set





    # ax2 = plt.subplot(212)
    # kalman1.kalmanfun(line_score)
    #
    # total_score.append(line_score)
    # outfile.write(str(total_score))
    # kalman.kalman_filter(line_score)







