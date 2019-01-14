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
# def gd(x,std=1,exp=0):
#     left=1/(math.sqrt(2*math.pi)*std)
#     right=math.exp(-math.pow(x-exp,2)/(2*math.pow(std,2)))
#     result = round(left*right,5)
#     return result

# 初始化list
# def get_init_list(len):
#     init_list = []
#     for i in range(len):
#         init_list.append(0)
#     return init_list

# left为中间左边个数，right为右边个数
# def comp_gd(left,right):
#     gd_result = []
#     std = 1
#     step = std*3/win_size
#     for num in range(left+right+1):
#         gd_result.append(gd((num-left)*step,std)*step)  #step = sd*2/k
#     return gd_result



# 计算一个点的范围值
# def get_point_var(point_pos,read_len):
#     point_var_score = get_init_list(read_len)
#
#     if(point_pos < win_size):                # *的左侧长度小于win_size
#         part_score = comp_gd(point_pos,win_size)  # 只是变异周围win_size范围内的值
#         for num in range(point_pos+win_size+1):
#             point_var_score[num] = part_score[num]   # 整个read范围内的值
#
#     if (point_pos > read_len-win_size-1):    # *的右侧长度小于win_size
#         part_score = comp_gd(win_size,read_len-point_pos-1)
#         for num in range(point_pos-win_size,read_len):
#             point_var_score[num] = part_score[num-(point_pos-win_size)]
#
#     if(point_pos >= win_size and point_pos <= read_len-win_size-1):
#         part_score = comp_gd(win_size,win_size)
#         for num in range(point_pos-win_size,point_pos+win_size+1):
#             point_var_score[num] = part_score[num-(point_pos-win_size)]
#
#     return point_var_score


# def write_list(in_list):
#     outfile = open("line_score", "w")
#     for num in in_list:
#         outfile.write(str(num))
#         outfile.write("\n")
#     outfile.close()

# def get_total_score(readInfo):
#     total_score = []
#
#     cir_num = 0
#     for read_elem in readInfo:
#         cir_num += 1
#         point_pos = -1  # 字符在read上的位置
#         line_score_list = []
#         for point in read_elem.pattern:
#             point_pos += 1
#             if (point == "*"):
#                 point_score_list = get_point_var(point_pos, len(read_elem.pattern))
#                 line_score_list.append(point_score_list)
#
#         line_score = list(map(sum, zip(*line_score_list)))  # 有的没四舍五入
#
#         # write_list(line_score)
#         total_score.append(line_score)
#
#         print(cir_num)
#         if(cir_num>5):
#             break
#     return total_score

# def get_filter_score(line_score):
#     total_filter_score = []
#     filter_score = filters.gaussian_filter(line_score,filter_var)
#
#     return total_filter_score

# def threshold_choose(total_filter_score, readInfo):
#     read_id = 0
#     total_var_pos = []
#
#     for (elem, read_elem) in zip(total_filter_score, readInfo):
#         read_id += 1
#         ave = elem.mean()
#         var = elem.var()
#         threshold = ave +  thre_ratio*var
#         site_pos = -1     # site位置从0开始计数
#         start_pos, end_pos = -1, -1
#         line_pos = []
#         for site_score in elem:
#             site_pos += 1
#             if(site_score>threshold and start_pos==-1):
#                 start_pos = site_pos
#             if(start_pos!=-1 and site_score<threshold and end_pos==-1):
#                 end_pos = site_pos
#             if(start_pos!=-1 and end_pos!= -1):
#                 line_pos.append(posinfo(start_pos,end_pos))
#                 start_pos, end_pos = -1, -1
#         for num in range(len(line_pos)):
#             total_var_pos.append(indelinfo(read_id,line_pos[num],read_elem))
#     return total_var_pos

def threshold_choose(line_score_filter):

    # ave = line_score_filter.mean()
    # var = line_score_filter.var()
    # threshold = ave + thre_ratio * var
    threshold = thre_ratio
    site_pos = -1  # site位置从0开始计数
    start_pos, end_pos = -1, -1
    line_pos = []
    for site_score in line_score_filter:
        site_pos += 1
        if (site_score > threshold and start_pos == -1):
            start_pos = site_pos
        if (start_pos != -1 and site_score < threshold and end_pos == -1):
            end_pos = site_pos
        if (start_pos != -1 and end_pos != -1):
            line_pos.append(posinfo(start_pos, end_pos))
            start_pos, end_pos = -1, -1
    return line_pos

# def test(total_indel_set):
#     count = 0
#     for elem in total_indel_set:
#         if(abs(elem.start_pos-elem.end_pos)>400 and abs(elem.start_pos-elem.end_pos)<600):
#             count += 1
#     print(len(total_indel_set))
#     print(count)

# def get_indel_region(readInfo):
#
#     global win_size   # 设定变异碱基影响范围 -win_size到 +win_size
#     win_size = 10
#     # global thre_ratio
#     # thre_ratio = 4
#     # global filter_var
#     # filter_var = 30
#
#     total_score = get_total_score(readInfo)   #所有read的变异分数，二维数组
#     # total_filter_score = get_filter_score(total_score)
#     # total_indel_set = threshold_choose(total_filter_score, readInfo)
#
#     return total_score

def get_filter_result(scorename, readInfo,thre_val,var_val):
    global thre_ratio
    thre_ratio = thre_val
    global filter_var
    filter_var = var_val

    fin = open(scorename)
    total_indel_set = []
    read_id = 0
    count_fi = 0
    print("all_num:")
    print(len(readInfo))
    for (line,read_elem) in zip(fin, readInfo):
        read_id += 1
        line_score = []
        for elem in line.split():
            line_score.append(float(elem))
        line_score_filter = filters.gaussian_filter(line_score, filter_var)
        line_pos = threshold_choose(line_score_filter)
        for num in range(len(line_pos)):
            total_indel_set.append(indelinfo(read_id,line_pos[num],read_elem))
        # if(read_elem.start<20950623 and read_elem.start>20940623):
        #     print("find")
        # print(read_id)
        if(count_fi % 1000 == 0):
            print(count_fi)
        count_fi+=1
    return total_indel_set



    # ax2 = plt.subplot(212)
    # kalman1.kalmanfun(line_score)
    #
    # total_score.append(line_score)
    # outfile.write(str(total_score))
    # kalman.kalman_filter(line_score)







