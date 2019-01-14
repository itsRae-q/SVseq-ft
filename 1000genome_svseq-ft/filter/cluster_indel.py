#!usr/bin/env python
# coding=utf-8

import numpy as np
# 该类是一个列表，列表中每个元素是com_indel类
class cluster:

    def __init__(self):
        self.clus_list = []
        self.left_bp = 0
        self.right_bp = 0

    def add_com(self,com):
        self.clus_list.append(com)

    def get_final_bp(self):

        left_bp_set = []
        right_bp_set = []
        for com_indel in self.clus_list:
            left_bp_set.append(com_indel.ref_start_pos)
            right_bp_set.append(com_indel.ref_end_pos)

        # self.left_bp = np.min(np.array(left_bp_set))
        # self.right_bp = np.max(np.array(right_bp_set))
        self.left_bp = int(np.mean(np.array(left_bp_set)))
        self.right_bp = int(np.mean(np.array(right_bp_set)))






# def find_cluster(com_set):
#     clus = cluster()
#     clus.add_com(com_set[0])
#     del_num = [0]
#
#     #区间有重合的聚成一类，存放在cluster类中，该类是一个列表，列表中存储着聚成一类的com_indel类
#     for num in range(1,len(com_set)):
#         if(com_set[num].ref_start_pos<com_set[0].ref_end_pos and com_set[num].ref_end_pos>com_set[0].ref_start_pos and abs(com_set[num].length-com_set[0].length)<0.4*com_set[0].length):   # 此处0.4为一个阈值，即同一类indel的长度差异不能太大，方便以后做纠错
#             clus.add_com(com_set[num])
#             del_num.append(num)
#
#     #从com_set中删除聚成一类的com_indel
#     for num in del_num[::-1]:
#         del com_set[num]
#
#     return clus,com_set



# def get_com_cluster(com_indel_set):
#     com_clus_set = []
#     while(com_indel_set!=[]):
#         clus, com_indel_set = find_cluster(com_indel_set)
#         if(len(clus.clus_list) > 1):
#             com_clus_set.append(clus)
#
#     #按照最大范围取最终的左右断点
#     for clus_elem in com_clus_set:
#         clus_elem.get_final_bp()
#
#     return com_clus_set

def get_com_cluster(com_indel_set,min_support):
    com_indel_set.sort(key=lambda com_indel: com_indel.ref_start_pos)
    com_clus_set = []
    for elem in com_indel_set:
        if(com_clus_set==[]):
            clus = cluster()
            clus.add_com(elem)
            com_clus_set.append(clus)
        else:
            if(elem.ref_start_pos<com_clus_set[-1].clus_list[0].ref_end_pos and elem.ref_end_pos>com_clus_set[-1].clus_list[0].ref_start_pos and abs(elem.length-com_clus_set[-1].clus_list[0].length)<0.2*com_clus_set[-1].clus_list[0].length):

                com_clus_set[-1].add_com(elem)
            else:
                if (abs(elem.ref_start_pos - com_clus_set[-1].clus_list[0].ref_start_pos) < 300 and com_clus_set[-1].clus_list[0].length > 2 * elem.length):
                    continue
                if(len(com_clus_set[-1].clus_list)<min_support):  # 支持度小鱼2，删除
                    com_clus_set.pop()
                clus_new = cluster()
                clus_new.add_com(elem)
                com_clus_set.append(clus_new)

    if (len(com_clus_set[-1].clus_list) < min_support):  # 支持度小鱼2，最后一个在循环里没检测到
        com_clus_set.pop()

    for clus_elem in com_clus_set:
        clus_elem.get_final_bp()

    return com_clus_set

