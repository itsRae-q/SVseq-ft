#!usr/bin/env python
# coding=utf-8

#给定左右断点，找对应的read片段
def get_indel(read,left_bp,right_bp):

    count_mis = 0 #记录“-”的个数
    real_left, real_right = 0, 0  # 如果左右断点有一个没找到，则返回为空“”  这里应该是有问题的，应该需要修改最大范围左右断点的选择，考虑最终的左右断点超过该read的开始和结束点的情况
    for num in range(len(read.tstring)):
        if(read.tstring[num] == "-"):
            count_mis += 1
        if(num+read.start-count_mis==left_bp):   # 找加了-的左端点
            real_left = left_bp - read.start + count_mis
        if(num+read.start-count_mis==right_bp):   #找加了-的右断点
            real_right = right_bp - read.start + count_mis
            break
    return read.qstring[real_left:real_right]


def get_compare_indel(com_clus_set,readInfo):

    for clus_elem in com_clus_set:
        for indel_elem in clus_elem.clus_list:
            new_read = get_indel(readInfo[indel_elem.read_id-1],clus_elem.left_bp,clus_elem.right_bp)
            indel_elem.add_read(new_read)