#!usr/bin/env python
# coding=utf-8

class simu_info:

    def __init__(self,del_pos,del_len):

        self.del_pos = int(del_pos)
        self.del_len = int(del_len)


def get_simu_info():
    fin = open("indelinformation.txt")
    simu_info_set = []
    for line in fin:
        if(line.split()[0]=="Type:"):
            simu_type = line.split()[1]
        if(line.split()[0]=="DelPos:"):
            simu = simu_info(line.split()[1],line.split()[3])
            simu_info_set.append(simu)
    return simu_info_set

def compare_simu_final(simu_info_set,result):
    correct_num = 0
    pos_threshold = 60
    len_threshold = 0.1
    correct_elem_set = []
    wrong_elem_set = []
    for elem in result:
        minus_list = [abs(x.del_pos-elem.del_pos) for x in simu_info_set]
        min_num = min(minus_list)
        min_index = minus_list.index(min_num)  # 与该元素差值最小的file_two中元素的下标
        if(min_num<pos_threshold and abs(elem.del_len-simu_info_set[min_index].del_len)<elem.del_pos*len_threshold):
            correct_elem_set.append(elem)
        else:
            wrong_elem_set.append(elem)
    return correct_elem_set,wrong_elem_set

def get_result_info():
    fin=open("correct")
    result_set = []
    for elem in fin:
        result_set.append(simu_info(elem.split()[0],elem.split()[1]))
    return result_set


def analyze_result():

    simu = get_simu_info()
    result = get_result_info()
    correct_set, wrong_set = compare_simu_final(simu,result)

    fcorrect = open("correct_info","w")
    fcorrect.write(str(len(correct_set))+"\n")
    for elem in correct_set:
        fcorrect.write(str(elem.del_pos)+" "+str(elem.del_len)+"\n")
    fcorrect.close()

    fwrong = open("wrong_info","w")
    fwrong.write(str(len(wrong_set))+"\n")
    for elem in wrong_set:
        fwrong.write(str(elem.del_pos)+" "+str(elem.del_len)+"\n")
    fcorrect.close()

    print(len(correct_set))
    print("11")


analyze_result()
# fcorrect = open("correct_info","w")
# correct_sort = []
# for elem in correct_set:
#     correct_sort.append(elem.del_pos)
# correct_sort.sort()
# fcorrect.write(str(len(correct_sort))+"\n")
# for elem in correct_sort:
#     fcorrect.write(elem[1]+" "+str(elem[0])+"\n")
# fcorrect.close()
#
# fwrong = open("wrong_info","w")
# wrong_sort = []
# for elem in wrong_set:
#     wrong_sort.append(elem.del_pos)
# wrong_sort.sort()
# fwrong.write(str(len(wrong_sort))+"\n")
# for elem in wrong_sort:
#     fwrong.write(elem[1]+" "+str(elem[0])+"\n")
# fwrong.close()