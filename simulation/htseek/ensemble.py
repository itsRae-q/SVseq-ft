#!usr/bin/env python
# coding=utf-8

from itertools import islice

class result_info:

    def __init__(self,line_info):
        # self.id_num = line_info[0]  #id号
        if(line_info[1]=="Ref"):
            self.indel_type = line_info[1]+line_info[2]  # 变异类型
            self.res_del_pos = int(line_info[3])  # 删除的开始位置
            self.res_del_len = int(line_info[4])  # 删除的长度
            self.res_one_pos = int(line_info[7])  # 片段一在ref原来的位置
            self.res_one_len = int(line_info[8])  # 片段一的长度
            self.res_two_pos = int(line_info[9])  # 片段二在ref原来的位置
            self.res_two_len = int(line_info[10])  # 片段二的长度，仿真结果有误
        else:
            self.indel_type = line_info[1]   #变异类型
            self.res_del_pos = int(line_info[2]) #删除的开始位置
            self.res_del_len = int(line_info[3])  #删除的长度
            self.res_one_pos = int(line_info[6])  #片段一在ref原来的位置
            self.res_one_len = int(line_info[7])  #片段一的长度
            self.res_two_pos = int(line_info[8])  #片段二在ref原来的位置
            self.res_two_len = int(line_info[9])  #片段二的长度，仿真结果有误

class simu_info:

    def __init__(self,type,del_pos,del_len):
        self.type = type
        self.del_pos = int(del_pos)
        self.del_len = int(del_len)

def get_file_info(filename):
    result_set = []  #存放一个文件的结果信息，每个元素是一个result_info类
    try:
        fin = open(filename)
    except OSError as reason:
        print(filename+"don't exist")
        return result_set
    else:
        for elem in islice(fin,1,None):   # 跳过第一行
            result_set.append(result_info(elem.split()))
        return result_set

def combine_two_file(file_one,file_two):
    add_new_set = [] #存放file_one中有但是file_two中没有的
    for elem in file_one:
        minus_list = [abs(x.res_del_pos-elem.res_del_pos) for x in file_two] # file_one里的该元素减去file_two里所有元素差值的列表
        min_num = min(minus_list)
        min_index = minus_list.index(min_num)  # 与该元素差值最小的file_two中元素的下标
        if(min_num<400 ):#and elem.indel_type==file_two[min_index].indel_type
            continue
        else:
            add_new_set.append(elem)

    return file_two+add_new_set

def get_ensembel(file_set):
    temp_file = file_set[0]
    for file_elem in file_set[1:]:
        temp_file = combine_two_file(file_elem,temp_file)
    return temp_file


def get_file_set():
    file_set = []   #把k较小的放在前面结果的断点跟准确

    # file_set.append(get_file_info("5+50+6"))
    if(get_file_info("5+50+6")!=[]):
        file_set.append(get_file_info("5+50+6"))

    for p_num in [5,4,2]: #4,2
        for k_num in [50,70,100]:#,70,100
            for l_num in [7,8,9,10,11,12,13]:
                # if(p_num==2 and k_num==100 and l_num==7):
                #     continue
                # if (p_num == 4 and k_num == 50 and l_num == 9):
                #     continue
                # if (p_num == 5 and k_num == 50 and l_num == 9):
                #     continue
                if(get_file_info(str(p_num)+"+"+str(k_num)+"+"+str(l_num))!=[]):
                    file_set.append(get_file_info(str(p_num)+"+"+str(k_num)+"+"+str(l_num)))

    # file_set.append(get_file_info("2+50+6"))
    if (get_file_info("2+50+6") != []):
        file_set.append(get_file_info("2+50+6"))

    return file_set

def write_final(final):

    final_sort = []
    for elem in final:
        final_sort.append((elem.res_del_pos, elem.indel_type, elem.res_del_len, elem.res_one_pos, elem.res_one_len,elem.res_two_pos, elem.res_two_len))
    final_sort.sort()
    fout = open("final_result", "w")
    fout.write(str(len(final_sort)) + "\n")
    print(final_sort)
    for elem in final_sort:
        fout.write(str(elem[1]) + " " + str(elem[0]) + " " + str(elem[2]) + " " + str(elem[3]) + " " + str(elem[4]) + " " + str(elem[5]) + " " + str(elem[6]) + "\n")
    fout.close()

def get_simu_info():
    fin = open("indelinformation.txt")
    simu_info_set = []
    for line in fin:
        if(line.split()[0]=="Type:"):
            simu_type = line.split()[1]
        if(line.split()[0]=="DelPos:"):
            simu = simu_info(simu_type,line.split()[1],line.split()[3])
            simu_info_set.append(simu)
    return simu_info_set

def compare_simu_final(simu_info_set,final):
    correct_num = 0
    correct_elem_set = []
    wrong_elem_set = []
    for elem in final:
        minus_list = [abs(x.del_pos-elem.res_del_pos) for x in simu_info_set]
        min_num = min(minus_list)
        min_index = minus_list.index(min_num)  # 与该元素差值最小的file_two中元素的下标
        if(len(elem.indel_type.split("&"))==2 and (elem.indel_type.split("&")[1]=="RefR" or elem.indel_type.split("&")[1]=="RefF")):
            temp_type = elem.indel_type.split("&")[1]+"&"+elem.indel_type.split("&")[0]
        else:
            temp_type = elem.indel_type
        if(min_num<50 and temp_type==simu_info_set[min_index].type):
            correct_elem_set.append(elem)
        else:
            wrong_elem_set.append(elem)
    return correct_elem_set,wrong_elem_set




file_set = get_file_set()
final = get_ensembel(file_set)
write_final(final)
simu_info_set = get_simu_info()
correct_set, wrong_set = compare_simu_final(simu_info_set,final)

fcorrect = open("correct","w")
correct_sort = []
for elem in correct_set:
    correct_sort.append((elem.res_del_pos, elem.indel_type))
correct_sort.sort()
fcorrect.write(str(len(correct_sort))+"\n")
for elem in correct_sort:
    fcorrect.write(elem[1]+" "+str(elem[0])+"\n")
fcorrect.close()

fwrong = open("wrong","w")
wrong_sort = []
for elem in wrong_set:
    wrong_sort.append((elem.res_del_pos, elem.indel_type))
wrong_sort.sort()
fwrong.write(str(len(wrong_sort))+"\n")
for elem in wrong_sort:
    fwrong.write(elem[1]+" "+str(elem[0])+"\n")
fwrong.close()



# pos = []

# for elem in final:
#     pos.append(elem.res_del_pos)
# print(len(final))
# pos.sort()
# for i in pos:
#     print(i)
#
# a = set(pos)
# print(len(a))

# for elem in final:
#     if(elem.res_del_pos==81023):
#         print(elem.indel_type)
# for elem in final:
#     print(elem.res_del_pos,elem.indel_type)
# def get_simu_info():
#     fin = open("indelinformation.txt")
#
#     count = 1
#     for line in fin:
#         print(line.split())
#         count+=1
#         if(count>17):
#             break

    # for line in fin:
    #     line_info = [0]
    #     if(line.split()[0]=="Type"):
    #         line_info.append(line.split()[1])
    #     if (line.split()[0] == "Type"):
    #         line_info.append(line.split()[1])






# for elem in file_set:
#     for elem2 in elem:
#         if(elem2.res_del_pos==441720):
#             print(elem2.res_del_pos)
# print("end")

