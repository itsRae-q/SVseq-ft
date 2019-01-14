#!usr/bin/env python
# coding=utf-8

from get_indel_region import *
from extract_feature import *
from svm import *
from cluster_indel import *
import read
from compare_indel import *
from self_correct import *
# from analyse_result import *
import sys, getopt



import matplotlib.pyplot as plt
import numpy as np
import datetime


# def write_total_score(total_score):
#     fout = open("win_size=10","w")
#     for line_score in total_score:
#         for elem in line_score:
#             fout.write(str(elem)+" ")
#         fout.write("\n")
#     fout.close()

# filename = "HG00731.pacbio-blasr-grch38-reheader.20180102.chr22_change.sam"
# readInfo = read.readsam(filename)
# # total_score = get_indel_region(readInfo)
# # write_total_score(total_score)
# print("read done")
# scorename = "win_size=10_score_chr22"
# total_indel_set = get_filter_result(scorename,readInfo)
#
# print(len(total_indel_set))
# feature_list_set = get_indel_feature(total_indel_set)
#
# print("get_feature")
# com_indel_set = get_com_indel(feature_list_set,total_indel_set,readInfo)
#
# print("get_com")
# com_clus_set = get_com_cluster(com_indel_set)
#
# print("get_clus")
# get_compare_indel(com_clus_set,readInfo)
# print(len(com_clus_set))
# count = 0
# pos_set = []
#
# for elem in com_clus_set:
#     pos_elem = (elem.left_bp,elem.right_bp-elem.left_bp)
#     pos_set.append(pos_elem)
#
# def takeFirst(elem):
#     return elem[0]
#
# pos_set.sort(key=takeFirst)
#
# print(pos_set)
#
# fout = open("pos_10-5-70","w")
# # cor = open("correct","w")
# # wro = open("wrong","w")
#
# for elem in pos_set:
#     fout.write(str(elem[0])+"   "+str(elem[1])+"\n")
#     # if(elem[1]<350 or elem[1]>1150):
#     #     wro.write(str(elem[0])+"   "+str(elem[1])+"\n")
#     # else:
#     #     cor.write(str(elem[0])+"   "+str(elem[1])+"\n")
# fout.close()

def main(argv):
   thre = ''
   var = ''
   min_support = ''  # 复杂indel的最小支持度
   filename = ''
   scorename = ''
   out_file = ''
   try:
      opts, args = getopt.getopt(argv,"ht:v:m:i:s:o:",["thre_value=","var_value=","min_support=","sam_file=","score_file=","out_file="])
   except getopt.GetoptError:
      print('connect.py -t <thre_value> - v <var_value> -m <min_support> - i <sam_file> - s <score_file> -o <out_file>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('connect.py -t <thre_value> - v <var_value> -i <sam_file> -s <score_file> -o <out_file>')
         sys.exit()
      elif opt in ("-t", "--thre_value"):
         thre = arg
      elif opt in ("-v", "--var_value"):
         var = arg
      elif opt in ("-m", "--min_support"):
         min_support = arg
      elif opt in ("-i", "--sam_file"):
         filename = arg
      elif opt in ("-s", "--score_file"):
         scorename = arg
      elif opt in ("-o", "--out_file"):
         out_file = arg


   thre_val = float(thre)
   var_val = int(var)
   readInfo = read.readsam(filename)
   # total_score = get_indel_region(readInfo)
   # write_total_score(total_score)
   print("read done")

   starttime = datetime.datetime.now()
   total_indel_set = get_filter_result(scorename, readInfo,thre_val,var_val)
   endtime_filter = datetime.datetime.now()
   print("filter_time:")
   print((endtime_filter - starttime).seconds)
   print("total_indel_len:")
   print(len(total_indel_set))
   # feature_list_set = get_indel_feature(total_indel_set)
   #
   # print("get_feature")
   com_indel_set = get_com_indel(total_indel_set, readInfo)
   print("com_indel_num:")
   print(len(com_indel_set))
   print("get_com")
   com_clus_set = get_com_cluster(com_indel_set,int(min_support))

   print("get_clus")
   get_compare_indel(com_clus_set, readInfo)
   print(len(com_clus_set))
   count = 0
   pos_set = []

   for elem in com_clus_set:
       pos_elem = (elem.left_bp, elem.right_bp - elem.left_bp)
       pos_set.append(pos_elem)

   def takeFirst(elem):
       return elem[0]

   pos_set.sort(key=takeFirst)

   print(pos_set)

   fout = open(out_file,"w")
   # cor = open("correct","w")
   # wro = open("wrong","w")

   for elem in pos_set:
      fout.write(str(elem[0]) + "   " + str(elem[1]) + "   " + str(elem[0] + elem[1]) + "\n")
       # if(elem[1]<350 or elem[1]>1150):
       #     wro.write(str(elem[0])+"   "+str(elem[1])+"\n")
       # else:
       #     cor.write(str(elem[0])+"   "+str(elem[1])+"\n")
   fout.close()





if __name__ == "__main__":
   main(sys.argv[1:])





# print("start_compare!")
#
#
# csv_set = get_csv(com_clus_set)
#
# fout = open("csv_15000","w")
# for elem in csv_set:
#     fout.write(str(elem.left_bp)+" "+str(elem.right_bp)+" "+elem.correct_str+"\n")
# fout.close()
#
# print("success!")
