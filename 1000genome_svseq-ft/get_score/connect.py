#!usr/bin/env python
# coding=utf-8

from get_indel_region import *
from extract_feature import *
from svm import *
from cluster_indel import *
import read
from compare_indel import *
from self_correct import *
from analyse_result import *
import sys,getopt



import matplotlib.pyplot as plt
import numpy as np

import datetime

starttime = datetime.datetime.now()


# def write_total_score(total_score):
#     fout = open("total_score_10","w")
#     for line_score in total_score:
#         for elem in line_score:
#             fout.write(str(elem)+" ")
#         fout.write("\n")
#     fout.close()

# filename = "HG00731.pacbio-blasr-grch38-reheader.20180102.chr22_change.sam"
# readInfo = read.readsam(filename)
# endtime_read = datetime.datetime.now()
# print("read_sam_time:")
# print((endtime_read - starttime).seconds)
# print(len(readInfo))
# total_score = get_indel_region(readInfo)
#
# endtime_comp = datetime.datetime.now()
# print("compute_time:")
# print((endtime_comp - starttime).seconds)
#
#
# write_total_score(total_score)
# endtime_write = datetime.datetime.now()
# print("write_time:")
# print((endtime_write - starttime).seconds)


def main(argv):
   filename = ''
   w_s = ''
   out_name = ''
   try:
      opts, args = getopt.getopt(argv,"hi:w:o:",["ifile=","win_size=","out_name="])
   except getopt.GetoptError:
      print('connect.py -i <sam_file> - w <win_size> -o <out_name>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('connect.py -i <sam_file> -n <win_size>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         filename = arg
      elif opt in ("-w", "--win_size"):
         w_s = arg
      elif opt in ("-o", "--out_name"):
         out_name = arg


   readInfo = read.readsam(filename)
   endtime_read = datetime.datetime.now()
   print("read_sam_time:")
   print((endtime_read - starttime).seconds)
   print(len(readInfo))
   total_score = get_indel_region(readInfo,int(w_s),filename,out_name)

   endtime_comp = datetime.datetime.now()
   print("compute_time:")
   print((endtime_comp - starttime).seconds)

   # write_total_score(total_score)
   # endtime_write = datetime.datetime.now()
   # print("write_time:")
   # print((endtime_write - starttime).seconds)


if __name__ == "__main__":
   main(sys.argv[1:])
