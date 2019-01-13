#!usr/bin/env python
# coding=utf-8

from get_indel_region import *
from extract_feature import *
from svm import *
from cluster_indel import *
import read
from compare_indel import *
from self_correct import *



import matplotlib.pyplot as plt
import numpy as np

import sys, getopt

def main(argv):
   sam_file = ''
   # win = ''
   # threshold = ''
   # filter_var = ''
   # min_sup = ''
   try:
      opts, args = getopt.getopt(argv,"hs:",["sam_file="])
   except getopt.GetoptError:
      print('connect.py  -s <sam>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('connect.py -s <sam_file>')
         sys.exit()
      elif opt in ("-s", "--sam_file"):
          sam_file = arg

   readInfo, reflength = read.readsam(sam_file)
   print("read_sam_done!")
   win = 10
   threshold = 0.3
   filter_var = 30
   min_sup = 1

   total_indel_set = get_indel_region(readInfo, win, threshold, filter_var)
   # print(len(total_indel_set))
   print("get_indel_done")

   feature_list_set = get_indel_feature(total_indel_set)
   com_indel_set = get_com_indel(feature_list_set, total_indel_set, readInfo)
   com_clus_set = get_com_cluster(com_indel_set, min_sup)
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

   # print(pos_set)

   fout_pos = open(sam_file+"_pos", "w")
   # cor = open("correct","w")
   # wro = open("wrong","w")

   for elem in pos_set:
       fout_pos.write(str(elem[0]) + "   " + str(elem[1]) + "\n")

   fout_pos.close()

   print("start_compare!")

   csv_set = get_csv(com_clus_set,sam_file)



   print(sam_file+"success!")


if __name__ == "__main__":
   main(sys.argv[1:])
