#!usr/bin/env python
# coding=utf-8

from itertools import groupby
import sys,getopt
from operator import itemgetter, attrgetter

# class PatternInfo:
#
#     def __init__(self,pattern_info):
#         self.pattern_list = [''.join(list(g)) for k, g in groupby(pattern_info, key=lambda x: x.isdigit())]

class Read_Seq:   #一条read信息包含在ref上开始的位置和ciga序列，read序列是公用的，在主read上

    def __init__(self,s_start,s_pattern):
        self.s_start = int(s_start)
        self.s_pattern = [''.join(list(g)) for k, g in groupby(s_pattern, key=lambda x: x.isdigit())]

class ReadInfo:

    def __init__(self,lineinfo):
        self.readName = lineinfo[0]
        self.all_readLength = len(lineinfo[9])
        self.refName = lineinfo[2]
        self.is_reverse = (False if int(lineinfo[1]) % 256 == 0 else True)
        self.is_valid = True   # 标记该条read是否有效，非正常情况的read为无效

        self.all_read = [Read_Seq(lineinfo[3], lineinfo[5])] + self.get_sup_read(lineinfo[11])   # 原始read信息+补充read信息
        self.qstring, self.pattern, self.tstring, self.start, self.end = self.get_rpr(lineinfo[9])  # lineinfo[9]是read序列

        self.readLength = len(self.qstring) - self.qstring.count("-")

    def get_sup_read(self,sup_info):   # 从补充信息中提取read，包含开始位置，ciga序列，是否反转
        sup_read_list = []
        if(sup_info.split(":")[0]=="SA"):   # 若不是SA，则无补充read
            for elem in sup_info.split(";")[:-1]:
                sup_info = elem.split(",")[1:]
                start_pos = int(sup_info[0])
                is_rev = (False if sup_info[1]=="+" else True)
                sup_pattern = sup_info[2]
                if(is_rev==self.is_reverse):
                    sup_read_list.append(Read_Seq(start_pos,sup_pattern))

        return sup_read_list


    def get_rpr(self,read_seq):   # get read,pattern, ref seq  对于第一个左边为S和最后一个右边为S的暂时不补充，对于删除的read比对到原来位置处的先不处理

        self.all_read.sort(key=attrgetter('s_start'))
        pattern = ""
        read = ""
        ref = ""
        start_pos = self.all_read[0].s_start
        end_pos = self.all_read[0].s_start-1
        for read_elem in self.all_read:
            if(read_elem.s_start <= end_pos):
                self.is_valid = False
                break
            if (read_elem.s_start - end_pos > 100000):   # 相距太远的情况不考虑，即有一部分比到其他位置了
                self.is_valid = False
                break
            else:            # 不同read不连续的部分填补*
                read += "U"*(read_elem.s_start-end_pos-1)
                ref += "U"*(read_elem.s_start-end_pos-1)
                pattern += "*"*(read_elem.s_start-end_pos-1)

            start_num = 0
            end_num = int(len(read_elem.s_pattern)/2.0)
            S_value_first = 0
            S_value_sec = 0

            if (read_elem.s_pattern[1] == "H"):
                H_value_first = int(read_elem.s_pattern[0])
                start_num += 1
                if (read_elem.s_pattern[3] == "S"):
                    S_value_first = int(read_elem.s_pattern[2])
                    start_num += 1
            else:
                if (read_elem.s_pattern[1] == "S"):
                    S_value_first = int(read_elem.s_pattern[0])
                    start_num += 1

            if (read_elem.s_pattern[-1] == "H"):
                H_value_sec = int(read_elem.s_pattern[-2])
                end_num -= 1
                if (read_elem.s_pattern[-3] == "S"):
                    S_value_sec = int(read_elem.s_pattern[-4])
                    end_num -= 1
            else:
                if (read_elem.s_pattern[-1] == "S"):
                    S_value_sec = int(read_elem.s_pattern[-2])
                    end_num -= 1

            # if(read_elem.s_pattern[1]=="S"):
            #     S_value_first = int(read_elem.s_pattern[0])
            #     start_num += 1
            # if(read_elem.s_pattern[-1]=="S"):
            #     S_value_sec = int(read_elem.s_pattern[-2])
            #     end_num -= 1

            del_num = 0
            ins_num = 0
            read_copy_pos = S_value_first
            for num in range(start_num, end_num):
                num_val = int(read_elem.s_pattern[2*num])
                str_val = read_elem.s_pattern[2*num+1]
                if(str_val=="M"):
                    pattern += num_val*"|"
                    read += read_seq[read_copy_pos:read_copy_pos+num_val]
                    ref += read_seq[read_copy_pos:read_copy_pos+num_val]
                    read_copy_pos += num_val
                elif(str_val=="I"):
                    pattern += num_val*"*"
                    read += read_seq[read_copy_pos:read_copy_pos+num_val]
                    ref += num_val*"-"
                    read_copy_pos += num_val
                    ins_num += num_val
                elif(str_val=="D"):
                    pattern += num_val*"*"
                    ref += num_val*"N"
                    read += num_val*"-"
                    del_num += num_val
                else:
                    print("UNKNOWN str_val")

            end_pos = read_elem.s_start + self.all_readLength + del_num - ins_num - 1 - S_value_first - S_value_sec

        return read,pattern,ref,start_pos,end_pos



# def read_sam(filename):
#     fin=open(filename)
#     readlist=[]
#     for line in fin:
#         if(line.split()[1]=="0"):
#             readinfo = ReadInfo(line.split())
#             readlist.append(readinfo)
#     fin.close()
#     return readlist

def change_sam(filename,file_out):
    fin = open(filename)
    fout = open(file_out,"w")
    # fout2 = open(filename+"_ori.sam","w")
    for line in fin:
        ciga_val = line.split()[5]
        str_seq = [x for x in ciga_val if not x.isdigit()]
        if ((line.split()[1] == "0" or line.split()[1] == "16" or line.split()[1] == "256" or line.split()[1] == "272") and str_seq[0]!="H" and str_seq[-1]!="H"):
            # fout2.write(line)
            readinfo = ReadInfo(line.split())
            if(readinfo.is_valid==True):
                fout.write(readinfo.readName+" "+str(readinfo.readLength)+" "+readinfo.refName+" "+str(readinfo.start)+" "+str(readinfo.end)+" "+readinfo.qstring+" "+readinfo.pattern+" "+readinfo.tstring+"\n")
    fout.close()
    # fout2.close()
    fin.close()







def main(argv):
   sam_file = ''
   file_out = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","file_out="])
   except getopt.GetoptError:
      print('connect.py -i <sam_file>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('read_sam.py -i <sam_file>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         sam_file = arg
      elif opt in ("-o", "--file_out"):
         file_out = arg

   change_sam(sam_file,file_out)




if __name__ == "__main__":
   main(sys.argv[1:])

