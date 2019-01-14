#!usr/bin/env python
# coding=utf-8

from itertools import groupby

# class PatternInfo:
#
#     def __init__(self,pattern_info):
#         self.pattern_list = [''.join(list(g)) for k, g in groupby(pattern_info, key=lambda x: x.isdigit())]


class ReadInfo:

    def __init__(self,lineinfo):
        self.readName = lineinfo[0]
        self.readLength = len(lineinfo[9])
        self.refName = lineinfo[2]
        self.start = int(lineinfo[3])
        self.pattern_list = [''.join(list(g)) for k, g in groupby(lineinfo[5], key=lambda x: x.isdigit())]
        self.qstring, self.pattern, self.tstring, self.ins_num, self.del_num, self.S_val_fir, self.S_val_sec = self.get_rpr(lineinfo[9])
        self.end = self.start + self.readLength + self.del_num - self.ins_num -1 - self.S_val_fir - self.S_val_sec

    def get_rpr(self,read_seq):   # get read,pattern, ref seq

        start_num = 0
        end_num = int(len(self.pattern_list)/2.0)
        S_value_first = 0
        S_value_sec = 0
        if(self.pattern_list[1]=="H"):
            H_value_first=int(self.pattern_list[0])
            start_num += 1
            if(self.pattern_list[3]=="S"):
                S_value_first = int(self.pattern_list[2])
                start_num += 1
        else:
            if(self.pattern_list[1]=="S"):
                S_value_first = int(self.pattern_list[0])
                start_num += 1

        if(self.pattern_list[-1]=="H"):
            H_value_sec=int(self.pattern_list[-2])
            end_num -= 1
            if(self.pattern_list[-3]=="S"):
                S_value_sec = int(self.pattern_list[-4])
                end_num -= 1
        else:
            if(self.pattern_list[-1]=="S"):
                S_value_sec = int(self.pattern_list[-2])
                end_num -= 1

        pattern = ""
        read = ""
        ref = ""
        del_num=0
        ins_num=0
        read_copy_pos = S_value_first
        for num in range(start_num, end_num):
            num_val = int(self.pattern_list[2*num])
            str_val = self.pattern_list[2*num+1]
            if(str_val=="="):
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
            elif(str_val=="X"):
                pattern += num_val*"*"
                read += read_seq[read_copy_pos:read_copy_pos+num_val]
                ref += num_val*"N"
                read_copy_pos += num_val
            else:
                print("UNKNOWN str_val")
        return read,pattern,ref,ins_num,del_num,S_value_first,S_value_sec



def read_sam(filename):
    fin=open(filename)
    readlist=[]
    for line in fin:
        if(line.split()[1]=="0" or line.split()[1]=="16"):
            readinfo = ReadInfo(line.split())
            readlist.append(readinfo)
    fin.close()
    return readlist

def change_sam(filename):
    fin = open(filename)
    fout = open("HG00731.pacbio-blasr-grch38-reheader.20180102.chr22_change.sam","w")
    for line in fin:
        if (line.split()[1] == "0"):
            readinfo = ReadInfo(line.split())
            fout.write(readinfo.readName+" "+str(readinfo.readLength)+" "+readinfo.refName+" "+str(readinfo.start)+" "+str(readinfo.end)+" "+readinfo.qstring+" "+readinfo.pattern+" "+readinfo.tstring+"\n")
    fout.close()
    fin.close()



change_sam("HG00731.pacbio-blasr-grch38-reheader.20180102.chr22.sam")
