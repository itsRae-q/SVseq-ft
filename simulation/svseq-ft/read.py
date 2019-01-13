#!usr/bin/env python
# coding=utf-8


class ReadInfo:

    def __init__(self, readlist):
        self.readName = readlist[0]
        self.readLength = readlist[1]
        self.refName = readlist[5]
        self.start = int(readlist[7])
        self.end = int(readlist[8])
        self.direct = readlist[9]
        self.qstring = readlist[16]
        self.pattern = readlist[17]
        self.tstring = readlist[18]
        if self.direct == '-':
            self.qstring = self.qstring[::-1]
            self.pattern = self.pattern[::-1]
            self.tstring = self.tstring[::-1]
            temp = list(self.qstring)
            for i in range(len(temp)):
                if temp[i] == 'A' or temp[i] == 'a':
                   temp[i] = 'T'
                elif temp[i] == 'G' or temp[i] == 'g':
                   temp[i] = 'C'
                elif temp[i] == 'C' or temp[i] == 'c':
                   temp[i] = 'G'
                elif temp[i] == 'T' or temp[i] == 't':
                   temp[i] = 'A'
            self.qstring = ''.join(temp)
            temp = list(self.tstring)
            for i in range(len(temp)):
                if temp[i] == 'A' or temp[i] == 'a':
                   temp[i] = 'T'
                elif temp[i] == 'G' or temp[i] == 'g':
                   temp[i] = 'C'
                elif temp[i] == 'C' or temp[i] == 'c':
                   temp[i] = 'G'
                elif temp[i] == 'T' or temp[i] == 't':
                   temp[i] = 'A'
            self.tstring = ''.join(temp)
            # self.start = int(readlist[6]) - int(readlist[8])
            # self.end = int(readlist[6]) - int(readlist[7])


def readsam(filename):
    fin = open(filename)
    readlist = []
    for line in fin:
        readinfo = ReadInfo(line.split())
        if len(readlist) == 0 or readlist[-1].readName != readinfo.readName:
            readlist.append(readinfo)
    fin.close()
    fin = open(filename)
    lens = fin.readline().split()[6]
    return readlist, lens   # lens是ref长度  readlist里每一个元素都是一个ReadInfo的类，该类存储了sam文件中的一行信息

# readsam("result_part.sam")