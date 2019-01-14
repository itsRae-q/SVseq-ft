#!usr/bin/env python
# coding=utf-8


class ReadInfo:

    def __init__(self, readlist):
        self.readName = readlist[0]
        self.readLength = readlist[1]
        self.refName = readlist[2]
        self.start = int(readlist[3])
        self.end = int(readlist[4])
        self.qstring = readlist[5]
        self.pattern = readlist[6]
        self.tstring = readlist[7]
def readsam(filename):
    fin = open(filename)
    readlist = []
    for line in fin:
        readinfo = ReadInfo(line.split())
        readlist.append(readinfo)
    fin.close()

    return readlist   # lens是ref长度  readlist里每一个元素都是一个ReadInfo的类，该类存储了sam文件中的一行信息

