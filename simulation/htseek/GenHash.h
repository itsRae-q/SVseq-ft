# pragma once
#include<vector>
#include<list>
#include<cmath>
#include<iostream>
#include<string>
#include<algorithm>
#include<openssl/md5.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"BPinfo.h"
#include"md5table.h"
using namespace std;
#define  min(a,b) ((a<b)?a:b) 
class locbp
{
public:
	locbp(unsigned long long i) :pos(i){};
	unsigned pos;
	bool isperfectmap = 0;
};


class bucket{
public:
	vector<locbp> bucketVct;
};
class finalpair{
public:
	BPinfo bpinfo;
	unsigned simnum;
	unsigned long distance;
	bool operator<(finalpair& right);
};
class MF{
public:
	unsigned MappNum=0;
	unsigned SimNum=0;
	unsigned startSeed=0;
	unsigned EndSeed=0; 
	bool start = 0;
	bool end = 0;
	vector<locbp> locVct;
};

class GenHash{
public:
	GenHash(string &ref, unsigned Kmer);
	~GenHash();
	/*设定100个BP出现变异的个数i*/
	void setSNP(unsigned i){ SNP = i; }
	void setkmerSNP(unsigned i){ KmerSNP = i; }
	/*进行复杂比对,用于比对无先验信息的query
	返回值的first表示是否成功比对成功  second是成功之后生成的BPinfo
	*/
	pair<bool, BPinfo> mappedstring(string query);
	/*进行单片段比对
	返回值的first表示是否成功比对成功  second是成功之后生成的BPinfo
	*/
	pair<bool, BPinfo> singleQuery(string query);

	
private:
	unsigned long(*md5table)[5] = new  unsigned long[1000][5];
	unsigned long simhash(string s);
	string Ref;
	unsigned hashfunction(string seg);//使用的是BKDRhash函数
	vector<bucket> bucketArray;
	unsigned long long MAXSIZE;
	unsigned kmerLenth;
	unsigned SNP = 2;//平均100Bp能容忍多少变异
	unsigned KmerSNP = 1;
	unsigned hammingE=10;
	//工具函数
	void splitseed(string query, vector<string>& seedVct);//string分割成seed
	void reversecom(string& str);//求反向互补
	int edistance(const string source, const string target);//求编辑距离


	//构造MF的相关函数
	void getLocVct(vector<string>& seedVct, vector<bucket>& LocVct);
	void selectPerfectLoction(vector<bucket>& LocVct, unsigned e);
	void initMF(vector<bucket>& LocVct,unsigned e,vector<MF>& MFvct);
	void growthMF(vector<bucket>& LocVct, vector<bucket>& Allvct, list<MF>& MFVct, vector<string>& seedVct, unsigned e, unsigned kmere);

	//寻找中间断点的相关函数
	unsigned simplefindBP(BPinfo& temp, bool direct, string kmerRead, string refRead, unsigned pos);
	int shortfindBP(BPinfo& temp, string kmerRead, string leftRead, string rightRead, unsigned kmere);
	int longfindBP(BPinfo& temp, string kmerRead, string leftRead, string rightRead, unsigned kmere);
	unsigned char_to_unsigned(char c);
	unsigned long finall(char* src);
	void middle(const unsigned char* src, char* des);
	void mymd5(unsigned const char* src, unsigned long* des);
	void init_cache(unsigned long cache[1000][5]);
};
int HammingDistance(unsigned long l1, unsigned long l2);

