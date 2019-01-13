# pragma once
#define IN 
#define OUT
#include<iostream>
#include<algorithm>
#include<string>
#include<fstream>
#include<vector>
#include<list>
#include<time.h>
#include<sstream>
#include"csv.h"
#include"GenHash.h"
#include"BPinfo.h"
using namespace std;
class parameter{
public:
	static unsigned K ;//对SNP聚类的邻接距离，同时也是比对时对query修剪的基数
	static unsigned XCUT ;//对query修剪的倍数；XCUT*K就是每侧修剪的碱基数量
	static unsigned SUPPORTNUM ;//复杂indel中的SNP的最少支持量
	static double BIASRATIO ;//复杂indel中的INS和DEL的最高偏差量
	static double SNPRATIO ;//复杂indel中的SNP的最低比率
	static unsigned KMERLENGTH ;//kmer的长度
	static unsigned SNP ;//每100BP中容许的测序以及SNP异常的数值
	static unsigned KSNP;//在每个kmer中，平均每KSNP个BP中容许1个SNP
	static unsigned reflength;//ref长度
	static string output;
};



/*读取ref文件,注意第一行是头信息要跳过*/
long readRef(IN OUT string& ref, IN string refilename = "ref.txt");
/*求出str的反向互补序列*/
void reversecom(IN string& str); 

void write_csv(vector<csv> &clusterVct, string csv_file);

/*将每个csv进行mapstring，并构造出BPinfoVct*/
void getBPinfoVct(IN vector<csv>& clustervCT, OUT vector<BPinfo> BPvct,IN string& ref,time_t& startime,long countN);
/*将剪切掉的片段做相似比对，如果片段之间的编辑距离小于SNP我们就保留这个片段*/
void BidirectionalGrowth(IN OUT BPinfo& temp,csv& indel,string & ref);

void Writefile(vector<BPinfo>& BPvct, double& t,long countN, string & ref,vector<string> strvct,string filename = "FINDresult.txt");
//计算编辑距离
int edistance(const string source, const string target);
int shortfindBP(BPinfo& temp, string kmerRead, string leftRead, string rightRead, unsigned kmere);
int longfindBP(BPinfo& temp, string kmerRead, string leftRead, string rightRead, unsigned kmere);
