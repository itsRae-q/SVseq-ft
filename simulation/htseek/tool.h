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
	static unsigned K ;//��SNP������ڽӾ��룬ͬʱҲ�Ǳȶ�ʱ��query�޼��Ļ���
	static unsigned XCUT ;//��query�޼��ı�����XCUT*K����ÿ���޼��ļ������
	static unsigned SUPPORTNUM ;//����indel�е�SNP������֧����
	static double BIASRATIO ;//����indel�е�INS��DEL�����ƫ����
	static double SNPRATIO ;//����indel�е�SNP����ͱ���
	static unsigned KMERLENGTH ;//kmer�ĳ���
	static unsigned SNP ;//ÿ100BP������Ĳ����Լ�SNP�쳣����ֵ
	static unsigned KSNP;//��ÿ��kmer�У�ƽ��ÿKSNP��BP������1��SNP
	static unsigned reflength;//ref����
	static string output;
};



/*��ȡref�ļ�,ע���һ����ͷ��ϢҪ����*/
long readRef(IN OUT string& ref, IN string refilename = "ref.txt");
/*���str�ķ��򻥲�����*/
void reversecom(IN string& str); 

void write_csv(vector<csv> &clusterVct, string csv_file);

/*��ÿ��csv����mapstring���������BPinfoVct*/
void getBPinfoVct(IN vector<csv>& clustervCT, OUT vector<BPinfo> BPvct,IN string& ref,time_t& startime,long countN);
/*�����е���Ƭ�������Ʊȶԣ����Ƭ��֮��ı༭����С��SNP���Ǿͱ������Ƭ��*/
void BidirectionalGrowth(IN OUT BPinfo& temp,csv& indel,string & ref);

void Writefile(vector<BPinfo>& BPvct, double& t,long countN, string & ref,vector<string> strvct,string filename = "FINDresult.txt");
//����༭����
int edistance(const string source, const string target);
int shortfindBP(BPinfo& temp, string kmerRead, string leftRead, string rightRead, unsigned kmere);
int longfindBP(BPinfo& temp, string kmerRead, string leftRead, string rightRead, unsigned kmere);
