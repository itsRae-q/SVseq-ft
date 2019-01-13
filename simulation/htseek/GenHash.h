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
	/*�趨100��BP���ֱ���ĸ���i*/
	void setSNP(unsigned i){ SNP = i; }
	void setkmerSNP(unsigned i){ KmerSNP = i; }
	/*���и��ӱȶ�,���ڱȶ���������Ϣ��query
	����ֵ��first��ʾ�Ƿ�ɹ��ȶԳɹ�  second�ǳɹ�֮�����ɵ�BPinfo
	*/
	pair<bool, BPinfo> mappedstring(string query);
	/*���е�Ƭ�αȶ�
	����ֵ��first��ʾ�Ƿ�ɹ��ȶԳɹ�  second�ǳɹ�֮�����ɵ�BPinfo
	*/
	pair<bool, BPinfo> singleQuery(string query);

	
private:
	unsigned long(*md5table)[5] = new  unsigned long[1000][5];
	unsigned long simhash(string s);
	string Ref;
	unsigned hashfunction(string seg);//ʹ�õ���BKDRhash����
	vector<bucket> bucketArray;
	unsigned long long MAXSIZE;
	unsigned kmerLenth;
	unsigned SNP = 2;//ƽ��100Bp�����̶��ٱ���
	unsigned KmerSNP = 1;
	unsigned hammingE=10;
	//���ߺ���
	void splitseed(string query, vector<string>& seedVct);//string�ָ��seed
	void reversecom(string& str);//���򻥲�
	int edistance(const string source, const string target);//��༭����


	//����MF����غ���
	void getLocVct(vector<string>& seedVct, vector<bucket>& LocVct);
	void selectPerfectLoction(vector<bucket>& LocVct, unsigned e);
	void initMF(vector<bucket>& LocVct,unsigned e,vector<MF>& MFvct);
	void growthMF(vector<bucket>& LocVct, vector<bucket>& Allvct, list<MF>& MFVct, vector<string>& seedVct, unsigned e, unsigned kmere);

	//Ѱ���м�ϵ����غ���
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

