#ifndef __TOOL__
#define __TOOL__
#include"configinfo.h"
#include"randomgen.h"
#include"cindel.h"
#include"cindel5AND3.h"
#include"cindel5OR3.h"
#include"cindelRef.h"
#include"cindelRefand5OR3.h"
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#define IN
#define OUT
using namespace std;
//�ַ�����split����
void split(IN string const& strSrc, IN string const& strDelimiters, OUT vector<string >& vctDst);

//��ȡ�����ļ�(ref�ļ���.,���ñ�����,���ͱ���,������Ϣ)
void readconfig(OUT vector<configinfo> &conVec, OUT string  &refFile, IN string configName = "config.txt");

//��ȡref
void readreference(IN string refName, OUT string&  ref);

//����cindel,���Ҵ洢deletion����Ϣ
void constructcindel(IN vector<configinfo> &conVec, IN string& ref, OUT vector<cindel*>& cindelVec);

//���첻ͬ���͵Ĳ�����Ϣ
void constructinsertlinfo(IN OUT vector<cindel*>& cindelVec, IN string& ref);
void inREF(IN OUT cindel* &p);
void in5AND3(IN OUT cindel* &p, IN unsigned max, IN vector<cindel*>& cindelVec);
void in5OR3(IN OUT cindel* &p, IN unsigned max, IN vector<cindel*>& cindelVec);
void inRefand5OR3(IN OUT cindel* &p, IN unsigned max, IN vector<cindel*>& cindelVec);

//���������Ϣ
void addseqinfo(IN OUT vector<cindel*>& cindelVec, IN  string& ref);

//����ref
void generateref(IN vector<cindel*>& cindelVec, IN  OUT string& ref);

//����������д���ļ���
void writeinfo(IN vector<cindel*>& cindelVec, IN  string& ref, IN string refsimfile = "refsim.fa", string cindelinfofile = "indelinformation.txt");

//ִ�в���
void process(string configname = "config.txt");

void reversecom(string&  str);
#endif

