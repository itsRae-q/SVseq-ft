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
//字符串的split函数
void split(IN string const& strSrc, IN string const& strDelimiters, OUT vector<string >& vctDst);

//读取配置文件(ref文件名.,倒置变异率,类型比例,分组信息)
void readconfig(OUT vector<configinfo> &conVec, OUT string  &refFile, IN string configName = "config.txt");

//读取ref
void readreference(IN string refName, OUT string&  ref);

//构建cindel,并且存储deletion的信息
void constructcindel(IN vector<configinfo> &conVec, IN string& ref, OUT vector<cindel*>& cindelVec);

//构造不同类型的插入信息
void constructinsertlinfo(IN OUT vector<cindel*>& cindelVec, IN string& ref);
void inREF(IN OUT cindel* &p);
void in5AND3(IN OUT cindel* &p, IN unsigned max, IN vector<cindel*>& cindelVec);
void in5OR3(IN OUT cindel* &p, IN unsigned max, IN vector<cindel*>& cindelVec);
void inRefand5OR3(IN OUT cindel* &p, IN unsigned max, IN vector<cindel*>& cindelVec);

//添加序列信息
void addseqinfo(IN OUT vector<cindel*>& cindelVec, IN  string& ref);

//构造ref
void generateref(IN vector<cindel*>& cindelVec, IN  OUT string& ref);

//将仿真配置写入文件中
void writeinfo(IN vector<cindel*>& cindelVec, IN  string& ref, IN string refsimfile = "refsim.fa", string cindelinfofile = "indelinformation.txt");

//执行步骤
void process(string configname = "config.txt");

void reversecom(string&  str);
#endif

