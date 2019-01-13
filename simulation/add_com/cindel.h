#ifndef __CINDEL__
#define __CINDEL__
#define CUT "-------------------------------------------------------"
#include<string>
#include<fstream>
using namespace std;
const char REF = 'a';
const char _5OR3 = 'b';
const char _5AND3= 'c';
const char REFAND53 = 'd';
class cindel{
public: 
	cindel(string::size_type delsize,unsigned int delpos):delsize(delsize),delpos(delpos){}
	void setdelseq(string seq){ delseq = seq;}
	unsigned int getdelpos()const{ return delpos; }
	string getdelseq(){ return delseq; }
	int getoffset(){ return offset; }
	//判断随机生成的insertion位点是否在deletion范围内
	bool isdeletionRange(unsigned int pos, string::size_type  size, unsigned int max);
	string::size_type getdelsize(){ return delsize; }
	//一系列的纯虚函数
	virtual string getType() = 0;//获得Type值
	virtual void setseq(string&) = 0;//用来保存deletion和insertion的seq
	//用来生成insertion大小和位点
	virtual void setinvertioninfo(unsigned int pos, string::size_type size, bool isinvert, unsigned int pos2, string::size_type size2, bool isinvert2) = 0;
	virtual void writefile(ofstream& out)=0;//写文件
	virtual void generatorcindel(string& ref,  int & offset) = 0;//对ref撒变异
	//输出一个随机值d，根据indel类型比率返回应为何种类型
	static char gettype(double d);
	//静态变量 用来保存倒置变异率以及indel类型比率
	static double invertrate;
	static double refrate;
	static double _5OR3rate;
	static double _5AND3rate;
	static double refand5OR3rate;
protected:
	void reversecomplement(string&  seq);//求反向互补
	string delseq;//
	string::size_type delsize;
	unsigned int delpos;
	int offset;//insertion与deletion的偏移量
};
#endif