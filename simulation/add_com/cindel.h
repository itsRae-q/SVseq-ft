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
	//�ж�������ɵ�insertionλ���Ƿ���deletion��Χ��
	bool isdeletionRange(unsigned int pos, string::size_type  size, unsigned int max);
	string::size_type getdelsize(){ return delsize; }
	//һϵ�еĴ��麯��
	virtual string getType() = 0;//���Typeֵ
	virtual void setseq(string&) = 0;//��������deletion��insertion��seq
	//��������insertion��С��λ��
	virtual void setinvertioninfo(unsigned int pos, string::size_type size, bool isinvert, unsigned int pos2, string::size_type size2, bool isinvert2) = 0;
	virtual void writefile(ofstream& out)=0;//д�ļ�
	virtual void generatorcindel(string& ref,  int & offset) = 0;//��ref������
	//���һ�����ֵd������indel���ͱ��ʷ���ӦΪ��������
	static char gettype(double d);
	//��̬���� �������浹�ñ������Լ�indel���ͱ���
	static double invertrate;
	static double refrate;
	static double _5OR3rate;
	static double _5AND3rate;
	static double refand5OR3rate;
protected:
	void reversecomplement(string&  seq);//���򻥲�
	string delseq;//
	string::size_type delsize;
	unsigned int delpos;
	int offset;//insertion��deletion��ƫ����
};
#endif