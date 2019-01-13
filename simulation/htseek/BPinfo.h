# pragma once
#include<string>
#include<fstream>
using namespace std;
class BPinfo{
public:
	unsigned   LeftBP = 0;
	unsigned   RightBP = 0;
	bool isCombination= 0;
	char SEG1is=0;
	bool SEG1isR=0;
	char SEG2is = 0;
	bool SEG2isR = 0;
	unsigned SEG1startpos = 0;
	unsigned SEG1endpos = 0;
	unsigned SEG2startpos = 0;
	unsigned SEG2endpos = 0;
	string str1;
	string str2;
	unsigned   MidBP = 0;
	string type;
	void printfile(ofstream& fout,long countN){
		
		fout << "Type:" << type << endl;
		fout << "breakpoint:" << LeftBP << "-" << RightBP << endl;
		if (isCombination)
		{
			fout << "segment1:" << SEG1startpos+countN << "-" << SEG1endpos+countN << endl;
			fout << "segment2:" << SEG2startpos+countN << "-" << SEG2endpos+countN << endl;
		}
		else
		{
			fout << "segment:" << SEG1startpos+countN << "-" << SEG1endpos+countN << endl;
		}
		fout << "-----------------------------------------------------" << endl;
	}
};
