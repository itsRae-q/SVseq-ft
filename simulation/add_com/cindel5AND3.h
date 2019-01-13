#ifndef __CINDEL5AND3__
#define __CINDEL5AND3__
#include"cindel.h"
class cindel5AND3 :public cindel
{
public:
	cindel5AND3(unsigned int delPos, string::size_type delSize) :cindel(delSize, delPos),type("5AND3"){}
	void setinvertioninfo(unsigned int pos1, string::size_type size1, bool is1invert, unsigned int pos2, string::size_type size2, bool is2invert);
	void setseq(string& ref);
	string getType(){ return type; }
	void generatorcindel(string& ref,  int & offset);
	virtual void writefile(ofstream& out);
private:
	string seq1;
	string::size_type size1;
	unsigned int  pos1;
	bool is1invertion;
	string seq2;
	string::size_type size2;
	unsigned int  pos2;
	bool is2invertion;
	string type;
};
#endif