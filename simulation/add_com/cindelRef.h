#ifndef __CINDELREF__
#define __CINDELREF__
#include"cindel.h"
class cindelRef:public cindel 
{
public:
	cindelRef(unsigned int delPos, string::size_type delSize) :cindel(delSize, delPos), type("REF"){}
	void setinvertioninfo(unsigned int pos, string::size_type size, bool isinvert,unsigned int pos2, string::size_type size2, bool isinvert2);
	void setseq(string& ref);
	string getType(){ return type; }
	void generatorcindel(string& ref, int & offset);
	void writefile(ofstream& fout);
private:
	string inseq;
	string::size_type insize;
	unsigned int  inpos;
	bool isinvertion;
    string type;
};
#endif