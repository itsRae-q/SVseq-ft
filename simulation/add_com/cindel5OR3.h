#ifndef __CINDEL5OR3__
#define __CINDEL5OR3__
#include"cindel.h"
class cindel5OR3 :public cindel
{
public:
	cindel5OR3(unsigned int delPos, string::size_type delSize) :cindel(delSize, delPos),type("5OR3"){}
	void setinvertioninfo(unsigned int pos, string::size_type size, bool isinvert, unsigned int pos2, string::size_type size2, bool isinvert2);
	void setseq(string& ref);
	string getType(){ return type; }
	void generatorcindel(string& ref,  int & offset);
	void writefile(ofstream& fout);
private:
	string inseq;
	string::size_type insize;
	unsigned int  inpos;
	bool isinvertion;
	string type;
};
#endif