#ifndef __CINDELREFAND5OR3__
#define __CINDELREFAND5OR3__
#include"cindel.h"
class cindelRefand5OR3 :public cindel
{
public:
	cindelRefand5OR3(unsigned int delPos, string::size_type delSize) :cindel(delSize, delPos),type("Refand5OR3"){}
	void setinvertioninfo(unsigned int rpos, string::size_type rsize, bool isrinvert, unsigned int pos, string::size_type size, bool isinvert);
	void setseq(string& ref);
	string getType(){ return type; }
	void generatorcindel(string& ref,  int & offset);
	void writefile(ofstream& out);
private:
	string refseq;
	string::size_type refsize;
	unsigned int  refpos;
	bool isrefinvertion;
	string _53seq;
	string::size_type _53size;
	unsigned int  _53pos;
	bool is53invertion;
	string type;
};
#endif