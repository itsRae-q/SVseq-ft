#include"cindelRefand5OR3.h"
void cindelRefand5OR3::setinvertioninfo(unsigned int rpos, string::size_type rsize, bool isrinvert, unsigned int _pos, string::size_type _size, bool _isinvert){
	refpos = rpos;
	refsize = rsize;
	_53size = _size;
	_53pos = _pos;
	isrefinvertion = isrinvert;
	is53invertion = _isinvert;
	string tmp;
	if (isrinvert)tmp = "Ref R";
	else tmp = "Ref F";
	tmp = tmp + "&";
	if (rpos < _53pos)tmp = tmp + "3'";
	else tmp = tmp + "5'";
	if (_isinvert)tmp = tmp + "R";
	else tmp = tmp + "F";
	type = tmp;
}
void cindelRefand5OR3::setseq(string& ref){
	delseq = ref.substr(delpos, delsize);
	string rseq = ref.substr(refpos, refsize);
	if (isrefinvertion) reversecomplement(rseq);
	refseq = rseq;
	string seq = ref.substr(_53pos, _53size);
	if (is53invertion) reversecomplement(seq);
	_53seq = seq;
	offset = refsize + _53size - delsize;
}
void cindelRefand5OR3::generatorcindel(string& ref,  int & offset){
	string seq = "";
	if (refpos < _53pos)seq = refseq+_53seq;
	else seq = _53seq + refseq;
	auto newdelpos = delpos + offset;
	ref.replace(newdelpos, delsize, seq);
}
void cindelRefand5OR3::writefile(ofstream& fout){
	fout << CUT << endl;
	fout << "Type: " << this->getType() << endl;
	fout << "DelPos:" << this->delpos << "	" << "DeletionSize:" << this->delsize << endl;
	fout << this->delseq << endl;
	fout << "RefInsertion  Pos:" << this->refpos << "	" << "RefInsertion  Size:" << this->refsize << "	";
	if (isrefinvertion)
		fout << "invert" << endl;
	else
		fout << "noninvert" << endl;
	fout << this->refseq << endl;
	fout << "5'OR3'Insertion  Pos" << this->_53pos << "	" << "5'OR3'Insertion  Size" << this->_53size<< "	";
	if (is53invertion)
		fout << "invert" << endl;
	else
		fout << "noninvert" << endl;
	fout << this->_53seq << endl;


}