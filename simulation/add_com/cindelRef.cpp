#include"cindelRef.h"
void cindelRef::setinvertioninfo(unsigned int pos, string::size_type size, bool isinvert, unsigned int pos2=0, string::size_type size2=0, bool isinvert2=0){
	inpos = pos;
	insize = size;
	isinvertion = isinvert;
	if (isinvertion)
		type = "Ref R";
	else
		type = "Ref F";
}
void cindelRef::setseq(string& ref)
{
	delseq = ref.substr(delpos, delsize);
	inseq = ref.substr(inpos, insize);
	if (isinvertion)reversecomplement(inseq);
	offset =insize - delsize;
}
void cindelRef::generatorcindel(string& ref, int & offset){
	auto newdelpos = delpos + offset;
	ref.replace(newdelpos,delsize,inseq);
	
}
void cindelRef::writefile(ofstream& fout){
	fout << CUT << endl;
	fout << "Type: " << this->getType() << endl;
	fout << "DelPos:" << this->delpos << "	" << "DeletionSize:" << this->delsize << endl;
	fout << this->delseq << endl;
	fout << "Insertion Pos:" << this->inpos << "	" << "Insertion  Size:" << this->insize << "	";
	if (isinvertion)
		fout << "invert" << endl;
	else
		fout << "noninvert" << endl;
	fout << this->inseq << endl;

}