#include"cindel5AND3.h"
void cindel5AND3::setinvertioninfo(unsigned int pos1, string::size_type size1, bool is1invert, unsigned int pos2, string::size_type size2, bool is2invert){
	if (pos1 > pos2)
	{
		auto tmp = pos2;
		pos2 = pos1;
		pos1 = tmp;
	}
	this->pos1 = pos1;
	this->size1 = size1;
	this->is1invertion = is1invert;
	this->pos2 = pos2;
	this->size2 = size2;
	this->is2invertion = is2invert;
	string tmp = "";
	if (pos1 < delpos)
	{
		tmp = tmp + "5'";
		if (is1invertion)tmp = tmp + "R&";
		else tmp = tmp + "F&";
	}
	else
	{
		tmp = tmp + "3'";
		if (is1invertion)tmp = tmp + "R&";
		else tmp = tmp + "F&";
	}
	if (pos2 < delpos)
	{
		tmp = tmp + "5'";
		if (is2invertion)tmp = tmp + "R";
		else tmp = tmp + "F";
	}
	else
	{
		tmp = tmp + "3'";
		if (is2invertion)tmp = tmp + "R";
		else tmp = tmp + "F";
	}
	type = tmp;
}
void cindel5AND3::setseq(string& ref){
	delseq = ref.substr(delpos, delsize);
    seq1 = ref.substr(pos1, size1);
	if (is1invertion) reversecomplement(seq1);
	delseq = ref.substr(delpos, delsize);
	seq2 = ref.substr(pos2, size2);
	if (is2invertion) reversecomplement(seq2);
	offset = size1 + size2 - delsize;
}
void cindel5AND3::generatorcindel(string& ref,  int & offset){
	string seq = seq1+seq2;
	auto newdelpos = delpos + offset;
	ref.replace(newdelpos, delsize, seq);
}
void cindel5AND3::writefile(ofstream& fout){
	fout << CUT << endl;
	fout << "Type: " << this->getType() << endl;
	fout<<"DelPos:"<<this->delpos<<"	"<<"DeletionSize:" << this->delsize << endl;
	fout << this->delseq << endl;
	fout << "Insertion I Pos:" << this->pos1 << "	" << "Insertion I Size:" << this->size1<<"	";
	if (is1invertion)
		fout << "invert" << endl;
	else
		fout << "noninvert"<<endl;
	fout << this->seq1 << endl;
	fout << "Insertion II Pos:" << this->pos2 << "	" << "Insertion II Size:" << this->size2<<"	";
	if (is2invertion)
		fout << "invert" << endl;
	else
		fout << "noninvert" << endl;
	fout << this->seq2 << endl;

}