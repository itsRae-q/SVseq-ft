#include"cindel.h"
double cindel::invertrate = 0.0;
double cindel::refrate = 0.0;
double cindel::_5OR3rate = 0.0;
double cindel::_5AND3rate = 0.0;
double cindel::refand5OR3rate = 100.0;
char cindel::gettype(double d){
	if (0 <= d &&d< refrate)return REF;
	if (refrate <= d&&d < _5OR3rate)return _5OR3;
	if (_5OR3rate <= d&&d < _5AND3rate)return _5AND3;
	if (_5AND3rate <= d&&d <= 100)return REFAND53;
	return 0;
}
void cindel::reversecomplement(string&  str){
	unsigned int i;
	for (i = 0; i < str.size(); ++i)
	{
		switch (str[i])
		{
		case 'A': str[i] = 'T'; break;
		case 'C': str[i] = 'G'; break;
		case 'T': str[i] = 'A'; break;
		case 'G': str[i] = 'C'; break;
		case 'a': str[i] = 'T'; break;
		case 'c': str[i] = 'G'; break;
		case 't': str[i] = 'A'; break;
		case 'g': str[i] = 'C'; break;
		}
	}
	reverse(str.begin(), str.end());
}
bool cindel::isdeletionRange(string::size_type pos, string::size_type size, string::size_type max){
	string::size_type tmp;
	if (pos > max - size)return false;//即选择的范围越界
	tmp = delpos - 50;
	if (delpos<50)	tmp = 0;
	bool a = (pos + size) < tmp;
	tmp = delpos + delsize + 50;
	bool b = pos>tmp;
	return a || b;
}