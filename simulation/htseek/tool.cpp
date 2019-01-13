#include"tool.h"
unsigned parameter::K = 100;
unsigned parameter::XCUT = 1;//��query�޼��ı�����XCUT*K����ÿ���޼��ļ������
unsigned parameter::SUPPORTNUM = 15;//����indel�е�SNP������֧����
double parameter::BIASRATIO = 0.1;//����indel�е�INS��DEL�����ƫ����
double parameter::SNPRATIO = 0.4;//����indel�е�SNP����ͱ���
unsigned parameter::KMERLENGTH = 10;//kmer�ĳ���
unsigned parameter::SNP = 15;//ÿ100BP������Ĳ����Լ�SNP�쳣����ֵ
unsigned parameter::KSNP = 5;//��ÿ��kmer�У�ƽ��ÿKSNP��BP������1��SNP
unsigned parameter::reflength = 1000000;
string parameter::output="FINDRESULT.txt";


void Writefile(vector<BPinfo>& BPvct,double& t,long countN,string& ref,vector<string> strvct,string filename){
	unsigned  count = 1;
	ofstream fout(filename);
	fout <<"count:"<<BPvct.size()<<"	";
	fout << "time count:" <<t/1000<<"s"<< endl;
	for (int i = 0;i<BPvct.size();++i)
	{
		auto var = BPvct[i];
		fout << count<<"	";
		++count;
		fout << var.type << "	";
		fout << var.LeftBP << "	" << var.RightBP - var.LeftBP + 1 << "	";
		fout << strvct[i] << "	";
		fout << ref.substr(var.LeftBP-countN, var.RightBP - var.LeftBP + 1)<<"	";
		if (var.isCombination)
		{
			fout << var.SEG1startpos+countN << "	" << var.SEG1endpos - var.SEG1startpos + 1 << "	";
			fout << var.SEG2startpos+countN << "	" << var.SEG2endpos - var.SEG2endpos + 1 << "	";
			fout << ref.substr(var.SEG1startpos, var.SEG1endpos - var.SEG1startpos + 1) << "	";
			fout << ref.substr(var.SEG2startpos, var.SEG2endpos - var.SEG2startpos + 1) << "	";
			
		}
		else
		{
			fout << var.SEG1startpos+countN << "	" << var.SEG1endpos - var.SEG1startpos + 1 << "	";
			fout << 0 << "	" << 0<<"	";
			fout << ref.substr(var.SEG1startpos, var.SEG1endpos - var.SEG1startpos + 1) << "	";
			fout << "*" << "	";
		}

		fout <<" "<<endl;
	}
}

long readRef(IN OUT string& ref,IN string refilename){
	ifstream in;
	in.open(refilename);
	if (!in.is_open())
	{
		cerr << "ref file not exist" << endl;
	}
	string temp;
	getline(in, temp);
	while (getline(in, temp))
	{
		ref += temp;
	}
        parameter::reflength=ref.size();
	transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
	long countN = 0;
	while (ref[countN] == 'N'){
		++countN;
	}
	ref.erase(0, countN);
	auto last = ref.end()-1;
	while (*last == 'N'){
		ref.erase(last);
		--last;
	}
	cout << "ref read complete" << endl;
	return countN;
}
void reversecom(string& str){
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


void write_csv(vector<csv> &clusterVct, string csv_file)
{
	ifstream in;
	in.open(csv_file);
	if (!in.is_open())
	{
		cerr << "csv is not exist" << endl;

	}
	string temp;

	while (getline(in, temp))
	{
		istringstream is(temp);
		csv csv_elem;
		is >> csv_elem.breakpointLEFT >> csv_elem.breakpointRIGHT >> csv_elem.querystring;
		clusterVct.push_back(csv_elem);

	}
}



void getBPinfoVct(IN vector<csv>& clusterVct, OUT vector<BPinfo> BPvct, IN string& ref, time_t& startime,long countN)
{
	parameter p;
	GenHash hashtable(ref, p. KMERLENGTH);//HASH����
	hashtable.setSNP(p.SNP);
	hashtable.setkmerSNP(p.KSNP);
	vector<string> strvct;
	for (int icsv = 0; icsv < clusterVct.size(); ++icsv){
	
		string query = clusterVct[icsv].querystring;
		//��query�������ɾ�����Ի�ô���Ƭ��
		if (query.size() < 4 * p.K*p.XCUT)continue;
		query.erase(0, p.K*p.XCUT);
		query.erase(query.size() - p.K*p.XCUT, p.K*p.XCUT);
		if (query.size() <=p.KMERLENGTH)continue;
		pair<bool, BPinfo> mapinfo = hashtable.mappedstring(query);
		if (mapinfo.first==0)continue;
		else
		{
			BPinfo temp = mapinfo.second;
			temp.LeftBP = clusterVct[icsv].breakpointLEFT;
			temp.RightBP = clusterVct[icsv].breakpointRIGHT;
			//�趨type����
			if (temp.isCombination == 0){
				if (temp.SEG1endpos < temp.LeftBP){
					if (temp.SEG1isR == 0)temp.type = "3'F";
					else temp.type = "3'R";
				}
				else if (temp.SEG1startpos>temp.RightBP){
					if (temp.SEG1isR == 0)temp.type = "5'F";
					else temp.type = "5'R";
				}
				else temp.type = "RefR";
			}
			if (temp.isCombination == 1){
				//SEG1
				if (temp.SEG1endpos < temp.LeftBP){
					if (temp.SEG1isR == 0)temp.type = "3'F";
					else temp.type = "3'R";
				}
				else if (temp.SEG1startpos>temp.RightBP){
					if (temp.SEG1isR == 0)temp.type = "5'F";
					else temp.type = "5'R";
				}
				else if (temp.SEG1isR == 0)temp.type = "RefF";
					 else  temp.type = "RefR";
				//SEG2
				if (temp.SEG2endpos < temp.LeftBP){
					if (temp.SEG2isR == 0)temp.type += "&3'F";
					else temp.type += "&3'R";
				}
				else if (temp.SEG2startpos>temp.RightBP){
					if (temp.SEG2isR == 0)temp.type += "&5'F";
					else temp.type += "&5'R";
				}
				else if (temp.SEG2isR == 0)
						 temp.type += "&RefF";
					else temp.type += "&RefR";
		
			}
			//��ÿ��BPinfo����������չ
			BidirectionalGrowth(temp, clusterVct[icsv], ref);
			BPvct.push_back(temp);
			strvct.push_back(clusterVct[icsv].querystring);
		}
	}	
	time_t end_time = clock();
	auto timed=difftime(end_time, startime);
	cout << '\n' << difftime(end_time, startime)<<endl;
	Writefile(BPvct,timed,countN,ref,strvct,parameter::output);

	cout << "complete";
}

void BidirectionalGrowth(IN OUT BPinfo& temp, csv& indel, string& ref){
	parameter p;
	unsigned cutlength = p.XCUT*p.K;
	if (temp.isCombination == 0){
		if (temp.SEG1isR == 0)
		{
			//���
			string refsegment = ref.substr(indel.breakpointLEFT, cutlength);
			string querysegment = indel.querystring.substr(0, cutlength);
			string breakpointsegment = ref.substr(temp.SEG1startpos - cutlength, cutlength);
			if (p.XCUT*p.K / p.KSNP >edistance(refsegment, querysegment))
			{//�������Ƭ�εı��첢���ߣ���ô���ǾͲ�������ϸ�Ļ���
				temp.LeftBP = indel.breakpointLEFT+cutlength;
			}
			else if (p.XCUT*p.K / p.KSNP >edistance(breakpointsegment, querysegment)){
				temp.LeftBP = indel.breakpointLEFT;
				temp.SEG1startpos -= cutlength;
			}
			else
			{//������е�Ƭ�εı����ʹ��ߣ���ô���Ǿ�Ҫ������ϸ������
				auto count = longfindBP(temp, querysegment, refsegment, breakpointsegment, 100);
				temp.LeftBP = indel.breakpointLEFT;
				temp.SEG1startpos -= cutlength;
				temp.LeftBP += count;
				temp.SEG1startpos += count;
			}
			//�Ҳ�
			refsegment = ref.substr(indel.breakpointRIGHT - cutlength + 1, cutlength);
			querysegment = indel.querystring.substr(indel.querystring.size() - cutlength, cutlength);
			breakpointsegment = ref.substr(temp.SEG1endpos + 1, cutlength);
			if (p.XCUT*p.K / p.KSNP > edistance(refsegment, querysegment))
			{//�������Ƭ�εı��첢���ߣ���ô���ǾͲ�������ϸ�Ļ���
				temp.RightBP = indel.breakpointRIGHT-cutlength;
			}
			else if (p.XCUT*p.K / p.KSNP > edistance(breakpointsegment, querysegment)){
				temp.RightBP = indel.breakpointRIGHT;
				temp.SEG1endpos += cutlength;
			}
			else
			{//������е�Ƭ�εı����ʹ��ߣ���ô���Ǿ�Ҫ������ϸ������
				auto count = longfindBP(temp, querysegment, refsegment, breakpointsegment, 100);
				temp.RightBP = indel.breakpointRIGHT - cutlength + count;
				temp.SEG1endpos += count;
			}
			return;
		}
		else
		{
			//���
			string refsegment = ref.substr(indel.breakpointLEFT, cutlength);
			string querysegment = indel.querystring.substr(0, cutlength);
			string breakpointsegment = ref.substr(temp.SEG1endpos + 1, cutlength);
			reversecom(breakpointsegment);
			if (p.XCUT*p.K / p.KSNP > edistance(refsegment, querysegment))
			{//�������Ƭ�εı��첢���ߣ���ô���ǾͲ�������ϸ�Ļ���
				temp.LeftBP = indel.breakpointLEFT + cutlength;
			}
			else if (p.XCUT*p.K / p.KSNP > edistance(breakpointsegment, querysegment)){
				temp.LeftBP = indel.breakpointLEFT;
				temp.SEG1endpos += cutlength;
			}
			else
			{//������е�Ƭ�εı����ʹ��ߣ���ô���Ǿ�Ҫ������ϸ������
				auto count = longfindBP(temp, querysegment, refsegment, breakpointsegment, 100);
				temp.LeftBP = indel.breakpointLEFT;
				temp.SEG1endpos += cutlength;
				temp.LeftBP += count;
				temp.SEG1endpos -= count;
			}
			//�Ҳ�
			refsegment = ref.substr(indel.breakpointRIGHT - cutlength + 1, cutlength);
			querysegment = indel.querystring.substr(indel.querystring.size() - cutlength, cutlength);
			breakpointsegment = ref.substr(temp.SEG1startpos - cutlength, cutlength);
			if (p.XCUT*p.K / p.KSNP > edistance(refsegment, querysegment))
			{//�������Ƭ�εı��첢���ߣ���ô���ǾͲ�������ϸ�Ļ���
				temp.RightBP = indel.breakpointRIGHT + cutlength;
			}
			else if (p.XCUT*p.K / p.KSNP > edistance(breakpointsegment, querysegment)){
				temp.RightBP = indel.breakpointRIGHT;
				temp.SEG1startpos -= cutlength;
			}
			if (p.XCUT*p.K / p.KSNP < edistance(refsegment, querysegment))

			{//�������Ƭ�εı��첢���ߣ���ô���ǾͲ�������ϸ�Ļ���
				temp.RightBP = indel.breakpointRIGHT;
				temp.SEG1startpos -= cutlength;
			}
			else
			{//������е�Ƭ�εı����ʹ��ߣ���ô���Ǿ�Ҫ������ϸ������
				auto count = longfindBP(temp, querysegment, refsegment, breakpointsegment, 100);
				temp.RightBP = indel.breakpointRIGHT - cutlength + count;
				temp.SEG1startpos -= cutlength;
				temp.SEG1startpos += count;
			}
			return;
		}
	}
	//���������Ƭ�ν������
	if (temp.isCombination == 1)
	{
		if (temp.SEG1isR == 0)
		{
			string refsegment = ref.substr(indel.breakpointLEFT+1, cutlength);
			string querysegment = indel.querystring.substr(0, cutlength);
			string breakpointsegment = ref.substr(temp.SEG1startpos - cutlength, cutlength);
			if (p.XCUT*p.K / p.KSNP > edistance(refsegment, querysegment))
			{//�������Ƭ�εı��첢���ߣ���ô���ǾͲ�������ϸ�Ļ���
				temp.LeftBP = indel.breakpointLEFT + cutlength;
			}
			else if (p.XCUT*p.K / p.KSNP > edistance(breakpointsegment, querysegment)){
				temp.LeftBP = indel.breakpointLEFT;
				temp.SEG1startpos -= cutlength;
			}
			else
			{//������е�Ƭ�εı����ʹ��ߣ���ô���Ǿ�Ҫ������ϸ������
				auto count = longfindBP(temp, querysegment, refsegment, breakpointsegment, 100);
				temp.LeftBP = indel.breakpointLEFT;
				temp.SEG1startpos -= cutlength;
				temp.LeftBP += count;
				temp.SEG1startpos += count;
			}
		}
		else
		{//Ƭ��һ��R�Ļ�
			string refsegment = ref.substr(indel.breakpointLEFT, cutlength);
			string querysegment = indel.querystring.substr(0, cutlength);
			string breakpointsegment = ref.substr(temp.SEG1endpos + 2, cutlength);
			reversecom(breakpointsegment);
			if (p.XCUT*p.K / p.KSNP > edistance(refsegment, querysegment))
			{//�������Ƭ�εı��첢���ߣ���ô���ǾͲ�������ϸ�Ļ���
				temp.LeftBP = indel.breakpointLEFT + cutlength;
			}
			else if (p.XCUT*p.K / p.KSNP > edistance(breakpointsegment, querysegment)){
				temp.LeftBP = indel.breakpointLEFT;
				temp.SEG1endpos += cutlength;
			}
			else
			{//������е�Ƭ�εı����ʹ��ߣ���ô���Ǿ�Ҫ������ϸ������
				auto count = longfindBP(temp, querysegment, refsegment, breakpointsegment, 100);
				temp.LeftBP = indel.breakpointLEFT;
				temp.SEG1endpos += cutlength;
				temp.LeftBP += count;
				temp.SEG1endpos -= count;
			}
		}
		if (temp.SEG2isR == 0)
		{
			string refsegment = ref.substr(indel.breakpointRIGHT - cutlength + 2, cutlength);
			string querysegment = indel.querystring.substr(indel.querystring.size() - cutlength, cutlength);
			string breakpointsegment = ref.substr(temp.SEG2endpos + 1, cutlength);
			if (p.XCUT*p.K / p.KSNP > edistance(refsegment, querysegment))
			{//�������Ƭ�εı��첢���ߣ���ô���ǾͲ�������ϸ�Ļ���
				temp.RightBP = indel.breakpointRIGHT - cutlength;
			}
			else if (p.XCUT*p.K / p.KSNP > edistance(breakpointsegment, querysegment)){
				temp.RightBP = indel.breakpointRIGHT;
				temp.SEG2endpos += cutlength;
			}
			else
			{//������е�Ƭ�εı����ʹ��ߣ���ô���Ǿ�Ҫ������ϸ������
				auto count = longfindBP(temp, querysegment, refsegment, breakpointsegment, 100);
				temp.RightBP = indel.breakpointRIGHT - cutlength + count;
				temp.SEG2endpos += count;
			}
		}
		else{
			string	refsegment = ref.substr(indel.breakpointRIGHT - cutlength , cutlength);
			string	querysegment = indel.querystring.substr(indel.querystring.size() - cutlength-2, cutlength);
			string	breakpointsegment = ref.substr(temp.SEG2startpos - cutlength+2, cutlength);
			reversecom(breakpointsegment);
			if (p.XCUT*p.K / p.KSNP > edistance(refsegment, querysegment))
			{//�������Ƭ�εı��첢���ߣ���ô���ǾͲ�������ϸ�Ļ���
				temp.RightBP = indel.breakpointRIGHT + cutlength;
			}
			else if (p.XCUT*p.K / p.KSNP > edistance(breakpointsegment, querysegment)){
				temp.RightBP = indel.breakpointRIGHT;
				temp.SEG2startpos -= cutlength;
			}
			else
			{//������е�Ƭ�εı����ʹ��ߣ���ô���Ǿ�Ҫ������ϸ������
				auto count = longfindBP(temp, querysegment, refsegment, breakpointsegment, 100);
				temp.RightBP = indel.breakpointRIGHT - cutlength + count;
				temp.SEG2startpos -= cutlength;
				temp.SEG2startpos += count;
			}
		}
		return;
	}
}
int edistance(const string source, const string target){
	
	//step 1  

	int n = source.length();
	int m = target.length();
	if (m == 0) return n;
	if (n == 0) return m;
	//Construct a matrix  
	typedef vector< vector<int> >  Tmatrix;
	Tmatrix matrix(n + 1);
	for (int i = 0; i <= n; i++)  matrix[i].resize(m + 1);

	//step 2 Initialize  

	for (int i = 1; i <= n; i++) matrix[i][0] = i;
	for (int i = 1; i <= m; i++) matrix[0][i] = i;

	//step 3  
	for (int i = 1; i <= n; i++)
	{
		const char si = source[i - 1];
		//step 4  
		for (int j = 1; j <= m; j++)
		{

			const char dj = target[j - 1];
			//step 5  
			int cost;
			if (si == dj){
				cost = 0;
			}
			else{
				cost = 1;
			}
			//step 6  
			const int above = matrix[i - 1][j] + 1;
			const int left = matrix[i][j - 1] + 1;
			const int diag = matrix[i - 1][j - 1] + cost;
			matrix[i][j] = min(above, min(left, diag));

		}
	}//step7  
	return matrix[n][m];

}
int longfindBP(BPinfo& temp, string kmerRead, string leftRead, string rightRead, unsigned kmere){
	unsigned low = 0, high = kmerRead.size() - 1;
	unsigned bppos = 0;
	while (low < kmerRead.size())
	{
		if (leftRead[low] != kmerRead[low])break;
		++low;
	}
	while (high >= 0){
		if (rightRead[high] != kmerRead[high])break;
		--high;
	}
	if (low >= high)
	{
		//λ�÷����˽���,�ϵ�ֱ���ж�
		bppos = high + (low - high) / 2;
		return bppos;
	}
	else//û��������,˵��������SNP����,��С��Χ�����Ǿ�ת������shortfindbp��Ѱ��һ�������Ķϵ�
	{
		bppos = low;
		int t = shortfindBP(temp, kmerRead.substr(low, high - low + 1), leftRead.substr(low, high - low + 1), rightRead.substr(low, high - low + 1), kmere);
		if (t == -1)
		{
			return  t;
		}
		return bppos + t;
	}
}
int	shortfindBP(BPinfo& temp, string kmerRead, string leftRead, string rightRead, unsigned kmere){
	unsigned mine = 50000;
	int p = 0; int i = 0;
	for (; i < kmerRead.size(); ++i){
		string a = kmerRead.substr(0, i);
		string b = kmerRead.substr(i);
		string refa = leftRead.substr(0, i);
		string refb = rightRead.substr(i);
		unsigned e = edistance(a, refa) + edistance(b, refb);
		if (mine > e)
		{
			mine = e; p = i;
		}
	}
	unsigned e = edistance(kmerRead, leftRead);
	if (mine >= e)
	{
		mine = e; p = i;
	}
	if (mine > 2 * kmere) return -1;
	return p;

}
