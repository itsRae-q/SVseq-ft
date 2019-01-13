#include"cindeltool.h"
void split(IN string const& strSrc, IN string const& strDelimiters, OUT vector<string >& vctDst)
{
	string::size_type pos1, pos2;
	pos2 = 0;
	while (true)
	{
		pos1 = strSrc.find_first_not_of(strDelimiters, pos2);
		if (pos1 == string::npos) break;
		pos2 = strSrc.find_first_of(strDelimiters, pos1 + 1);
		if (pos2 == string::npos)
		{
			vctDst.push_back(strSrc.substr(pos1));
			break;
		}
		vctDst.push_back(strSrc.substr(pos1, pos2 - pos1));
		++pos2;
	}
}
void readconfig(OUT vector<configinfo> &conVec, OUT string&  refFile, IN string configName ){
	ifstream in;
	in.open(configName);
	if (!in.is_open())
	{
		cerr << "文件打开错误";
	}
	//读取ref文件信息
	getline(in, refFile);
	//读取inversion变异率
	string temp;
	getline(in, temp);
	auto pos = temp.find_first_of(':');
	temp.erase(0, pos + 1);
	cindel::invertrate = stod(temp);
	//读取类型比例
	getline(in, temp);
	pos = temp.find_first_of('=');
	temp.erase(0, pos + 1);
	vector<string> rate;
	split(temp,":", rate);
	cindel::refrate=stod(rate[0]);
	cindel::_5OR3rate = stod(rate[1]) + cindel::refrate;
	cindel::_5AND3rate = stod(rate[2]) + cindel::_5OR3rate;
	if (cindel::_5AND3rate+stod(rate[3])!=100)
	{
		cerr << "错误:类型比例之和不为100" << endl;
	}
	//读取分组indel信息
	string  configInfo;
	getline(in, configInfo);
	
	while (getline(in, configInfo))
	{
		string::size_type  pos = configInfo.find_first_of(',');
		unsigned int low = stoul(configInfo.substr(0, pos));
		configInfo.erase(0, pos + 1);
		pos = configInfo.find_first_of(':');
		unsigned int high = stoul(configInfo.substr(0, pos));
		configInfo.erase(0, pos + 1);
		if (low > high){ cerr << "范围错误"; }
		unsigned int count = stoul(configInfo);
		conVec.push_back(configinfo(low, high, count));
	}
}
void readreference(IN string refName, OUT string&  ref){
	ifstream in;
	in.open(refName);
	if (!in.is_open())
	{
		cerr << "读取ref文件错误,请检查文件名或文件是否存在"<<endl;
	}
	string temp;
	while (getline(in, temp))
	{
		ref += temp;
	}
	transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
	cout << "成功读取ref" << endl;
}
void constructcindel(IN vector<configinfo> &conVec, IN string& ref, OUT vector<cindel*>& cindelVec){
	randomgen randomGen;
	auto max = ref.size();
	cindel* p;
	for (auto ite = conVec.begin(); ite < conVec.end(); ++ite)
	{
		unsigned count = ite->count;
		while (count != 0)
		{
			bool b = true;
			string::size_type pos = 0, size = 0;
			for (size_t i = 0; i < 1000000; ++i)
			{
				b = true;
				pos = randomGen.randompos(ref.size());
				size = randomGen.randomdeletionsize(ite->low, ite->high);
				for each (auto var in cindelVec)
				{
					if (!var->isdeletionRange(pos, size, max)){
						b = false;
						break;
					}
				}
				if (b)break;

			}
			if (b)
			{
				char TYPE=cindel::gettype(randomGen.randomtype());

				switch (TYPE)
				{
				case REF:
					p = new cindelRef(pos, size);
					cindelVec.push_back(p); break;
				case _5OR3:
					p = new cindel5OR3(pos, size);
					cindelVec.push_back(p); break;
				case _5AND3:
					 p = new cindel5AND3(pos, size);
					cindelVec.push_back(p); break;
				case REFAND53:
					 p = new cindelRefand5OR3(pos, size);
					cindelVec.push_back(p); break;
				default:cerr << "构造cindel类型时出错!!"<<endl; break;
				}
				--count;
			}
			else
			{
				cerr << "找不到合适的删除位点请增大ref或者减少indel的个数";
				return;
			}
		}
	}
	cout << "初步构造cindel成功,成功写入了deltion信息" << endl;
}
void constructinsertlinfo(IN OUT vector<cindel*>& cindelVec, IN string& ref)
{
	auto max = ref.size();
	//按照delpos进行排序
	sort(cindelVec.begin(), cindelVec.end(), [](const cindel* a, const  cindel* b){return a->getdelpos() < b->getdelpos(); });
	//根据类型信息对不同类型选择不同的插入方法;
	string TYPE = "";
	auto MAX = ref.size();
	for (auto ite = cindelVec.begin();ite!=cindelVec.end();++ite)
	{
		string TYPE = (*ite)->getType();
		if (TYPE == "REF")inREF(*ite);
		if (TYPE == "5AND3")in5AND3(*ite, MAX,cindelVec );
		if (TYPE == "5OR3")in5OR3(*ite, MAX, cindelVec);
		if (TYPE == "Refand5OR3")inRefand5OR3(*ite, MAX, cindelVec);
	}
	cout << "构造insertion信息成功" << endl;
}
void inREF(IN OUT cindel* &p){
	randomgen randomGen;
	unsigned int pos1, pos2;
		//保证片段是在deletion上取得的
		pos1 = randomGen.randomdeletionsize(p->getdelpos(), p->getdelpos() + p->getdelsize());
		pos2 = randomGen.randomdeletionsize(p->getdelpos(), p->getdelpos() + p->getdelsize());
		if (pos1 > pos2)
		{
			auto temp = pos2;
			pos2 = pos1;
			pos1 = temp;
		}
		string::size_type size = pos2 - pos1;
		bool invert = randomGen.isInversion(cindel::invertrate);
		p->setinvertioninfo(pos1, size, invert, 0, 0, 0);
	}
void in5AND3(IN OUT cindel* &p, IN unsigned max, IN vector<cindel*>& cindelVec){
	randomgen randomGen;
	auto insertionsize = (*p).getdelsize() + randomGen.randominsertiondeviation();
	unsigned size1 = 0, size2 = 0;
	while (true){
		size1 = randomGen.randompos(insertionsize);
		size2 = insertionsize - size1;
		unsigned min = 0;
		if (size1 >= size2)min = size2;
		else min = size1;
		if (min>30&&(double)min/insertionsize>0.1)break;
	}
	bool b; unsigned pos1, pos2; bool invert1, invert2;
	for (size_t i = 0; i < 100000000; ++i)
	{
		b = true;
		pos1 = randomGen.randompos(max);
		pos2 = randomGen.randompos(max);
		for each (auto info in cindelVec)
		{
			if (!(info->isdeletionRange(pos1, size1, max) || info->isdeletionRange(pos2, size2, max))){
				b = false;
				break;
			}
		}
		if (b)break;
	}
	if (b)
	{
		 invert1 = randomGen.isInversion(cindel::invertrate);
		 invert2 = randomGen.isInversion(cindel::invertrate);
		p->setinvertioninfo(pos1, size1, invert1, pos2, size2, invert2);
	}
	else
	{
		cerr << "找不到合适的插入位点请增大ref或者减少indel的个数";
		return;
	}
}
void in5OR3(IN OUT cindel* &p, IN unsigned max, IN vector<cindel*>& cindelVec){
	randomgen randomGen;
	auto insertionsize = (*p).getdelsize() + randomGen.randominsertiondeviation();
	bool b; unsigned pos;
	for (size_t i = 0; i < 100000000; ++i)
	{
		b = true;
		pos = randomGen.randompos(max);
		for each (auto info in cindelVec)
		{
			if (!info->isdeletionRange(pos, insertionsize, max)){
				b = false;
				break;
			}
		}
		if (b)break;
	}
	if (b)
	{
		bool invert = randomGen.isInversion(cindel::invertrate);
		p->setinvertioninfo(pos, insertionsize, invert,0,0,0);
	}
	else
	{
		cerr << "找不到合适的插入位点请增大ref或者减少indel的个数";
		return;
	}
}
void inRefand5OR3(IN OUT cindel* &p, IN unsigned max, IN vector<cindel*>& cindelVec){
	randomgen randomGen;
	auto insertionsize = (*p).getdelsize() + randomGen.randominsertiondeviation();
	unsigned size1 = 0, size2 = 0;
	while (true){
		size1 = randomGen.randompos(insertionsize);
		size2 = insertionsize - size1;
		unsigned min = 0;
		if (size1 >= size2)min = size2;
		else min = size1;
		if (min > 30 && (double)min / insertionsize > 0.1)break;
	}
	bool b; unsigned pos1, pos2; bool invert1, invert2=0;
	pos1 = randomGen.randomdeletionsize(p->getdelpos(), p->getdelpos() +p->getdelsize()- size1);
	invert1 = randomGen.isInversion(cindel::invertrate);
	for (size_t i = 0; i < 100000000; ++i)
	{
		b = true;
		pos2 = randomGen.randompos(max);
		for each (auto info in cindelVec)
		{
			if (!info->isdeletionRange(pos2, insertionsize, max)){
				b = false;
				break;
			}
		}
		if (b)break;
	}
	if (b)
	{
		bool invert2 = randomGen.isInversion(cindel::invertrate);
		p->setinvertioninfo(pos1, size1, invert1, pos2, size2, invert2);
	}
	else
	{
		cerr << "找不到合适的插入位点请增大ref或者减少indel的个数";
		return;
	}
}
void addseqinfo(IN OUT vector<cindel*>& cindelVec, IN  string& ref){
	for (auto ite = cindelVec.begin(); ite != cindelVec.end();++ite)
		(*ite)->setseq(ref);
	cout << "添加序列信息成功" << endl;
}
void generateref(IN vector<cindel*>& cindelVec, IN  OUT string& ref){
	int offset = 0;
	for (auto ite = cindelVec.begin(); ite != cindelVec.end(); ++ite)
	{
		(*ite)->generatorcindel(ref, offset);
		offset += (*ite)->getoffset();
	}

}
void writeinfo(IN vector<cindel*>& cindelVec, IN  string& ref, IN string refsimfile , string cindelinfofile){
	ofstream fout1, fout2;
	fout1.open(refsimfile);
	fout1 << ">chr_19" << endl;

	fout1 << ref;
	fout2.open(cindelinfofile);
	for each (auto var in cindelVec)
		var->writefile(fout2);
	cout << "文件输出成功!!!!!!!!!" << endl;
	for each (auto var in cindelVec)
		delete var;
}
void reversecom(string&  str){
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
void process(string configname ){
	vector<configinfo> convec;
	vector<cindel*> cindelvec;

	string refname, ref;
	readconfig(convec, refname, configname);
	readreference(refname, ref);
	constructcindel(convec, ref, cindelvec);
	constructinsertlinfo(cindelvec, ref);
	addseqinfo(cindelvec, ref);
	generateref(cindelvec, ref);
	writeinfo(cindelvec, ref);
	system("pause");
}



