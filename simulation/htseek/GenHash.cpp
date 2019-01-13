#include"GenHash.h"
bool finalpair::operator<(finalpair& right){
	
	if (this->distance> right.distance)return true;
	else if (this->distance < right.distance)return false;
	else return false;
}
GenHash::GenHash(string& ref, unsigned Kmer) :bucketArray(pow(4, Kmer)),kmerLenth(Kmer){
	MAXSIZE = pow(4, Kmer);
	for (unsigned long long i = 0; i <= ref.size() - Kmer; ++i){
		auto pos = hashfunction(ref.substr(i, Kmer));
		bucketArray[pos].bucketVct.push_back(i);
	}
	Ref = ref;
        init_cache(md5table);
}
GenHash::~GenHash(){
	delete[] md5table;
}
void GenHash::splitseed(string query, vector<string>& seedVct){
	while (query.size()>=kmerLenth)
	{
		seedVct.push_back(query.substr(0,kmerLenth));
		query.erase(0, kmerLenth);
	}
}
void GenHash::reversecom(string& str){
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
unsigned GenHash::hashfunction(string seg){
	{
		unsigned int seed = 131; // the magic number, 31, 131, 1313, 13131
		unsigned int hash = 0;
		unsigned char *p = (unsigned char *)seg.c_str();
		while (*p)
			hash = hash*seed+(*p++);
		return hash % MAXSIZE;
	}
}
pair<bool,BPinfo> GenHash::mappedstring(string query){
         if(query.size()==0){
            pair<bool, BPinfo> re;
	    re.first = false;
	    re.second = BPinfo();
	    return  re;
         }
	vector<finalpair> finalpairvct;
	
	unsigned delLenth = query.size() % kmerLenth;
	query.erase(query.size() - delLenth, delLenth);
	
	unsigned e = query.size() / 100 * SNP;
	unsigned kmere = kmerLenth / KmerSNP;
	if (query.size() % 100 != 0)++e;
	
	vector<string> SeedVct;
	splitseed(query, SeedVct);
	
	vector<bucket> LocVct;
	getLocVct(SeedVct, LocVct);
	  if(LocVct.empty()){
            pair<bool, BPinfo> re;
	    re.first = false;
	    re.second = BPinfo();
	    return  re;
         }
	auto perf = LocVct;
      
	selectPerfectLoction(perf, e);
	
	list<MF> MFVct;
	growthMF(perf, LocVct, MFVct, SeedVct, e, kmere);
	
	auto minsim1 = MFVct.end();
	for (auto i = MFVct.begin(); i != MFVct.end(); ++i)
	{
		if (i->start&&i->end)
		{
			if (minsim1 == MFVct.end()){ minsim1 = i; }
			if (minsim1->SimNum > i->SimNum)minsim1 = i;
		}
	}
	if (minsim1 != MFVct.end())
		{
			BPinfo temp;
			temp.SEG1startpos = minsim1->locVct[0].pos;
			temp.SEG1endpos =minsim1->locVct.back().pos + kmerLenth+delLenth;
			temp.str1 = query;
			finalpair fp;
			fp.bpinfo = temp;
			fp.simnum = minsim1->SimNum;
			long mid=temp.RightBP - temp.LeftBP;
			fp.distance = abs(mid-(long)temp.SEG1startpos);
			finalpairvct.push_back(fp);
		}
	vector<MF> MFstart, MFend;
	for (auto i = MFVct.begin(); i != MFVct.end(); ++i)
	if (i->start)MFstart.push_back(*i);
	else MFend.push_back(*i);
	for (unsigned i = 0; i < MFstart.size(); ++i)
	{
		for (unsigned j = 0; j < MFend.size(); ++j)
		{
			if (MFstart[i].EndSeed + 1 == MFend[j].startSeed)
			{
				if (MFstart[i].locVct.back().isperfectmap&&MFend[j].locVct[0].isperfectmap)
				{
					
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG1startpos = MFstart[i].locVct[0].pos;
					temp.SEG1endpos = MFstart[i].locVct.back().pos + kmerLenth - 1;
					temp.SEG2startpos = MFend[j].locVct[0].pos;
					temp.SEG2endpos = MFend[j].locVct.back().pos + kmerLenth - 1;
					temp.LeftBP = temp.SEG1startpos;
					temp.RightBP = temp.SEG2endpos;
					temp.MidBP = temp.SEG2startpos;
					temp.str1 = query.substr(0, temp.SEG1endpos - temp.SEG1startpos + 1);
					query.erase(0, temp.SEG1endpos - temp.SEG1startpos + 1);
					temp.str2 = query;
					temp.SEG2endpos += delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = MFstart[i].SimNum + MFend[j].SimNum;
					long mid = temp.LeftBP+(temp.RightBP - temp.LeftBP)/2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
				if (MFstart[i].locVct.back().isperfectmap || MFend[j].locVct[0].isperfectmap){
					
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG1startpos = MFstart[i].locVct[0].pos;
					temp.SEG1endpos = MFstart[i].locVct.back().pos + kmerLenth - 1;
					temp.SEG2startpos = MFend[j].locVct[0].pos;
					temp.SEG2endpos = MFend[j].locVct.back().pos + kmerLenth - 1;
					temp.LeftBP = temp.SEG1startpos;
					temp.RightBP = temp.SEG2endpos;
					temp.MidBP = temp.SEG2startpos;
					if (MFstart[i].locVct.back().isperfectmap)
					{
						string kmerread = SeedVct[MFstart[i].EndSeed-1];
						string refread = Ref.substr(temp.SEG1endpos+1, kmerLenth);
						auto count=simplefindBP(temp, 1, kmerread, refread, temp.MidBP);
						temp.MidBP += count;
						temp.SEG1endpos += count;
						temp.SEG2startpos += count;
						temp.str1 = query.substr(0, temp.SEG1endpos - temp.SEG1startpos + 1);
						query.erase(0, temp.SEG1endpos - temp.SEG1startpos + 1);
						temp.str2 = query;

					}
					else
					{
						string kmerread = SeedVct[MFend[j].EndSeed-1];
						string refread = Ref.substr(temp.SEG2startpos-kmerLenth, kmerLenth);
						auto i=simplefindBP(temp, 0, kmerread, refread, temp.MidBP);
						temp.MidBP -= i;
						temp.SEG1endpos -= i;
						temp.SEG2startpos -= i;
						temp.str1 = query.substr(0, temp.SEG1endpos - temp.SEG1startpos + 1);
						query.erase(0, temp.SEG1endpos - temp.SEG1startpos + 1);
						temp.str2 = query;

					}
					temp.SEG2endpos += delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = MFstart[i].SimNum + MFend[j].SimNum;
					long mid = temp.LeftBP+(temp.RightBP - temp.LeftBP)/2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
				else
				{
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG1startpos = MFstart[i].locVct[0].pos;
					temp.SEG1endpos = MFstart[i].locVct.back().pos + kmerLenth - 1 - kmere;
					temp.SEG2startpos = MFend[j].locVct[0].pos - kmere;
					temp.SEG2endpos = MFend[j].locVct.back().pos + kmerLenth - 1;
					temp.LeftBP = temp.SEG1startpos;
					temp.RightBP = temp.SEG2endpos;
					temp.MidBP = temp.SEG2startpos;
					string kmerread = SeedVct[MFstart[i].EndSeed-1].substr(kmerLenth - kmere, kmere) + SeedVct[MFend[j].startSeed-1].substr(0, kmere);
					string left = Ref.substr(temp.SEG1endpos + 1, 2 * kmere);
					string right = Ref.substr(temp.SEG2startpos, 2 * kmere);
					auto count=shortfindBP(temp, kmerread, left, right,kmere);
					if (count == -1)continue; 
					temp.SEG1endpos += count;
					temp.SEG2startpos += count;
					temp.str1 = query.substr(0, temp.SEG1endpos - temp.SEG1startpos + 1);
					query.erase(0, temp.SEG1endpos - temp.SEG1startpos + 1);
					temp.str2 = query;
					temp.SEG2endpos += delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = MFstart[i].SimNum + MFend[j].SimNum;
					long mid = temp.LeftBP+(temp.RightBP - temp.LeftBP)/2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
			}
			if (MFstart[i].EndSeed + 1 == MFend[j].startSeed - 1)
			{
				BPinfo temp;
				temp.isCombination = 1;
				temp.SEG1startpos = MFstart[i].locVct[0].pos;
				temp.SEG1endpos = MFstart[i].locVct.back().pos + kmerLenth - 1;
				temp.SEG2startpos = MFend[j].locVct[0].pos - kmerLenth;
				temp.SEG2endpos = MFend[j].locVct.back().pos + kmerLenth - 1;
				temp.LeftBP = temp.SEG1startpos;
				temp.RightBP = temp.SEG2endpos;
				temp.MidBP = temp.SEG2startpos;
				string left = Ref.substr(temp.SEG1endpos + 1, kmerLenth);
				string right = Ref.substr(temp.SEG2startpos, kmerLenth);
				auto count=longfindBP(temp, SeedVct[MFstart[i].EndSeed-1], left, right,kmere);
				if (count == -1)continue;
				temp.SEG1endpos += count;
				temp.SEG2startpos +=count;
				temp.str1 = query.substr(0, temp.SEG1endpos - temp.SEG1startpos + 1);
				query.erase(0, temp.SEG1endpos - temp.SEG1startpos + 1);
				temp.str2 = query;
				temp.SEG2endpos += delLenth;
				finalpair fp;
				fp.bpinfo = temp;
				fp.simnum = MFstart[i].SimNum + MFend[j].SimNum;
				long mid = temp.LeftBP+(temp.RightBP - temp.LeftBP)/2;
				fp.distance = abs(mid - (long)temp.SEG1startpos);
				finalpairvct.push_back(fp);
				continue;
			}
		}
	}
	
	string rquery = query;
	reversecom(rquery);
	vector<string> rSeedVct;
	splitseed(rquery, rSeedVct);
	vector<bucket> rLocVct;
	getLocVct(rSeedVct, rLocVct);
	auto rperf = rLocVct;
         if(rLocVct.empty()){
            pair<bool, BPinfo> re;
	    re.first = false;
	    re.second = BPinfo();
	    return  re;
         }
	selectPerfectLoction(rperf, e);
	list<MF> rMFVct;
	growthMF(rperf, rLocVct, rMFVct, rSeedVct, e, kmere);
	/*判断是不是5R或者3R,或者倒置
	*
	*
	*
	*
	*/
	auto minsim = rMFVct.end();
	for (auto i = rMFVct.begin(); i != rMFVct.end(); ++i)
	{
		if (i->start&&i->end)
		{
			if (minsim == rMFVct.end()){ minsim = i; }
			if (minsim->SimNum > i->SimNum)minsim = i;
		}
	}
		if (minsim != rMFVct.end())
		{//如果头尾都到了边界,那么我们就认为这是一个最长的ref
			BPinfo temp;
			temp.SEG1startpos = minsim->locVct[0].pos;
			temp.SEG1endpos = minsim->locVct.back().pos + kmerLenth;
			temp.str1 = rquery;
			temp.SEG1isR = 1;
			temp.SEG1startpos -= delLenth;
			finalpair fp;
			fp.bpinfo = temp;
			fp.simnum = minsim->SimNum;
			long mid = temp.LeftBP+(temp.RightBP - temp.LeftBP)/2;
			fp.distance = abs(mid - (long)temp.SEG1startpos);
			finalpairvct.push_back(fp);
		
		}
	/********************************************
	 ********************************************
	 ** 判断是不是反向互补的两个片段拼接组合而成的**/
	vector<MF> rMFstart, rMFend;
	for (auto i = rMFVct.begin(); i != rMFVct.end(); ++i)
	if (i->start)rMFstart.push_back(*i);
	else rMFend.push_back(*i);
	for (unsigned i = 0; i < rMFstart.size(); ++i)
	{
		for (unsigned j = 0; j < rMFend.size(); ++j)
		{
			if (rMFstart[i].EndSeed + 1 == rMFend[j].startSeed)
			{
				if (rMFstart[i].locVct.back().isperfectmap&&rMFend[j].locVct[0].isperfectmap)
				{
					//说明恰好kmer的分界就是BP
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG2startpos = rMFstart[i].locVct[0].pos;
					temp.SEG2endpos = rMFstart[i].locVct.back().pos + kmerLenth - 1;
					temp.SEG1startpos = rMFend[j].locVct[0].pos;
					temp.SEG1endpos = rMFend[j].locVct.back().pos + kmerLenth - 1;
					temp.LeftBP = temp.SEG2startpos;
					temp.RightBP = temp.SEG1endpos;
					temp.MidBP = temp.SEG1startpos;
					temp.str2 = rquery.substr(0, temp.SEG2endpos - temp.SEG2startpos + 1);
					rquery.erase(0, temp.SEG2endpos - temp.SEG2startpos + 1);
					temp.str1 = rquery;
					temp.SEG2isR = 1;
					temp.SEG1isR = 1;
					temp.SEG2startpos -= delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = rMFstart[i].SimNum + rMFend[j].SimNum;
					long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
				if (rMFstart[i].locVct.back().isperfectmap || rMFend[j].locVct[0].isperfectmap){
					//若其中一个完美比对另一个相似比对
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG2startpos = rMFstart[i].locVct[0].pos;
					temp.SEG2endpos = rMFstart[i].locVct.back().pos + kmerLenth - 1;
					temp.SEG1startpos = rMFend[j].locVct[0].pos;
					temp.SEG1endpos = rMFend[j].locVct.back().pos + kmerLenth - 1;
					temp.LeftBP = temp.SEG2startpos;
					temp.RightBP = temp.SEG1endpos;
					temp.MidBP = temp.SEG1startpos;
					if (rMFstart[i].locVct.back().isperfectmap)
					{//左边完美右边相似
						string kmerread = rSeedVct[rMFstart[i].EndSeed-1];
						string refread = Ref.substr(temp.SEG2endpos + 1, kmerLenth);
						auto count=simplefindBP(temp, 1, kmerread, refread, temp.MidBP);
						temp.MidBP += count;
						temp.SEG2endpos += count;
						temp.SEG1startpos += count;
						temp.str2 = rquery.substr(0, temp.SEG2endpos - temp.SEG2startpos + 1);
						rquery.erase(0, temp.SEG2endpos - temp.SEG2startpos + 1);
						temp.str1 = rquery;


					}
					else
					{//左边相似右边完美
						string kmerread = rSeedVct[rMFstart[i].EndSeed-1];
						string refread = Ref.substr(temp.SEG1startpos-kmerLenth, kmerLenth);
						auto i=simplefindBP(temp, 0, kmerread, refread, temp.MidBP);
						temp.MidBP -= i;
						temp.SEG2endpos -= i;
						temp.SEG1startpos -= i;
						temp.str2 = rquery.substr(0, temp.SEG2endpos - temp.SEG2startpos + 1);
						rquery.erase(0, temp.SEG2endpos - temp.SEG2startpos + 1);
						temp.str1 = rquery;

					}
					temp.SEG1isR = 1;
					temp.SEG2isR = 1;
					temp.SEG2startpos -= delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = rMFstart[i].SimNum + rMFend[j].SimNum;
					long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
				else
				{
					//两个都是sim比对
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG2startpos = rMFstart[i].locVct[0].pos;
					temp.SEG2endpos = rMFstart[i].locVct.back().pos + kmerLenth - 1 - kmere;
					temp.SEG1startpos = rMFend[j].locVct[0].pos - kmere;
					temp.SEG1endpos = rMFend[j].locVct.back().pos + kmerLenth - 1;
					temp.LeftBP = temp.SEG2startpos;
					temp.RightBP = temp.SEG1endpos;
					temp.MidBP = temp.SEG1startpos;
					string kmerread = rSeedVct[rMFstart[i].EndSeed - 1].substr(kmerLenth - kmere, kmere) + SeedVct[rMFend[j].startSeed - 1].substr(0, kmere);
					string left = Ref.substr(temp.SEG2endpos + 1, 2 * kmere);
					string right = Ref.substr(temp.SEG1startpos, 2 * kmere);
					auto count=shortfindBP(temp, kmerread, left, right,kmere);
					if (count == -1)continue; //比对失败,变异太高
					temp.SEG2endpos += count;
					temp.SEG1startpos += count;
					temp.str1 = rquery.substr(0, temp.SEG2endpos - temp.SEG2startpos + 1);
					rquery.erase(0, temp.SEG2endpos - temp.SEG2startpos + 1);
					temp.str1 = rquery;
					temp.SEG1isR = 1;
					temp.SEG2isR = 1;
					temp.SEG2startpos -= delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = rMFstart[i].SimNum + rMFend[j].SimNum;
					long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
			}
			//下面讨论另一种特殊情况.断点刚好卡在kmer的中间点附近导致kmer完全无法比对
			if (rMFstart[i].EndSeed + 1 == rMFend[j].startSeed - 1)
			{
				BPinfo temp;
				temp.isCombination = 1;
				temp.SEG2startpos = rMFstart[i].locVct[0].pos;
				temp.SEG2endpos = rMFstart[i].locVct.back().pos + kmerLenth - 1;
				temp.SEG1startpos = rMFend[j].locVct[0].pos - kmerLenth;//特殊处理
				temp.SEG1endpos = rMFend[j].locVct.back().pos + kmerLenth - 1;
				temp.LeftBP = temp.SEG2startpos;
				temp.RightBP = temp.SEG1endpos;
				temp.MidBP = temp.SEG2startpos;
				string left = Ref.substr(temp.SEG2endpos + 1, kmerLenth);
				string right = Ref.substr(temp.SEG1startpos, kmerLenth);
				auto count=longfindBP(temp, rSeedVct[rMFstart[i].EndSeed-1], left, right,kmere);
				if (count == -1)continue; //比对失败,变异太高
				temp.SEG2endpos += count;
				temp.SEG1startpos += count;
				temp.str1 = rquery.substr(0, temp.SEG2endpos - temp.SEG2startpos + 1);
				rquery.erase(0, temp.SEG2endpos - temp.SEG2startpos + 1);
				temp.str2 = rquery;
				temp.SEG1isR = 1;
				temp.SEG2isR = 1;
				temp.SEG2startpos -= delLenth;
				finalpair fp;
				fp.bpinfo = temp;
				fp.simnum = rMFstart[i].SimNum + rMFend[j].SimNum;
				long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
				fp.distance = abs(mid - (long)temp.SEG1startpos);
				finalpairvct.push_back(fp);
				continue;
			}
		}
	}
	/****************
	接下来比较复杂,就是两个片段一个从正向上打下来,一个从反向上打下来
	左正右反


	**********************/
	for (unsigned i = 0; i < MFstart.size(); ++i)
	{
		for (unsigned j = 0; j < rMFstart.size(); ++j)
		{
			/*if (MFstart[i].locVct.size() + rMFstart[j].locVct.size() == SeedVct.size()+1)
			{
				MFstart[i].locVct.pop_back();
				rMFstart[j].locVct.pop_back();
			}*/
			if (MFstart[i].locVct.size() + rMFstart[j].locVct.size() == SeedVct.size())
			{
				if (MFstart[i].locVct.back().isperfectmap&&rMFstart[j].locVct.back().isperfectmap)
				{
					//说明恰好kmer的分界就是BP
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG1startpos = MFstart[i].locVct[0].pos;
					temp.SEG1endpos = MFstart[i].locVct.back().pos + kmerLenth - 1;
					temp.SEG2startpos = rMFstart[j].locVct[0].pos;
					temp.SEG2endpos = rMFstart[j].locVct.back().pos + kmerLenth - 1;
					temp.LeftBP = temp.SEG1startpos;
					temp.RightBP = temp.SEG2endpos;
					temp.MidBP = temp.SEG2startpos;
					temp.str1 = Ref.substr(temp.SEG1startpos, temp.SEG1endpos - temp.SEG1startpos + 1);
					temp.str2 = Ref.substr(temp.SEG2startpos, temp.SEG2endpos - temp.SEG2startpos + 1);
					temp.SEG1isR = 0;
					temp.SEG2isR = 1;
					temp.SEG2startpos -= delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = MFstart[i].SimNum + rMFstart[j].SimNum;
					long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
				if (MFstart[i].locVct.back().isperfectmap || rMFstart[j].locVct.back().isperfectmap)
				{
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG1startpos = MFstart[i].locVct[0].pos;
					temp.SEG1endpos = MFstart[i].locVct.back().pos + kmerLenth - 1;
					temp.SEG2startpos = rMFstart[j].locVct[0].pos;
					temp.SEG2endpos = rMFstart[j].locVct.back().pos + kmerLenth - 1;
					temp.LeftBP = temp.SEG1startpos;
					temp.RightBP = temp.SEG2endpos;
					temp.MidBP = temp.SEG2startpos;
					temp.SEG1isR = 0;
					temp.SEG2isR = 1;
					if (MFstart[i].locVct.back().isperfectmap){
						//如果左边片段是完美比对的 
						string kmerread = SeedVct[MFstart[i].EndSeed-1];
						string refread = Ref.substr(temp.SEG1endpos+1, kmerLenth);
						auto count=simplefindBP(temp, 1, kmerread, refread, temp.MidBP);
						temp.SEG1endpos += count;
						temp.SEG2endpos -= count;
						temp.str1 = Ref.substr(temp.SEG1startpos, temp.SEG1endpos - temp.SEG1startpos + 1);
						temp.str2 = Ref.substr(temp.SEG1startpos, temp.SEG2endpos - temp.SEG1startpos + 1);
					}
					else
					{//右边片段是完美比对的
						cout<<rSeedVct.size()<<"     "<<rMFstart[j].EndSeed<<"&"<<endl;
						string kmerread = rSeedVct[rMFstart[j].EndSeed-1];
						string refread = Ref.substr(temp.SEG2endpos + 1, kmerLenth);
						auto count = simplefindBP(temp, 1, kmerread, refread, temp.MidBP);
						temp.SEG1endpos -= count;
						temp.SEG2endpos += count;
						temp.str1 = Ref.substr(temp.SEG1startpos, temp.SEG1endpos - temp.SEG1startpos + 1);
						temp.str2 = Ref.substr(temp.SEG2startpos, temp.SEG2endpos - temp.SEG2startpos + 1);
					}
					temp.SEG2startpos -= delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = MFstart[i].SimNum + rMFstart[j].SimNum;
					long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
				else
				{
					//两个都是sim的情况
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG1startpos = MFstart[i].locVct[0].pos;
					temp.SEG1endpos = MFstart[i].locVct.back().pos + kmerLenth - 1 - kmere;
					temp.SEG2startpos = rMFstart[j].locVct[0].pos;
					temp.SEG2endpos = rMFstart[j].locVct.back().pos + kmerLenth - 1+kmere;
					temp.LeftBP = temp.SEG1startpos;
					temp.RightBP = temp.SEG2endpos;
					temp.MidBP = temp.SEG2startpos;
					string kmerread = SeedVct[MFstart[i].EndSeed - 1].substr(kmerLenth - kmere, kmere);
					kmerread += SeedVct[MFstart[i].EndSeed-1].substr(0, kmere);
					string left = Ref.substr(temp.SEG1endpos + 1, 2 * kmere);
					string right = Ref.substr(temp.SEG1endpos - 2 * kmere + 1, 2 * kmere);
					reversecom(right);
					auto count=shortfindBP(temp, kmerread, left, right,kmere);
					if (count == -1)continue; //比对失败,变异太高
					//是不是有问题??????????????????????????
					temp.SEG1endpos += count;
					temp.SEG2endpos -=count;
					temp.str1 = Ref.substr(temp.SEG1startpos, temp.SEG1endpos - temp.SEG1startpos + 1);
					temp.str2 = Ref.substr(temp.SEG2startpos, temp.SEG2endpos - temp.SEG2startpos + 1);
					temp.SEG1isR = 0;
					temp.SEG2isR = 1;
					temp.SEG2startpos -= delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = MFstart[i].SimNum + rMFstart[j].SimNum;
					long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
				
			}
			if (MFstart[i].EndSeed + rMFstart[j].EndSeed == SeedVct.size()-1)
			{
				//如果是最复杂的情况,也就是两个片段左边正向,右边反向,并且之间还空余着一个kmer
				BPinfo temp;
				temp.isCombination = 1;
				temp.SEG1startpos = MFstart[i].locVct[0].pos;
				temp.SEG1endpos = MFstart[i].locVct.back().pos + kmerLenth - 1;
				temp.SEG2startpos = rMFstart[j].locVct[0].pos;
				temp.SEG2endpos = rMFstart[j].locVct.back().pos + 2*kmerLenth - 1;//特殊处理
				temp.LeftBP = temp.SEG1startpos;
				temp.RightBP = temp.SEG2endpos;
				temp.MidBP = temp.SEG2startpos;
				string left = Ref.substr(temp.SEG1endpos + 1, kmerLenth);
				string right = Ref.substr(temp.SEG2endpos-kmere+1, kmerLenth);
				reversecom(right);
				auto count=longfindBP(temp, SeedVct[MFstart[i].EndSeed-1], left, right,kmere);
				if (count == -1)continue; //比对失败,变异太高
				temp.SEG1endpos += count;
				temp.SEG2endpos -= kmerLenth-count;
				temp.str1 = Ref.substr(temp.SEG1startpos, temp.SEG1endpos - temp.SEG1startpos + 1);
				temp.str2 = Ref.substr(temp.SEG2startpos, temp.SEG2endpos - temp.SEG2startpos + 1);
				temp.SEG1isR = 0;
				temp.SEG2isR = 1;
				temp.SEG2startpos -= delLenth;
				finalpair fp;
				fp.bpinfo = temp;
				fp.simnum = MFstart[i].SimNum + rMFstart[j].SimNum;
				long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
				fp.distance = abs(mid - (long)temp.SEG1startpos);
				finalpairvct.push_back(fp);
				continue;
			}
		}
	}
	/*
	最后一种情况
	左反向右正向
	*/
	for (unsigned i = 0; i < rMFend.size(); ++i)
	{
		for (unsigned j = 0; j < MFend.size(); ++j)
		{
			/*if (rMFend[i].locVct.size() + MFend[j].locVct.size() == SeedVct.size() + 1)
			{
				MFstart[i].locVct.erase(MFstart[i].locVct.begin());
				rMFstart[j].locVct.erase(rMFstart[j].locVct.begin());
			}*/
			if (rMFend[i].locVct.size() + MFend[j].locVct.size() == SeedVct.size())
			{
				if (rMFend[i].locVct[0].isperfectmap&& MFend[j].locVct[0].isperfectmap)
				{
					//说明恰好kmer的分界就是BP
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG1startpos = rMFend[i].locVct[0].pos;
					temp.SEG1endpos = rMFend[i].locVct.back().pos + kmerLenth - 1;
					temp.SEG2startpos = MFend[j].locVct[0].pos;
					temp.SEG2endpos = MFend[j].locVct.back().pos + kmerLenth - 1;
					temp.LeftBP = temp.SEG1startpos;
					temp.RightBP = temp.SEG2endpos;
					temp.MidBP = temp.SEG2startpos;
					temp.str1 = Ref.substr(temp.SEG1startpos, temp.SEG1endpos - temp.SEG1startpos + 1);
					temp.str2 = Ref.substr(temp.SEG2startpos, temp.SEG2endpos - temp.SEG2startpos + 1);
					temp.SEG1isR = 1;
					temp.SEG2isR = 0;
					temp.SEG2endpos += delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = rMFend[i].SimNum + MFend[j].SimNum;
					long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
				if (rMFend[i].locVct[0].isperfectmap || MFend[j].locVct[0].isperfectmap)
				{
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG1startpos = rMFend[i].locVct[0].pos;
					temp.SEG1endpos = rMFend[i].locVct.back().pos + kmerLenth - 1;
					temp.SEG2startpos = MFend[j].locVct[0].pos;
					temp.SEG2endpos = MFend[j].locVct.back().pos + kmerLenth - 1;
					temp.LeftBP = temp.SEG1startpos;
					temp.RightBP = temp.SEG2endpos;
					temp.MidBP = temp.SEG2startpos;
					temp.SEG1isR = 1;
					temp.SEG2isR = 0;
					if (rMFend[i].locVct[0].isperfectmap){
						//如果左边片段是完美比对的 
						string kmerread = SeedVct[rMFend[i].startSeed-2];
						string refread = Ref.substr(temp.SEG1startpos-kmerLenth, kmerLenth);
						auto count = simplefindBP(temp, 0, kmerread, refread, temp.MidBP);//但是实际上我们做的是向左扩展
						temp.SEG1startpos -= count;
						temp.SEG2startpos += count;
						temp.str1 = Ref.substr(temp.SEG1startpos, temp.SEG1endpos - temp.SEG1startpos + 1);
						temp.str2 = Ref.substr(temp.SEG1startpos, temp.SEG2endpos - temp.SEG1startpos + 1);
					}
					else
					{//右边片段是完美比对的
						string kmerread = SeedVct[MFend[j].startSeed-2];
						string refread = Ref.substr(temp.SEG2startpos -kmerLenth , kmerLenth);
						auto count = simplefindBP(temp, 0, kmerread, refread, temp.MidBP);
						temp.SEG1startpos += count;
						temp.SEG2startpos -= count;
						temp.str1 = Ref.substr(temp.SEG1startpos, temp.SEG1endpos - temp.SEG1startpos + 1);
						temp.str2 = Ref.substr(temp.SEG2startpos, temp.SEG2endpos - temp.SEG2startpos + 1);
					}
					temp.SEG2endpos += delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = rMFend[i].SimNum + MFend[j].SimNum;
					long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}
				else
				{
					//两个都是sim的情况
					BPinfo temp;
					temp.isCombination = 1;
					temp.SEG1startpos = rMFend[i].locVct[0].pos+kmere;
					temp.SEG1endpos = rMFend[i].locVct.back().pos + kmerLenth-1;
					temp.SEG2startpos = MFend[j].locVct[0].pos-kmere;
					temp.SEG2endpos = MFend[j].locVct.back().pos + kmerLenth - 1 ;
					temp.LeftBP = temp.SEG1startpos;
					temp.RightBP = temp.SEG2endpos;
					temp.MidBP = temp.SEG2startpos;
					string kmerread = SeedVct[MFend[j].startSeed - 1].substr(0, kmere);
					kmerread = SeedVct[MFend[j].startSeed].substr(kmerLenth-kmere, kmere)+kmerread;
					string left = Ref.substr(temp.SEG1startpos - 2 * kmere, 2 * kmere);
					string right = Ref.substr(temp.SEG2startpos - 2 * kmere, 2 * kmere);
					reversecom(left);
					auto count = shortfindBP(temp, kmerread, left, right, kmere);
					if (count == -1)continue; //比对失败,变异太高
					temp.SEG1startpos -= count;
					temp.SEG2startpos += count;
					temp.str1 = Ref.substr(temp.SEG1startpos, temp.SEG1endpos - temp.SEG1startpos + 1);
					temp.str2 = Ref.substr(temp.SEG2startpos, temp.SEG2endpos - temp.SEG2startpos + 1);
					temp.SEG1isR = 1;
					temp.SEG2isR = 0;
					temp.SEG2endpos += delLenth;
					finalpair fp;
					fp.bpinfo = temp;
					fp.simnum = rMFend[i].SimNum + MFend[j].SimNum;
					long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
					fp.distance = abs(mid - (long)temp.SEG1startpos);
					finalpairvct.push_back(fp);
					continue;
				}

			}
			//周四要完成的!!!!
			if (rMFend[i].locVct.size() + MFend[j].locVct.size() == SeedVct.size() - 1)
			{
				//如果是最复杂的情况,也就是两个片段左边正向,右边反向,并且之间还空余着一个kmer
				BPinfo temp;
				temp.isCombination = 1;
				temp.SEG1startpos = rMFend[i].locVct[0].pos;
				temp.SEG1endpos = rMFend[i].locVct.back().pos + kmerLenth - 1;
				temp.SEG2startpos = MFend[j].locVct[0].pos-kmerLenth;//特殊处理
				temp.SEG2endpos = MFend[j].locVct.back().pos +kmerLenth - 1;
				temp.LeftBP = temp.SEG1startpos;
				temp.RightBP = temp.SEG2endpos;
				temp.MidBP = temp.SEG2startpos;
				string left = Ref.substr(temp.SEG1startpos-kmerLenth,kmerLenth);
				reversecom(left);
				string right = Ref.substr(temp.SEG2startpos - kmerLenth, kmerLenth);
				auto count = longfindBP(temp, SeedVct[MFend[j].startSeed-2], left, right, kmere);
				if (count == -1)continue; //比对失败,变异太高
				//这里是有问题的
				temp.SEG1startpos -= count;
				temp.SEG2startpos +=kmerLenth-count;
				temp.str1 = Ref.substr(temp.SEG1startpos, temp.SEG1endpos - temp.SEG1startpos + 1);
				temp.str2 = Ref.substr(temp.SEG2startpos, temp.SEG2endpos - temp.SEG2startpos + 1);
				temp.SEG1isR = 1;
				temp.SEG2isR = 0;
				temp.SEG2endpos += delLenth;
				finalpair fp;
				fp.bpinfo = temp;
				fp.simnum = rMFend[i].SimNum + MFend[j].SimNum;
				long mid = temp.LeftBP + (temp.RightBP - temp.LeftBP) / 2;
				fp.distance = abs(mid - (long)temp.SEG1startpos);
				finalpairvct.push_back(fp);
				continue;
			}
		}
	}
	if (finalpairvct.empty()){
		pair<bool, BPinfo> re;
		re.first = false;
		re.second = BPinfo();
		return  re;
	}
	pair<bool, BPinfo> re;
	re.first = true;
	re.second = max(finalpairvct.begin(), finalpairvct.begin())->bpinfo;
	return re;
}
void GenHash::getLocVct(vector<string>& seedVct, vector<bucket>& LocVct){
	for (unsigned i = 0; i < seedVct.size(); ++i){
		LocVct.push_back(bucketArray[hashfunction(seedVct[i])]);
	}
}
void GenHash::selectPerfectLoction(vector<bucket>& LocVct, unsigned e){
	//if (LocVct.size() < 3)
	//{
	//	std::cerr << "seedÊýÁ¿Ð¡ÓÚ3£¬³ÌÐòÍË³ö" << endl;
	//	std::exit(0);
	//}
	unsigned  prepos, nowpos, nextpos;
        
	for (unsigned i = 1; i < LocVct.size() - 1; ++i)//Íš¹ýÕâžöÑ­»·ÎÒÃÇÄÜ¹»³õ²œ±êŒÇ³öÎÒÃÇµÄÍêÃÀÎ»ÖÃ
	{  
                
		for (unsigned j = 0; j < LocVct[i].bucketVct.size(); ++j)
		{
			nowpos = LocVct[i].bucketVct[j].pos;
			for (unsigned k1 = 0; k1 < LocVct[i - 1].bucketVct.size(); ++k1)
			{
				auto a = LocVct[i].bucketVct[j].pos;
				auto b = LocVct[i - 1].bucketVct[k1].pos;
				if (a < b)
				{
					auto temp = a;
					a = b;
					b = temp;
				}
				if (a+e >= kmerLenth + b&&a - b <= kmerLenth + e)
				{	//Èç¹ûÎªÕæ£¬ËµÃ÷ÕâÁœžöÎ»ÖÃÊÇÓÐÏà»¥ÏìÓŠµÄ
					prepos = LocVct[i - 1].bucketVct[k1].pos;
					for (unsigned k2 = 0; k2 < LocVct[i + 1].bucketVct.size(); ++k2)
					{
						auto a = LocVct[i].bucketVct[j].pos;
						auto b = LocVct[i + 1].bucketVct[k2].pos;
						if (a < b)
						{
							auto temp = a;
							a = b;
							b = temp;
						}
						if (a +e >= kmerLenth + b&&a - b <= kmerLenth + e)
						{
							nextpos = LocVct[i + 1].bucketVct[k2].pos;
							//µœÕâÎªÖ¹£¬ÎÒÃÇÒÑŸ­Ñ¡Ôñ³öÁËºÏÊÊµÄÈýÁ¬Î»ÖÃ£¬ÒòŽËÒª°ÑÖÐŒäµÄÎ»ÖÃ¶šÒå³ÉÎªÒ»žöperfectmaploc
							LocVct[i].bucketVct[j].isperfectmap = 1;
							if (i - 1 == 0)//Èç¹ûi-1ÊÇÒ»žö0£¬ŸÍ±íÊŸÈýÁ¬Î»ÖÃµÄµÚÒ»žöÎ»ÖÃÊÇµÚÒ»²ã£¬Ò²Òªœ«ËüÉèÖÃ³ÉÎªÒ»žöperfectmaploc
							{
								LocVct[i - 1].bucketVct[k1].isperfectmap = 1;
							}
							if (i + 1 == LocVct.size() - 1)//Èç¹ûi+1ÊÇ×îºóÒ»²ãlayer
							{
								LocVct[i + 1].bucketVct[k2].isperfectmap = 1;
							}
						}
					}
				}
			}
		}
	}
	//œÓ×ÅÉŸ³ýµôÆäÖÐËùÓÐ²»ÍêÃÀµÄlocbp
	for (auto i = LocVct.begin(); i != LocVct.end(); ++i)
	{
		for (auto j = i->bucketVct.begin(); j != i->bucketVct.end();)
		{
			if (j->isperfectmap == false)
				j = i->bucketVct.erase(j);
			else
				++j;
		}
	}
}
void GenHash::initMF(vector<bucket>& LocVct, unsigned e, vector<MF>& MFVct){
	//Íš¹ýÍêÃÀLocpos¹¹ÔìÒ»žö³õÊŒµÄMFÐòÁÐ
	for (unsigned i = 0; i < e; ++i){//iÓÃÀŽ¿ØÖÆµ±Ç°ÆðÊŒseedµÄ²ãÊý(²ãÊý²»¿ÉÄÜŽóÓÚe)
		for (auto j = LocVct[i].bucketVct.begin(); j != LocVct[i].bucketVct.end();)//jÓÃÀŽ±íÊŸµ±Ç°µÚÒ»žöseedµÄÑ¡È¡
		{
			//Éú³ÉÒ»žö³õÊŒµÄMF
			MF temp;
			temp.startSeed = i + 1; temp.EndSeed = i + 1; temp.MappNum += 1; temp.locVct.push_back(*j);
			j = LocVct[i].bucketVct.erase(j);
			//³¢ÊÔŽÓÏÂÒ»²ã¿ªÊŒ³õÊŒ»¯À©Õ¹MF
			for (unsigned L = i + 1; L < e; ++L)
			{
				//Éè¶šÒ»žöboolÓÃÀŽÅÐ¶šÕâžöÕâ²ãÓÐÃ»ÓÐŒÓÈëÍêÃÀÎ»ÖÃµœMFÖÐ
				bool tag = 0;
				//ÓÃk±íÊŸÉžÑ¡µ±Ç°²ãÊýµÄseed
				for (auto k = LocVct[L].bucketVct.begin(); k != LocVct[L].bucketVct.end();)
				{
					auto a = temp.locVct.back().pos;
					auto b = k->pos;
					if (a < b)
					{
						auto temp = a;
						a = b;
						b = temp;
					}
					if (a - b >= kmerLenth - e&&a - b <= kmerLenth + e){

						//Èç¹ûÁœžöÍêÃÀÎ»ÖÃÊÇÏà»¥ÏìÓŠµÄ,ÄÇÃŽÎÒÃÇŸÍœ«ÕâžöÎ»ÖÃÒ²ŒÓÈëµ±Ç°µÄMFÖÐ
						++temp.EndSeed; ++temp.MappNum; temp.locVct.push_back(*k);
						k = LocVct[L].bucketVct.erase(k);
						tag = 1;
						break;
					}
					else
						++k;
				}
				if (tag == 0)
				{//Èç¹ûÔÚL²ã²¢Ã»ÓÐÕÒµœºÏÊÊµÄÍêÃÀLocPos ÎÒÃÇŸÍÌø³öÑ­»·
					break;
				}
			}
			MFVct.push_back(temp);
		}
	}
}
void GenHash::growthMF(vector<bucket>& PerLocVct, vector<bucket>& Allvct, list<MF>& MFVct, vector<string>& seedVct, unsigned e, unsigned kmere){
	unsigned MAX = PerLocVct.size();//µÃµœÁË×îŽó²ãÊý
	//ÏÂÃæ¶ÔPerLocVctÖÐµÄÃ¿Ò»žöÎ»ÖÃœøÐÐË«ÏòÀ©Õ¹
	for (int Layer = 0; Layer < PerLocVct.size(); ++Layer)
	{
		
		for (auto i = PerLocVct[Layer].bucketVct.begin(); i != PerLocVct[Layer].bucketVct.end();)
		{
			MF temp; ++temp.MappNum; temp.locVct.push_back(i->pos);
			temp.startSeed = temp.EndSeed = Layer + 1;
			temp.locVct[0].isperfectmap = 1;
			i = PerLocVct[Layer].bucketVct.erase(i);
			/*********************************/
			//Ñ¡ÔñLayer²ãµÄµÚižöµãœøÐÐÇ°ÏòÀ©Õ¹
			for (int j = Layer - 1; j >= 0; --j){
				bool tag = 0;
				//jÊÇ²ãÊýµÄÏÂ±ê,Ê×ÏÈÔÚPerLocÖÐÑ¡ÔñÊÇ·ñÓÐºÏÊÊµÄÎ»ÖÃ
				for (auto k = PerLocVct[j].bucketVct.begin(); k != PerLocVct[j].bucketVct.end(); ++k)
				{
					auto a = temp.locVct[0].pos;
					auto b = k->pos;
					if (a < b)
					{
						continue;
					}
					if (a -b>=kmerLenth-2*kmere&&a - b <= kmerLenth + 2*kmere){
						locbp LOC(k->pos);
						LOC.isperfectmap = 1;
						temp.locVct.insert(temp.locVct.begin(), LOC);
						++temp.MappNum;
						--temp.startSeed;
						tag = 1;
						PerLocVct[j].bucketVct.erase(k);
						break;
					}
				}
				if (tag == 0){
					//Èç¹ûtagÊÇ0ÎÒÃÇÐèÒªœøÐÐALLVctœøÐÐŸ«×ŒÑéÖ€
					for (auto k = 0; k < Allvct[j].bucketVct.size(); ++k){
						auto a = temp.locVct[0].pos;
						auto b = Allvct[j].bucketVct[k].pos;
						if (a < b)
						{
							continue;
						}
						if (a - b >= kmerLenth -2* kmere&&a - b <= kmerLenth + 2*kmere){
							locbp LOC(Allvct[j].bucketVct[k].pos);
							LOC.isperfectmap = 1;
							temp.locVct.insert(temp.locVct.begin(), LOC);
							++temp.MappNum;
							--temp.startSeed;
							tag = 1;
							break;
						}
					}
				}
				if (tag == 0){
					//Èç¹ûtag»¹ÊÇ0,ÄÇÃŽÎÒÃÇŸÍÒªÑ¡ÔñœøÐÐsim±È¶Ô
					auto Pos = temp.locVct[0].pos - kmerLenth;
					string refKmer = Ref.substr(Pos, kmerLenth);
					string seedKmer = seedVct[j];
					if (edistance(seedKmer, refKmer) <= kmere)
					{
						//Èç¹û±àŒ­ŸàÀëÐ¡ÓÚe,ÎÒÃÇŸÍÅÐ¶šÕâžöseedÊÇÒ»žöÏàËÆseed
						locbp LOC(Pos);
						LOC.isperfectmap = 0;
						temp.locVct.insert(temp.locVct.begin(), LOC);
						temp.SimNum++;
						--temp.startSeed;
						tag = 1;
					}
				}
				if (tag == 0)break;//ËµÃ÷Ç°Ïò±È¶ÔÖÕÖ¹
			}
			/***************************************/
			//ÏÂÃæÊÇºóÏòÀ©Õ¹,ºÍÇ°ÏòÀ©Õ¹ËŒÂ·ÀàËÆ
			for (auto j = Layer + 1; j < PerLocVct.size(); ++j){
				bool tag = 0;
				//jÊÇ²ãÊýµÄÏÂ±ê,Ê×ÏÈÔÚPerLocÖÐÑ¡ÔñÊÇ·ñÓÐºÏÊÊµÄÎ»ÖÃ
				for (auto k = PerLocVct[j].bucketVct.begin(); k != PerLocVct[j].bucketVct.end(); ++k)
				{
					auto a = temp.locVct.back().pos;
					auto b = k->pos;
					if (a > b)
					{
						continue;
					}
					if (b-a >= kmerLenth -2* kmere&&b-a <= kmerLenth + 2*kmere){
						locbp LOC(k->pos);
						LOC.isperfectmap = 1;
						temp.locVct.push_back(LOC);
						++temp.MappNum;
						++temp.EndSeed;
						tag = 1;
						PerLocVct[j].bucketVct.erase(k);
						break;
					}
				}
				if (tag == 0){
					//Èç¹ûtagÊÇ0ÎÒÃÇÐèÒªœøÐÐALLVctœøÐÐŸ«×ŒÑéÖ€
					for (auto k = 0; k < Allvct[j].bucketVct.size(); ++k){
						auto a = temp.locVct.back().pos;
						auto b = Allvct[j].bucketVct[k].pos;
						if (a < b)
						{
							continue;
						}
						if (b - a >= kmerLenth - 2*kmere&&b - a <= kmerLenth +2* kmere){
							locbp LOC(Allvct[j].bucketVct[k].pos);
							LOC.isperfectmap = 1;
							temp.locVct.push_back(LOC);
							++temp.MappNum;
							++temp.EndSeed;
							tag = 1;
							break;
						}
					}
				}
				if (tag == 0){
					//Èç¹ûtag»¹ÊÇ0,ÄÇÃŽÎÒÃÇŸÍÒªÑ¡ÔñœøÐÐsim±È¶Ô
					auto Pos = temp.locVct.back().pos + kmerLenth;
					string refKmer = Ref.substr(Pos, kmerLenth);
					string seedKmer = seedVct[j];
					if (edistance(seedKmer, refKmer) <= kmere)
					{
						//Èç¹û±àŒ­ŸàÀëÐ¡ÓÚe,ÎÒÃÇŸÍÅÐ¶šÕâžöseedÊÇÒ»žöÏàËÆseed
						locbp LOC(Pos);
						LOC.isperfectmap = 0;
						temp.locVct.push_back(LOC);
						temp.SimNum++;
						++temp.EndSeed;
						tag = 1;
					}
				}
				if (tag == 0)break;//ËµÃ÷ºóÏò±È¶ÔÖÕÖ¹
			}
			if (temp.startSeed == 1)temp.start = 1;
			if (temp.EndSeed == MAX)temp.end = 1;
			MFVct.push_back(temp);
		}
	}
	//È¥ÖØžŽµÄMFÇÒ±»°üº¬µÄMF
	for (auto ite = MFVct.begin(); ite != MFVct.end();++ite){
		for (auto ite2 = MFVct.begin(); ite2 != MFVct.end();)
		{
			if (ite == ite2){ ++ite2; continue; }
			if (ite->locVct[0].pos <= ite2->locVct[0].pos&&ite->locVct.back().pos >= ite2->locVct.back().pos)
			{
				ite2 = MFVct.erase(ite2); continue;
			}
			++ite2;
		}
	}
	
	//³¢ÊÔºÏ²¢ÓÐœ»Œ¯µÄMF
	for (auto ite = MFVct.begin(); ite != MFVct.end(); ++ite)
	{
		//Í·
		if (ite->start == 1 && ite->end == 0)
		{//Èç¹ûÍ·²¿ºÃ£¬Î²²¿²»ºÃ£¬ÎÒÃÇŸÍ³¢ÊÔœøÐÐœ»Œ¯ºÏ²¢¹€×÷
			for (auto ite2 = MFVct.begin(); ite2 != MFVct.end();)
			{
				if (ite2->startSeed > ite->EndSeed)
				{
					++ite2; continue;
				}//seedÎ»ÖÃÒ²Òªœ»²æ
				bool b1 = ite2->start == 0;//Ñ¡³ö²»œÓÎ²µÄMF
				bool b2 = ite2->locVct[0].pos > ite->locVct[0].pos&&ite2->locVct[0].pos < ite->locVct.back().pos;//ÅÐ¶ÏÊÇ·ñÏàœ»£¬ÓÉÓÚÖ®Ç°ÒÑŸ­È¥³ýÁË°üº¬ºÍÖØžŽ£¬ÒòŽËÎÒÃÇÕâÀï²»ÓÃÅÐ¶Ï×ó±ßpos
				int pos1 = ite2->locVct[0].pos;
				int pos2 = ite->locVct[ite2->startSeed - ite->startSeed].pos;
				int distance = abs(pos1 - pos2);
				bool b3 = distance < kmerLenth;
				if (b1&&b2&&b3)
				{//Èç¹ûÓÐÒ»žö·ÇÍ·œÖÆ¬¶ÎÓëÒ»žöÍ·œÓÆ¬¶ÎÓÐœ»Œ¯£¬ËµÃ÷Æ¬¶Î¿Ï¶šÊÇÄÜºÏ²¢µÄ
					//1.ÏÈ°Ñite¶àÓàµÄÆ¬¶ÎÉŸ³ý
					auto itetemp1 = ite->locVct.begin();
					auto seedNum = ite->startSeed;
					while (1)
					{
						if (seedNum == ite2->startSeed || itetemp1 == ite->locVct.end())break;
						++seedNum;
						++itetemp1;
					}
					ite->locVct.erase(itetemp1, ite->locVct.end());
					//2.²åÈëºóÃæµÄÆ¬¶Î
					ite->locVct.insert(ite->locVct.end(), ite2->locVct.begin(), ite2->locVct.end());
					//3.ÉŸ³ýºóÃæµÄÆ¬¶Î
					ite->EndSeed = ite2->EndSeed;
					ite2 = MFVct.erase(ite2);
					//4.ÐÞžÄÐÅÏ¢
					unsigned mapnum = 0, simnum = 0;
					for (auto var : ite->locVct)
					{
						if (var.isperfectmap == 1)++mapnum;
						else ++simnum;
					}
					ite->MappNum = mapnum; ite->SimNum = simnum;
					if (ite->EndSeed == MAX)ite->end = 1;
					ite2 = MFVct.begin();
					continue;
				}
				++ite2;
			}
		}//Í·Íê
		 //Î²
		if (ite->start == 0 && ite->end == 1){
					for (auto ite2 = MFVct.begin(); ite2 != MFVct.end();)
					{
						if (ite->startSeed>ite2->EndSeed)
						{
							++ite2; continue;
						}//seedÎ»ÖÃÒ²Òªœ»²æ
						bool b1 = ite2->end == 0;//Ñ¡³ö²»œÓÎ²µÄMF
						bool b2 = ite2->locVct.back().pos > ite->locVct[0].pos&&ite2->locVct.back().pos < ite->locVct.back().pos;//ÅÐ¶ÏÊÇ·ñÏàœ»£¬ÓÉÓÚÖ®Ç°ÒÑŸ­È¥³ýÁË°üº¬ºÍÖØžŽ£¬ÒòŽËÎÒÃÇÕâÀï²»ÓÃÅÐ¶Ï×ó±ßpos
						int pos1 = ite2->locVct.back().pos;
						int pos2 = ite->locVct[ite2->EndSeed - ite->startSeed].pos;
						int distance = abs(pos1-pos2);
						bool b3 = distance < kmerLenth;
						if (b1&&b2&&b3)
						{//Èç¹ûÓÐÒ»žö·ÇÍ·œÖÆ¬¶ÎÓëÒ»žöÍ·œÓÆ¬¶ÎÓÐœ»Œ¯£¬ËµÃ÷Æ¬¶Î¿Ï¶šÊÇÄÜºÏ²¢µÄ
							//ÔÚœøÐÐÒ»ŽÎÅÐ¶Ï
							//1.ÏÈ°Ñite¶àÓàµÄÆ¬¶ÎÉŸ³ý
							auto itetemp1 = ite->locVct.begin();
							auto seedNum = ite->startSeed;
							while (1)
							{
								if (seedNum == ite2->EndSeed || itetemp1 == ite->locVct.end())break;
								++seedNum;
								++itetemp1;
							}
							ite->locVct.erase(ite->locVct.begin(), itetemp1+1);
							//2.²åÈëÇ°ÃæµÄÆ¬¶Î
							ite->locVct.insert(ite->locVct.begin(), ite2->locVct.begin(), ite2->locVct.end());
							//3.ÉŸ³ýÇ°ÃæµÄÆ¬¶Î
							ite->startSeed = ite2->startSeed;
							ite2 = MFVct.erase(ite2);
							//4.ÐÞžÄÐÅÏ¢
							unsigned mapnum = 0, simnum = 0;
							for  (auto var : ite->locVct)
							{
								if (var.isperfectmap == 1)++mapnum;
								else ++simnum;
							}
							ite->MappNum = mapnum; ite->SimNum = simnum;
							if (ite->startSeed == 1)ite->start = 1;
							ite2 = MFVct.begin();
							continue;
						}
						++ite2;
				}
			}
		//Î²Íê
	}
	//œ»Œ¯ºÏ²¢Íê
	
	//³¢ÊÔºÏ²¢MF
	for (auto ite = MFVct.begin(); ite != MFVct.end(); ++ite){
		if (ite->start == 1 && ite->end == 0)
		{
			for (auto ite2 = MFVct.begin(); ite2 != MFVct.end();)
			{
				
				if (ite2->start == 0)
				{
					bool b1 = ite->EndSeed + 2 == ite2->startSeed;
					bool b = ite->EndSeed + 1 == ite2->startSeed;
					if (b1 == 0&&b==0){ ++ite2; continue; }
					int time;
					if (b1 == 1)time = 2;
					if (b == 1)time = 1;
					int pos1 = (int)ite->locVct.back().pos;
					int pos2 = (int)ite2->locVct.begin()->pos;
					if (pos2 < pos1){ ++ite2; continue; }
					bool b2 = pos2 + e >= time * kmerLenth + pos1&&pos2 - pos1 <= time * kmerLenth + e;
					if (b2 == 0){ ++ite2; continue; }
					//¿ÉÒÔºÏ²¢
					if (b1 == 1){
						ite->EndSeed = ite2->EndSeed;
						locbp temp(pos1 + kmerLenth);
						temp.isperfectmap = 0;
						ite->locVct.push_back(temp);
						for (int i = 0; i < ite2->locVct.size(); ++i)
							ite->locVct.push_back(ite2->locVct[i]);
						ite->MappNum += ite2->MappNum;
						ite->SimNum += ite2->SimNum + 1;
						if (ite->EndSeed == MAX)ite->end = 1;
						ite2 = MFVct.erase(ite2);
						
						ite2 = MFVct.begin();
					}
					if (b==1)
					{
						ite->EndSeed = ite2->EndSeed;
						for (int i = 0; i < ite2->locVct.size(); ++i)
							ite->locVct.push_back(ite2->locVct[i]);
						ite->MappNum += ite2->MappNum;
						ite->SimNum += ite2->SimNum;
						if (ite->EndSeed == MAX)ite->end = 1;
						ite2 = MFVct.erase(ite2);
						ite2 = MFVct.begin();
					}
				}
				else
				{
					++ite2; continue;
				}

			}
		}
		if (ite->start == 0 && ite->end == 1)
		{
			for (auto ite2 = MFVct.begin(); ite2 != MFVct.end();)
			{
				if (ite2->end == 0)
				{
					bool b1 = ite2->EndSeed + 2 == ite->startSeed;
					bool b = ite2->EndSeed + 1 == ite->startSeed;
					if (b1 == 0&&b==0){ ++ite2; continue; }
					int time;
					if (b1 == 1)time = 2;
					if (b == 1)time = 1;
					int pos1 = (int)ite2->locVct.back().pos;
					int pos2 = (int)ite->locVct.begin()->pos;
					if (pos2 < pos1){ ++ite2; continue; }
					bool b2 = pos2 + e >= time * kmerLenth + pos1&&pos2 - pos1 <= time * kmerLenth + e;
					if (b2 == 0){ ++ite2; continue; }
					//¿ÉÒÔºÏ
					if (b1 == 1){
						ite->startSeed = ite2->startSeed;
						locbp temp(pos1 + kmerLenth);
						temp.isperfectmap = 0;
						ite->locVct.insert(ite->locVct.begin(), temp);
						for (int i = ite2->locVct.size() - 1; i >= 0; --i)
							ite->locVct.insert(ite->locVct.begin(), ite2->locVct[i]);
						ite->MappNum += ite2->MappNum;
						ite->SimNum += ite2->SimNum + 1; 
						if (ite->startSeed == 1)ite->start = 1;
						ite2 = MFVct.erase(ite2);
						ite2 = MFVct.begin();
					}
					if (b == 1){
						ite->startSeed = ite2->startSeed;
						for (int i = ite2->locVct.size() - 1; i >= 0; --i)
							ite->locVct.insert(ite->locVct.begin(), ite2->locVct[i]);
						ite->MappNum += ite2->MappNum;
						ite->SimNum += ite2->SimNum;
						if (ite->startSeed == 1)ite->start = 1;
						ite2 = MFVct.erase(ite2);
						ite2 = MFVct.begin();
					}
				}
				else
				{
					++ite2; continue;
				}

			}
		}	
	}
	//ÉŸ³ý²»·ûºÏ¹æ¶šµÄMF
	for (auto ite = MFVct.begin(); ite != MFVct.end();){
		if ((ite->start == 0 && ite->end == 0)) ite = MFVct.erase(ite);//ÉŸ³ýÁËsimseedÊýÁ¿Òª>eµÄÌõŒþ
		else ++ite;
	}
}
//Çó±àŒ­ŸàÀë
int GenHash::edistance(const string source, const string target){
	//simhashÏÈÅÐ¶Ï
	int hammingdistance=HammingDistance(simhash(source), simhash(target));
	if (hammingdistance < hammingE) return 0;
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
//ÏÂÃæÈýžöÊÇÕÒÖÐŒä¶ÏµãµÄÈýÖÖ¹€Ÿßº¯Êý
unsigned GenHash::simplefindBP(BPinfo& temp, bool direct, string kmerRead, string refRead, unsigned pos){
	
	if (direct == 0){
		unsigned count = 0;
		//ËµÃ÷ÓÒ±ßµÄkmerÊÇÍêÃÀ±È¶ÔµÄ,ÐèÒªÏò×óÑÓÉì
		for (auto i = kmerLenth - 1; i >= 0; --i)
		{
			if (kmerRead[i] != refRead[i])break;
			++count;
			if (i == 0)break;
		}
		return count;
	}
	else
	{
		//ËµÃ÷×ó±ßµÄKmerÊÇÍêÃÀ±È¶ÔµÄ,ÐèÒªÓÒÑÓÉì
		unsigned i = 0;
		for (; i < kmerLenth; ++i)
		if (kmerRead[i] != refRead[i])break;
		return i;
	}
}
int GenHash::shortfindBP(BPinfo& temp,  string kmerRead ,string leftRead, string rightRead,unsigned kmere){
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
int GenHash::longfindBP(BPinfo& temp, string kmerRead, string leftRead, string rightRead,unsigned kmere){
	unsigned low = 0, high = kmerRead.size() - 1;
	unsigned bppos=0;
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
		//Î»ÖÃ·¢ÉúÁËœ»Ží,¶ÏµãÖ±œÓÅÐ¶Ï
		bppos = high + (low - high) / 2;
		return bppos;
	}
	else//Ã»·¢Éúœ»Ží,ËµÃ÷ÀïÃæÓÐSNPžÉÈÅ,ËõÐ¡·¶Î§ºóÎÒÃÇŸÍ×ª¶øµ÷ÓÃshortfindbpŸÍÑ°ÕÒÒ»žöÍêÃÀµÄ¶Ïµã
	{
		bppos = low;
		int t=shortfindBP(temp, kmerRead.substr(low, high - low+1 ), leftRead.substr(low, high - low+1 ), rightRead.substr(low, high - low+1 ),kmere);
		if (t==-1)
		{
			return  t;	 
		}
		return bppos+t;
	}
}
unsigned long GenHash::simhash(string s){
	unsigned long value = 0;
	int len = s.size();
	int array[64] = { 0 };
	unsigned long tmp = 1;
	unsigned long md_value = 0;
	int i, j, k;
	k = 0;
	for (j = 0; j<len; j++){
		switch (s[j]){
		case 'A':k = 0; break;
		case 'T':k = 1; break;
		case 'C':k = 2; break;
		case 'G':k = 3; break;
		case 'N':k = 4; break;
		}
		md_value = md5table[j][k];

		for (i = 0; i<64; i++){
			tmp = 1;
			if (md_value & (tmp << i))
				array[i]++;
			else
				array[i]--;
		}

	}
	for (i = 0; i<64; i++){
		if (array[i]>0){
			tmp = 1;
			value |= (unsigned long)(tmp << (63 - i));
		}
	}
	return value;
}
unsigned GenHash::char_to_unsigned(char c){
	unsigned tmp = 0;
	switch (c){
	case '0':tmp = 0; break;
	case '1':tmp = 1; break;
	case '2':tmp = 2; break;
	case '3':tmp = 3; break;
	case '4':tmp = 4; break;
	case '5':tmp = 5; break;
	case '6':tmp = 6; break;
	case '7':tmp = 7; break;
	case '8':tmp = 8; break;
	case '9':tmp = 9; break;
	case 'a':tmp = 10; break;
	case 'b':tmp = 11; break;
	case 'c':tmp = 12; break;
	case 'd':tmp = 13; break;
	case 'e':tmp = 14; break;
	case 'f':tmp = 15; break;

	}
	return tmp;

}
unsigned long GenHash::finall(char* src){
	unsigned long tmp = 0;
	int i;
	for (i = 0; i<16; i++)
		tmp = tmp * 16 + char_to_unsigned(src[i]);
	return tmp;
}
void GenHash::middle(const unsigned char* src, char* des){
	int i;
	for (i = 0; i<16; i++)
		des[i] = src[i + 8];

}
void GenHash::mymd5(unsigned const char* src, unsigned long* des){
	unsigned char md[16];
	int i;
	char tmp[3] = { '\0' };
	char buf[33] = { '\0' };
	MD5(src, strlen((const char*)src), md);
	for (i = 0; i<16; i++){
		sprintf(tmp, "%2.2x", md[i]);
		strcat(buf, tmp);
	}
	//	printf("%s\n",buf);
	char mid[33] = { '\0' };
	middle(src, mid);
	for (i = 4; i<12; i++){
		sprintf(tmp, "%2.2x", md[i]);
		strcat(mid, tmp);
	}
	//printf("%s\n",mid);
	*des = finall(mid);

}

void GenHash::init_cache(unsigned long cache[1000][5]){
	int i;
	for (i = 0; i<1000; i++){
		char src[4];
		sprintf(src, "%d", i);
		char a[5] = { 0 };
		char t[5] = { 0 };
		char c[5] = { 0 };
		char g[5] = { 0 };
		char n[5] = { 0 };
		sprintf(a, "%c", 'A');
		sprintf(t, "%c", 'T');
		sprintf(c, "%c", 'C');
		sprintf(g, "%c", 'G');
		sprintf(n, "%c", 'N');

		strcat(a, src);
		strcat(t, src);
		strcat(c, src);
		strcat(g, src);
		strcat(n, src);

		unsigned const char* aa = (unsigned const char*)a;
		unsigned const char* tt = (unsigned const char*)t;
		unsigned const char* cc = (unsigned const char*)c;
		unsigned const char* gg = (unsigned const char*)g;
		unsigned const char* nn = (unsigned const char*)n;
		mymd5(aa, &cache[i][0]);
		mymd5(tt, &cache[i][1]);
		mymd5(cc, &cache[i][2]);
		mymd5(gg, &cache[i][3]);
		mymd5(nn, &cache[i][4]);

	}
}

int HammingDistance(unsigned long l1, unsigned long l2){
	unsigned long x = l1^l2;
	int count = 0;
	while (x){
		x = x&(x - 1);
		++count;
	}
	return count;
}

