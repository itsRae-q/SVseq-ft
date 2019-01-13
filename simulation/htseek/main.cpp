#include"tool.h"
#include<time.h>
int main(int argc,char *argv[]){
	//read  parameter
	string refname;
        string csvname;
	refname=argv[1];
	csvname=argv[2];
        for(int i=3;i<argc;i+=2){
	     switch(argv[i][1]){
               case 'k':parameter::K=stoi(string(argv[i+1]));break;
               case 'b':parameter::BIASRATIO=stod(string(argv[i+1]));break;
               case 's':parameter::SNPRATIO=stod(string(argv[i+1]));break;
	       case 'l':parameter::KMERLENGTH=stoi(string(argv[i+1]));break;
	       case 'p':parameter::KSNP=stoi(string(argv[i+1]));break;
	       case 'o':parameter::output=argv[i+1];break;
	       case 'h':cout<<"samfile.sam  reference.sam  [-kbso]"<<endl;
			cout<<"-k distance between SNV[10]"<<endl;
                        cout<<"-b biasratio[0.3]"<<endl;
			cout<<"-s varratio[0.5]"<<endl;
			cout<<"-l kmerlength[12]"<<endl;
			cout<<"-o outfilename[FINDRESULT.txt]"<<endl;
			break;
               default:cout<<"parameter error.please -h read instruction"<<endl;
            }
	}
 	//read end
	time_t start_time = clock();
	
	vector<csv> cs;
	vector<BPinfo> bpinfo;
	string ref;
	//cout<<"input ref name"<<endl;
	//cin>>refname;  
	//cout<<"input sam name"<<endl;
	//cin>>samname;
        // unsigned K;
	////cout<<"input k"<<endl;
	//cin>>K;
	//parameter::K=K;
	long countN=readRef(ref,refname);
	write_csv(cs, csvname);
	getBPinfoVct(cs,bpinfo,ref,start_time,countN);
        
	cout<<endl<<parameter::KSNP<<endl;
	cout<<parameter::K<<endl;
	cout<<parameter::KMERLENGTH<<endl;
	return 0;
}
