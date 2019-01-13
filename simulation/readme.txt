Process using SVseq-ft:
	1.Use “add_com” to add mutations to ref.fasta, and the mutation information is specified in the config.txt file. 
	  The algorithm is based on c++11 under the windows platform.
	  
	2.The simulation and comparison process of the genetic data was completed using PBsim and Blasr to obtain the .sam file. 
	  The commands used by PBsim and Blasr are as follows:
	  --------------------------（pbsim和blasr命令）
	  
	3.The identification and sequence correction of the variant region is done using SVseq-ft (implemented in Python based on Python 3.6) where we will get csv file.
	  The command is as follows:
	           python connect.py -s *.sam
	
	4.Using htseek to complete the location and type determination of the original location of the complex indel, HTseek is based on C++11 under linux. 
	  HTseek uses the command as follows:
	  ./HTseek ref.fasta *.sam_csv -k  -l  -p  -o   where 'k' is the length of the seed
	  
	
	  