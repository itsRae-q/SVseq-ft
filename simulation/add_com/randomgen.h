#ifndef __RANDOMGENTOOL__
#define __RANDOMGENTOOL__
#include<random>
#include<ctime>
using namespace std;
class randomgen
{
public:
	double randomgen::randomtype();
	int randominsertiondeviation();
	unsigned int randomdeletionsize(unsigned min, unsigned max);
	unsigned int randompos(unsigned size);
	bool isInversion(const double ratio);
	static default_random_engine e;
};
#endif
