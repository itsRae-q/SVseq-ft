#ifndef __CONFIGINFO__
#define __CONFIGINFO__
class configinfo
{
public:
	configinfo(unsigned int low, unsigned int high, unsigned int count) :low(low), high(high), count(count)
	{}
	unsigned int low;
	unsigned int high;
	unsigned int count;
};
#endif
