#include "randomgen.h"
default_random_engine randomgen::e(time(0));
int randomgen::randominsertiondeviation(){
	static uniform_int_distribution<signed> u(-5, 5);
	return u(e);
}
double randomgen::randomtype(){
	static uniform_real_distribution<double> u(0, 100);
	return u(e);
}
unsigned int randomgen::randomdeletionsize(unsigned min, unsigned max){
	 uniform_int_distribution<signed> u(min, max);
	return u(e);
}
unsigned int randomgen::randompos(unsigned size){
	 uniform_int_distribution<signed> u(0, size-1);
	return u(e);
}
bool randomgen::isInversion(const double ratio){
	static uniform_real_distribution<double> u(0, 1);
	double tmp = u(e);
	return tmp <=ratio;//���������ɵ�С��<=inversion�ı��� �ͷ���true
}