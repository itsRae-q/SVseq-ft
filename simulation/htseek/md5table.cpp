#include"md5table.h"
unsigned char_to_unsigned(char c){
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
unsigned long finall(char* src){
	unsigned long tmp = 0;
	int i;
	for (i = 0; i<16; i++)
		tmp = tmp * 16 + char_to_unsigned(src[i]);
	return tmp;
}
void middle(const unsigned char* src, char* des){
	int i;
	for (i = 0; i<16; i++)
		des[i] = src[i + 8];

}
void mymd5(unsigned const char* src, unsigned long* des){
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

void init_cache(unsigned long cache[1000][5]){
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

