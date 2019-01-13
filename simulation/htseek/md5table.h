# pragma once
#include<openssl/md5.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
using namespace std;
unsigned char_to_unsigned(char c);
unsigned long finall(char* src);
void middle(const unsigned char* src, char* des);
void mymd5(unsigned const char* src, unsigned long* des);
void init_cache(unsigned long cache[1000][5]);



