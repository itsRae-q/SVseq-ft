main: main.o tool.o GenHash.o 
	g++ main.o tool.o GenHash.o -lcrypto -std=c++11 -g -o HTseek
main.o: main.cpp tool.h
	g++ -lcrypto -std=c++11 -g -c main.cpp
tool.o: tool.cpp tool.h 
	g++ -lcrypto -std=c++11 -g -c tool.cpp
GenHash.o: GenHash.cpp GenHash.h md5table.h
	g++ -lcrypto -std=c++11 -g -c GenHash.cpp

clean:
	rm *.o


