# Makefile

all: gtf ffa
	
gtf: gtf.cpp graph.cpp graphtv.cpp
	g++ -std=c++11 -o gtf gtf.cpp graph.cpp graphtv.cpp

ffa: ffa.cpp
	g++ -std=c++11 -o ffa ffa.cpp

clean:
	rm -f gtf ffa