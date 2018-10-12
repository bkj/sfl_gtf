# Makefile

all: gtf ffa
	
gtf: src/gtf.cpp src/graph.cpp src/graphtv.cpp
	g++ -std=c++11 -o gtf src/gtf.cpp src/graph.cpp src/graphtv.cpp -Isrc

ffa: src/ffa.cpp
	g++ -std=c++11 -o ffa src/ffa.cpp -Isrc

clean:
	rm -f gtf ffa