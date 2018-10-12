# gtf/Makefile

all :
	mkdir -p bin
	make ffa
	make gtf

ffa : ffagtf/min.cpp
	g++ -std=c++11 -o bin/ffa $^

gtf : graphtv/gtf.cpp graphtv/graph.cpp graphtv/graphtv.cpp
	g++ -std=c++11 -o bin/gtf $^

clean :
	rm -rf bin
