<<<<<<< HEAD
# gtf/Makefile

all :
	mkdir -p bin
	make ffa
	make gtf

ffa : ffagtf/ffa.cpp
	g++ -std=c++11 -o bin/ffa $^

gtf : graphtv/gtf.cpp graphtv/graph.cpp graphtv/graphtv.cpp
	g++ -std=c++11 -o bin/gtf $^

clean :
	rm -rf bin
=======
# Makefile

all: gtf ffa
	
gtf: src/gtf.cpp src/graph.cpp src/graphtv.cpp
	g++ -std=c++11 -o gtf src/gtf.cpp src/graph.cpp src/graphtv.cpp -Isrc

ffa: src/ffa.cpp
	g++ -std=c++11 -o ffa src/ffa.cpp -Isrc

clean:
	rm -f gtf ffa
>>>>>>> 904f1341a060394ed0f4d1b29a1f43d7d6dd953a
