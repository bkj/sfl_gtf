# Makefile

all: main
	
main: main.cpp graph.cpp graphtv.cpp
	g++ -std=c++11 -o main main.cpp graph.cpp graphtv.cpp
	
clean:
	rm -f main