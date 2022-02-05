all: PPRquery metric
PPRquery: SFMT.c main.cpp 
	g++ -march=core2 -pthread -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o PPRquery SFMT.c main.cpp
	
metric: SFMT.c metric.cpp
	g++ -march=core2 -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o metric SFMT.c metric.cpp
