CC=g++
CFLAGS=-O2 -std=c++0x

nps: grs stdgt

grs: src/genetic_score.bimbam.cpp
	$(CC) $(CFLAGS) -o grs src/genetic_score.bimbam.cpp 

stdgt: src/standardize_gt.bimbam.cpp
	$(CC) $(CFLAGS) -o stdgt src/standardize_gt.bimbam.cpp

clean:
	rm grs stdgt
