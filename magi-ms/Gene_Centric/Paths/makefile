CC=g++ -mcmodel=large
CFLAGS= -c -O3

Pathway_GeneCenter: color_coding.o PPI_Graph.o
	$(CC) -O3 -o Pathway_GeneCenter color_coding.o PPI_Graph.o

color_coding.o: color_coding.cpp color_coding.h
	$(CC) $(CFLAGS) color_coding.cpp

PPI_Graph.o: PPI_Graph.cpp PPI_Graph.h
	$(CC) $(CFLAGS) PPI_Graph.cpp
