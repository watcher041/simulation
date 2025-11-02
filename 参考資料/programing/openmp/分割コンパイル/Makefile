
compile: main.o output.o
	gcc -o test -fopenmp main.o output.o
.c.o:
	gcc -c -Wall -O3 -fopenmp $<
