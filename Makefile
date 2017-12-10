a.out: example.o
	gcc -L/usr/local/lib example.o -lgsl -lgslcblas -lm

example.o: example.c
	gcc -Wall -I/usr/local/include -c example.c

clean:
	rm -f *.o

clean-all: clean
	rm -f a.out
