a.out: mc_chain.o
	gcc -L/usr/local/lib mc_chain.o -lgsl -lgslcblas -lm -o mc_chain

mc_chain.o: mc_chain.c
	gcc -Wall -I/usr/local/include -c mc_chain.c

clean:
	rm -f *.o

clean-all: clean
	rm -f mc_chain
