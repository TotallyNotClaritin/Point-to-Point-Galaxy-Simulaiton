all: galaxysim

galaxysim:

	gcc -o galaxysim.o randn.c leapfrog.c nrutil.c galaxysim.c -fopenmp -g -lm

clean:
	rm *.o


# DO NOT DELETE

