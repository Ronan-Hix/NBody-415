all: NBody
NBody: NBody.c nrutil.c rk4.c
	gcc NBody.c rk4.c nrutil.c -o NBody -lm
clean:
	rm NBody
	rm *.out
