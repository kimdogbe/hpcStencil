stencil: stencil.c
	gcc -O2 -g -pg -static-libgcc -std=c99 -Wall $^ -o $@
