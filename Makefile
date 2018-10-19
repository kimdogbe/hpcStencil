stencil: stencil.c
	gcc -O3 -g -pg -static-libgcc -std=c99 -Wall $^ -o $@
