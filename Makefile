stencil: stencil.c
	gcc -g -pg -static-libgcc -std=c99 -Wall $^ -o $@
