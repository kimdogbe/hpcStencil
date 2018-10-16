stencil: stencil.c
	gcc -02 -g -pg -static-libgcc -std=c99 -Wall $^ -o $@
