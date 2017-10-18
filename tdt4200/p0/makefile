TARGET=julia
CC=gcc
CFLAGS=-Wall -std=c99
OBJECTS=julia.o bitmap.o


all: julia
julia: $(OBJECTS)
julia.o: julia.c julia.h bitmap.o
bitmap.o: bitmap.c bitmap.h

clean:
	rm -f julia *.o
run: julia
	./julia 
