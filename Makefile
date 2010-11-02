CC=gcc
CFLAGS=-g -Wall

cli-testing: main.c STABIL.c
	$(CC) $(CFLAGS) -D DEBUG -o STABIL main.c STABIL.c

cli: main.c STABIL.c
	$(CC) $(CFLAGS) -o STABIL main.c STABIL.c

all: stabil1 stabil2 stabcol

stabcol: germ1.c 
	$(CC) $(CFLAGS) -o stabcol germ1.c 

stabil1: rus1.c
	$(CC) $(CFLAGS) -o stabil1 rus1.c

stabil2: rus2.c
	$(CC) $(CFLAGS) -o stabil2 rus2.c

test: all
	-./stabcol input1; ./stabil1 input1; ./stabil2 input1

clean:
	-rm stabil1 stabil2 stabcol *.o *.out *.exe

