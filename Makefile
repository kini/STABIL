CC=gcc
## Release
#CFLAGS=-O2
## Development
CFLAGS=-g -Wall
## Debugging
#CFLAGS=-g -Wall -D DEBUG

cli: main.c STABIL.c
	$(CC) $(CFLAGS) -o STABIL main.c STABIL.c

old: rus2.c
	$(CC) $(CFLAGS) -o STABIL.old rus2.c

stabcol: germ1.c
	$(CC) $(CFLAGS) -o STABCOL germ1.c

all: cli old stabcol

test: all
	./STABIL.old 1.in; ./STABCOL 1.in; ./STABIL 1.in

clean:
	rm -f STABIL STABIL.old STABCOL *.o *.out *.exe

