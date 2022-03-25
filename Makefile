CC=gcc
CFLAGS=-W -Wall -O3
LIB=-lm

all: edges2k acqplot7amoon longav

edges2k: src/edges2k.c
	$(CC) $(CFLAGS) -o bin/$@ $^ $(LIB)

acqplot7amoon: src/acqplot7amoon.c
	$(CC) $(CFLAGS) -o bin/$@ $^ $(LIB)

longav: src/longav.c
	$(CC) $(CFLAGS) -o bin/$@ $^ $(LIB)
