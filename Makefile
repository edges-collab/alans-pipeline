CC=gcc
CFLAGS=-W -Wall -O3
LIB=-lm

all: edges2k acqplot7amoon longav

edges2k: src/edges2k.c
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

acqplot7amoon: src/acqplot7amoon.c
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

longav: src/longav.c
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)
