CC=gcc
CFLAGS=-W -Wall -O3
LIB=-lm

all: edges2k acqplot7amoon longav edges2k_fittpfix edges3 edges3_fittpfix reads1p1 corrcsv

edges2k: src/edges2k.c
	$(CC) $(CFLAGS) -o bin/$@ $^ $(LIB)

edges2k_fittpfix: src/edges2k.c
	$(CC) -DFITTP_FIX $(CFLAGS) -o bin/$@ $^ $(LIB)

edges3: src/edges3.c
	$(CC) $(CFLAGS) -o bin/$@ $^ $(LIB)

edges3_fittpfix: src/edges3.c
	$(CC) -DFITTP_FIX $(CFLAGS) -o bin/$@ $^ $(LIB)

acqplot7amoon: src/acqplot7amoon.c
	$(CC) $(CFLAGS) -o bin/$@ $^ $(LIB)

longav: src/longav.c
	$(CC) $(CFLAGS) -o bin/$@ $^ $(LIB)

reads1p1: src/reads1p1.c
	$(CC) $(CFLAGS) -o bin/$@ $^ $(LIB)

corrcsv: src/corrcsv.c
	$(CC) $(CFLAGS) -o bin/$@ $^ $(LIB)