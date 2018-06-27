CC=gcc
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -lm
BIN=./bin/
SRC=./src/c/

all: mergeMultiFasta uniseqDBCoverage format qnuclparserblast slen taxomaker filterParsedBlast

mergeMultiFasta: $(SRC)mergeMultiFasta.c
	$(CC) $(CFLAGS) $(SRC)mergeMultiFasta.c -o $(BIN)mergeMultiFasta

uniseqDBCoverage: $(SRC)uniseqDBCoverage.c
	$(CC) $(CFLAGS) $(SRC)uniseqDBCoverage.c -o $(BIN)uniseqDBCoverage

format: $(SRC)format.c
	$(CC) $(CFLAGS) $(SRC)format.c -o $(BIN)format

qnuclparserblast: $(SRC)qnuclparserblast.c
	$(CC) $(CFLAGS) $(SRC)qnuclparserblast.c -o $(BIN)qnuclparserblast

slen: $(SRC)slen.c
	$(CC) $(CFLAGS) $(SRC)slen.c -o $(BIN)slen

taxomaker: $(SRC)taxomaker.c
	$(CC) $(CFLAGS) $(SRC)taxomaker.c -o $(BIN)taxomaker

filterParsedBlast: $(SRC)filterParsedBlast.c
	$(CC) $(CFLAGS) $(SRC)filterParsedBlast.c -o $(BIN)filterParsedBlast
