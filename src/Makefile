CC=gcc
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64
BIN=../bin/

all: mergeMultiFasta uniseqDBCoverage format qnuclparserblast slen taxomaker filterParsedBlast

mergeMultiFasta: mergeMultiFasta.c
	$(CC) $(CFLAGS) mergeMultiFasta.c -o ../bin/mergeMultiFasta

uniseqDBCoverage: uniseqDBCoverage.c
	$(CC) $(CFLAGS) uniseqDBCoverage.c -o ../bin/uniseqDBCoverage

format: format.c
	$(CC) $(CFLAGS) format.c -o ../bin/format

qnuclparserblast: qnuclparserblast.c
	$(CC) $(CFLAGS) qnuclparserblast.c -o ../bin/qnuclparserblast

slen: slen.c
	$(CC) $(CFLAGS) slen.c -o $(BIN)slen

taxomaker: taxomaker.c
	$(CC) $(CFLAGS) taxomaker.c -o $(BIN)taxomaker

filterParsedBlast: filterParsedBlast.c
	$(CC) $(CFLAGS) filterParsedBlast.c -o $(BIN)filterParsedBlast
