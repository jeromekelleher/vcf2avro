
CFLAGS=-g -Wall -I${HOME}/.local/include
LDFLAGS=-L${HOME}/.local/lib -lavro

vcfcat: vcfcat.c
	gcc -o vcfcat ${CFLAGS} vcfcat.c ${LDFLAGS} 
