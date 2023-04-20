
CC  := gcc

CFLAGS=-Wall -O3 -std=c11  
CXXFLAGS=$(CFLAGS) -std=c++11
LIBS=-lz
PROG=tod

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

all:$(PROG)


tod:count.c merge.c map.c set.c  kvec.h khashl.h ketopt.h kseq.h kthread.h kstring.h ksort.h misc.h 
	$(CC) $(CFLAGS) -o $@ main.c  count.c merge.c match.c map.c set.c kthread.c kstring.c  $(LIBS) -lpthread

clean:
	rm -fr *.dSYM $(PROG)

