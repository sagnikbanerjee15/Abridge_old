SRCS = $(wildcard *.c)

PROGS = $(patsubst %.c,%,$(SRCS))

CFLAGS = -Ofast -g

INSTALLDIR ?= /usr/local/bin

all: $(PROGS)

%: %.c

install: $(PROGS)
	cp $(PROGS) $(INSTALLDIR)

clean: 
	rm -f $(PROGS) *.o