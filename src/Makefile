SRCS = $(wildcard *.c)

PROGS = $(patsubst %.c,%,$(SRCS))

CFLAGS = -Ofast -g -Wall

all: $(PROGS)

%: %.c

clean: 
	rm -f $(PROGS)
