TARGET=print
INCDIR=.
CC=gcc
CFLAGS= -Wall -Wextra -Wpedantic

all: $(TARGET)

$(TARGET): hello
	$(CC) $(CFLAGS) -o $(TARGET) hello.c

.PHONY: clean

clean:
	rm -f $(TARGET) *~ core

