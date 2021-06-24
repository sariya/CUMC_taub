TARGET=combinate_across
INCDIR=.
CC=gcc
CFLAGS= -Wall -Wextra -Wpedantic

all: $(TARGET)

$(TARGET):combinate_across
	$(CC) $(CFLAGS) -o $(TARGET) *.c

.PHONY: clean

clean:
	rm -f $(TARGET) *~ core

