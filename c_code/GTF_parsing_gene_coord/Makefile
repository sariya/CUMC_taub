CC=gcc
CFLAGS= -Wall -Wextra -Wpedantic

default=gtf_parse

gtf_parse: main_parse_gtf.o functions_parse_gtf.o
	$(CC) $(CFLAGS) -o $(default) functions_parse_gtf.o main_parse_gtf.o 

functions_parse_gtf.o: functions_parse_gtf.c functions_parse_gtf.h

	$(CC) $(CFLAGS) -c functions_parse_gtf.c

main_parse_gtf.o: main_parse_gtf.c
	$(CC) $(CFLAGS) -c main_parse_gtf.c

.PHONY: clean ### https://stackoverflow.com/a/2145605/2740831

clean: 
	$(RM) $(default) *.o 

	##$(RM) gtf_parse *.o *~

## all: $(TARGET)

## $(TARGET): hello
##	$(CC) $(CFLAGS) -o $(TARGET) hello.c

##.PHONY: clean

##clean:
##	rm -f $(TARGET) *~ core

