CC=g++
CFLAGS= -Wall -Wextra -Wpedantic -O3

default=gtf_parse

gtf_parse: main_parse_gtf.o functions_parse_gtf.o
	$(CC) $(CFLAGS) -o $(default) functions_parse_gtf.o main_parse_gtf.o 

functions_parse_gtf.o: functions_parse_gtf.cpp functions_parse_gtf.hpp
	$(CC) $(CFLAGS) -c functions_parse_gtf.cpp 

main_parse_gtf.o: main_parse_gtf.cpp
	$(CC) $(CFLAGS) -c main_parse_gtf.cpp 

.PHONY: 
	clean 

clean: 
	$(RM) $(default) *.o 

## all: $(TARGET)
#### $(RM) gtf_parse *.o *~

##clean:
## ### https://stackoverflow.com/a/2145605/2740831
##	rm -f $(TARGET) *~ core

## g++ -c main_parse_gtf.cpp ; g++ -c functions_parse_gtf.cpp ; g++ -o gtf_parse functions_parse_gtf.o main_parse_gtf.o