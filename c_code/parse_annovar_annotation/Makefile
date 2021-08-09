CC=g++
CFLAGS= -Wall -Wextra -Wpedantic -O3

default=parse_annovar_annot_GCTA

parse_annovar_annot_GCTA: main_parse_annotation.o parse_function.o
	$(CC) $(CFLAGS) -o $(default) parse_function.o main_parse_annotation.o

parse_function.o: parse_function.hpp parse_function.cpp
	$(CC) $(CFLAGS) -c parse_function.cpp

main_parse_annotation.o: main_parse_annotation.cpp
	$(CC) $(CFLAGS) -c main_parse_annotation.cpp

.PHONY: 
	clean 

clean: 
	$(RM) $(default) *.o 

## all: $(TARGET)
#### $(RM) gtf_parse *.o *~

##clean:
## ### https://stackoverflow.com/a/2145605/2740831
##	rm -f $(TARGET) *~ core
