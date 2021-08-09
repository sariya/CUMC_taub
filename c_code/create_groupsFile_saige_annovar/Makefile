## Date 07 26 2021

CC=g++
CFLAGS= -Wall -Wextra -Wpedantic -O3

default=saige_create_group_annovar

saige_create_group_annovar: main_create_groupsfile.o functions_create_saige.o
	$(CC) $(CFLAGS) -o $(default) functions_create_saige.o main_create_groupsfile.o

functions_create_saige.o: functions_create_saige.hpp functions_create_saige.cpp
	$(CC) $(CFLAGS) -c functions_create_saige.cpp

main_create_groupsfile.o: main_create_groupsfile.cpp
	$(CC) $(CFLAGS) -c main_create_groupsfile.cpp

.PHONY: clean ### https://stackoverflow.com/a/2145605/2740831

clean: 
	$(RM) $(default) *.o 

## all: $(TARGET)

## $(TARGET): hello
##	$(CC) $(CFLAGS) -o $(TARGET) hello.c

##.PHONY: clean 	##$(RM) gtf_parse *.o *~


##clean:
##	rm -f $(TARGET) *~ core

