## Date 05 17 2021

CC=gcc
CFLAGS= -Wall -Wextra -Wpedantic

default=create_group

create_group: create_group_file.o functions_group.o
	$(CC) $(CFLAGS) -o $(default) create_group_file.o functions_group.o

functions_group.o: functions_group.h functions_group.c

	$(CC) $(CFLAGS) -c functions_group.c

create_group_file.o: create_group_file.c
	$(CC) $(CFLAGS) -c create_group_file.c

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

