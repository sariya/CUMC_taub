## Date 06 15 2021

CC=gcc
CFLAGS= -Wall -Wextra -Wpedantic

default=create_group_gmmat

create_group: main_groups_gmmat.o functions_groups_gmmat.o
	$(CC) $(CFLAGS) -o $(default) main_groups_gmmat.o functions_groups_gmmat.o

functions_groups_gmmat.o: functions_groups_gmmat.h functions_groups_gmmat.c

	$(CC) $(CFLAGS) -c functions_groups_gmmat.c

main_groups_gmmat.o: main_groups_gmmat.c
	$(CC) $(CFLAGS) -c main_groups_gmmat.c

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

