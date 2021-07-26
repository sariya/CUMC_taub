## Date 07 26 2021

CC=gcc
CFLAGS= -Wall -Wextra -Wpedantic

default=saige_create_group_annovar

saige_create_group_annovar: main_create_groupsfile.o functions_create_saige.o
	$(CC) $(CFLAGS) -o $(default) functions_create_saige.o main_create_groupsfile.o

functions_create_saige.o: functions_create_saige.h functions_create_saige.c
	$(CC) $(CFLAGS) -c functions_create_saige.c

main_create_groupsfile.o: main_create_groupsfile.c
	$(CC) $(CFLAGS) -c main_create_groupsfile.c

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

