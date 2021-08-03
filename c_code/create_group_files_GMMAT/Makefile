## Date 06 15 2021

CC=g++
CFLAGS= -Wall -Wextra -Wpedantic -O3

default=create_group_gmmat

create_group: main_groups_gmmat.o functions_groups_gmmat.o
	$(CC) $(CFLAGS) -o $(default) main_groups_gmmat.o functions_groups_gmmat.o

functions_groups_gmmat.o: functions_groups_gmmat.hpp functions_groups_gmmat.cpp

	$(CC) $(CFLAGS) -c functions_groups_gmmat.cpp

main_groups_gmmat.o: main_groups_gmmat.cpp
	$(CC) $(CFLAGS) -c main_groups_gmmat.cpp

.PHONY: clean ### https://stackoverflow.com/a/2145605/2740831

clean: 
	$(RM) $(default) *.o 

## all: $(TARGET)
## ##$(RM) gtf_parse *.o *~

##.PHONY: clean

##clean:
##	rm -f $(TARGET) *~ core

