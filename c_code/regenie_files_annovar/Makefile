## Date 08 02 2021

CC=g++
CFLAGS= -Wall -Wextra -Wpedantic -O3 

default=regenie_files_annovar ## final name of the executable to be used

regenie_files_annovar: functions_regenie.o create_regenie_files_main.o
	$(CC) $(CFLAGS) -o $(default) functions_regenie.o create_regenie_files_main.o

functions_regenie.o: functions_regenie.cpp
	$(CC) $(CFLAGS) -c functions_regenie.cpp

create_regenie_files_main.o: create_regenie_files_main.cpp
	$(CC) $(CFLAGS) -c create_regenie_files_main.cpp

.PHONY: clean ### https://stackoverflow.com/a/2145605/2740831

clean: 
	$(RM) $(default) *.o 
