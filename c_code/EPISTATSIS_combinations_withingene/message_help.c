#include <stdio.h>
#include <unistd.h>
#include "message_help.h"

void print_help(int variable_count, char *argument_line[])
{
    int option;

    if (variable_count < 2)
    {
        printf("use \ncombinate -h\n");
    }
    else
    {

        char message[] = "./combinate -g genefile -o outputfile\n\
            -g genefile has four columns\noutputfile is output file\n\
            -o exiting files will not be overwritten";
        while ((option = getopt(variable_count, argument_line, "g:o:")) != -1)
        { //get option from the getopt() method
            switch (option)
            {
            case 'g':
                //  printf("Giveninput file: %s\n", optarg);
                break;
            case 'o':
                //printf("Given output Option: %s\n", optarg);
                break;
            case '?': //used for some unknown options
                 printf("use as \n%s\n", message);
                break;
            }
        }
    }
}