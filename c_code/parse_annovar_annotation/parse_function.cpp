
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "parse_function.hpp"

// 05 09 2021
/////////////

void print_usage()
{
    printf("Usage: ././parse_annovar_annot_GCTA -a annovar_file -o output_file\n");
}
///


void remove_trailingspaces(char *newline)
{
    /**
    * remove trailing new line character from the line 
    */
    size_t line_length = strlen(newline);

    if (line_length == 0)
    {
        printf("remove_trailingspaces issues with line trimming\n");
        return;
    }
    else
    {
        newline[line_length - 1] = '\0';
    }
}

///////////////remove_trailingspaces function ends


void assign_variables(int count_argc, char **argv_array, const char *short_array, const struct option long_opt_temp[], char **temp_annovar_file, char **temp_out_file)
{
    /**
     * We need a function to assign variables
    */

    int opt;
    while ((opt = getopt_long(count_argc, argv_array, short_array, long_opt_temp, NULL)) != -1)
    {
        switch (opt)
        {
        case 'a':
            *temp_annovar_file = optarg;
            break;

        case 'o':
            *temp_out_file = optarg;
            break;
        case 'h':
            print_usage();
            exit(0);
            //break;
        }
        //switch case ends
    }
}
////////



// ///////// search_gene Function ends