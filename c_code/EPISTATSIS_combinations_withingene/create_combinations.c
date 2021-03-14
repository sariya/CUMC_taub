#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "funcs_combinations.h"
#include "common_include_file.h"
#include "message_help.h"

//#define MAX_LEN 500

//////////////////////////////
int main(int argc, char *argv[])
{
    FILE *fptr;
    char c[MAX_LEN];
    int line_number = 0;
    //  gcc -Wpedantic -Wextra -Wall *.c -o combinate
    //./combinate  ggggg  testCHR1_combinations.txt

    data_store *start = NULL;
//    print_help(argc, argv);

    fptr = fopen(argv[1], "r");

    if (fptr == NULL)
    {
        printf("null pointer for file opening\n");
    }

    while (fgets(c, MAX_LEN, fptr))
    {
        if (line_number > 0) // make sure the line excluded header
        {
            if (strcmp(c, "\n") != 0)
            {
                remove_trailingspaces(c);
                get_gene_name(c);
                add_node(&start, c); //only if not new line
            }
        }
        line_number++;
    }

    fclose(fptr);

    printf("we have so many genes %d\n", check_length(start));
    //display_list(start);
    create_combination(start, argv[2]);
delete_linked_list(&start);

    return 0;
}
/////////////////////////////////////////