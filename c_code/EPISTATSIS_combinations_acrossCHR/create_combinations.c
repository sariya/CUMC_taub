#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "funcs_combinations.h"
#include "common_include_file.h"
#include "message_help.h"

//////////////////////////////
int main(int argc, char *argv[])
{
    FILE *fptr_first, *fptr_second; // use this for opening files first and second gene names
    char read_line[MAX_LEN];
    char chr_first[6] = "CHR";
    char chr_second[6] = "CHR";
    int line_number = 0; //use this for counter of header

    // ./combinate_across  small_CHR11.txt 11 small_CHR8.txt 8

    data_store *start_one = NULL;
    data_store *start_second = NULL; //hold gene names of other CHR
    //    print_help(argc, argv);

    if (argc != 5)
    {
        printf("incorrect number of arguments are %d\n", argc);
        printf("./combinate_across small_CHR11.txt 11 small_CHR8.txt 8\n");
        return -1;
    }

    concatenate_str(chr_first, argv[2]);
    concatenate_str(chr_second, argv[4]);

    fptr_first = fopen(argv[1], "r");

    if (fptr_first == NULL)
    {
        printf("null pointer for file opening first gene file\n");
        return -1;
    }

    while (fgets(read_line, MAX_LEN, fptr_first))
    {
        if (line_number > 0) // make sure the line excluded header
        {
            if (strcmp(read_line, "\n") != 0)
            {
                remove_trailingspaces(read_line);
                get_gene_name(read_line);
                add_node(&start_one, read_line); //only if not new line
            }
        }
        line_number++;
    }
    fclose(fptr_first);

    fptr_second = fopen(argv[3], "r");
    if (fptr_second == NULL)
    {
        printf("null pointer for file opening second gene file\n");
        return -1;
    }
    line_number = 0;
    while (fgets(read_line, MAX_LEN, fptr_second))
    {
        if (line_number > 0) // make sure the line excluded header
        {
            if (strcmp(read_line, "\n") != 0)
            {
                remove_trailingspaces(read_line);
                get_gene_name(read_line);
                add_node(&start_second, read_line); //only if not new line
            }
        }
        line_number++;
    }
    fclose(fptr_second);

    printf("we have so many genes %d\n", check_length(start_one));
    printf("we have so many genes %d\n", check_length(start_second));
    // display_list(start_one);
    // display_list(start_second);
    printf("The new chrs are %s\n", chr_first);
    printf("The new chrs are %s\n", chr_second);
    create_combination(start_one, start_second,chr_first, chr_second ,argv[5]);
    delete_linked_list(&start_one);
    delete_linked_list(&start_second);

    return 0;
}
/////////////////////////////////////////