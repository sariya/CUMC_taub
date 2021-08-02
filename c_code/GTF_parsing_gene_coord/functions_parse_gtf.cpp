#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "functions_parse_gtf.hpp"

void print_usage()
{
    printf("Usage: ./gtf_parse -g GTF_file -o output_file \n");
}
///Function ends

void assign_variables(int count_argc, char **argv_array, const char *short_array, const struct option long_opt_temp[], char **temp_gtf_file, char **temp_out_file)
{
    /**
     * We need a function to assign variables
    */

    int opt;
    while ((opt = getopt_long(count_argc, argv_array, short_array, long_opt_temp, NULL)) != -1)
    {
        switch (opt)
        {
        case 'g':
            *temp_gtf_file = optarg;
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
int check_chromosome(char *chr_string)
{
    /**
     * check if CHR exists 
     * We do not want patch names in the output
    */
    int flag_found_chr = 0;
    if (strcmp(chr_string, "1") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "2") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "3") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "4") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "5") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "6") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "7") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "8") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "9") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "10") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "11") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "12") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "13") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "14") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "15") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "16") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "17") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "18") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "19") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "20") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "21") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "22") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "Y") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "X") == 0)
    {
        flag_found_chr = 1;
    }
    if (strcmp(chr_string, "MT") == 0)
    {
        flag_found_chr = 1;
    }

    return flag_found_chr;
}
// ///


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

void remove_double_quotes(char **temp_gene_info)
{
    /**
     * 
    */
    if (*temp_gene_info == NULL)
    {
        printf("String received is NULL\n");
        return;
    }
    else
    {
        int str_length = strlen(*temp_gene_info); //store length
        int itr_quote = 0;                        //us to iterate and find quote

        while (itr_quote < str_length)
        {
            if (*(*temp_gene_info + itr_quote) == '"') // wait until you find a double quote.
            {
                *(*temp_gene_info + str_length - 1) = '\0';        //make last before character as NULL to remove "
                *temp_gene_info = *temp_gene_info + itr_quote + 1; //move pointer to the location after "
                break;
            }
            itr_quote++;
        }
        /// while ends
    }
}
////////Function ends
// char *mystrcat(char *dest, char *src)
// {
//     /** function to concatenate two strings. 
//      * https://stackoverflow.com/a/21881314/2740831
//     */
//     while (*dest)
//     {
//         dest++;
//     }
//     while (*dest++ = *src++)
//     {
//         ;
//     }

//     return --dest;
// }
////////////////////////////