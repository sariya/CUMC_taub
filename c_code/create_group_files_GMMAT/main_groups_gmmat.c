#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions_groups_gmmat.h"

// //////////////////////////////////////
int main(int argc, char *argv[])
{
    /**
     * //Date 06 15 2021
     * we will read file from the user and output it in a tab delimited manner
     * 
    */

    FILE *fp_annovar = fopen("input_annotation.txt", "r");

    if (fp_annovar == NULL)
    {
        fprintf(stderr, "Input file for annotation isn't correct. Exiting\n");
        exit(EXIT_FAILURE);
    }
    char *line = NULL;
    size_t len = 0;
    int check_line_counter = 0;

    const char *delimiters = "\t "; //use this for spliting lines

    char final_string[9000]; //ANKRD62P1-PARP4P3	chr22	17102891	G	T	1
    char string_annotation[7000]; //chr22	17102891	G	T	1

    while ((getline(&line, &len, fp_annovar)) != -1)
    {
        remove_trailingspaces(line); //remove trailing new line character

        if (check_line_counter > 1)
        {
            string_annotation[0] = '\0';
            char *token = NULL;
            int count_split = 1; //use this for counting split. as columns are fixed

            token = strtok(line, delimiters);

            while (token != NULL && count_split <= 7)
            {
                if (count_split == 1)
                {
                    mystrcat(string_annotation, "chr");
                    mystrcat(string_annotation, token);
                }
                //first column check ends

                if (count_split == 2)
                {
                    mystrcat(string_annotation, "\t");
                    mystrcat(string_annotation, token);
                }
                //second column check ends

                if (count_split == 4)
                {
                    if (strlen(token) == 1 && strcmp(token, ".") != 0)
                    {
                        mystrcat(string_annotation, "\t");
                        mystrcat(string_annotation, token);
                    }
                    else
                    {
                        //if we have an indel get away now
                        break;
                    }
                }
                //fourth column check ends
                if (count_split == 5)
                {
                    if (strlen(token) == 1 && strcmp(token, ".") != 0)
                    {
                        mystrcat(string_annotation, "\t");
                        mystrcat(string_annotation, token);
                    }
                    else
                    {
                        ///we have an indel
                        break;
                    }
                }
                //fifth column check ends

                if (count_split == 7)
                {
                    //iterate over the column with gene infomration 
                    if (strcmp(token, ".") != 0)
                    {

                        mystrcat(string_annotation, "\t");
                        mystrcat(string_annotation, "1");

                        const char *delimiters_gene = ","; //use this for spliting gene column
                        char *token_gene = NULL;
                        token_gene = strtok(token, delimiters_gene);

                        while (token_gene != NULL)
                        {
                            final_string[0] = '\0';
                            mystrcat(final_string, token_gene);
                            mystrcat(final_string, "\t");
                            //
                            mystrcat(final_string, string_annotation);
                            printf("%s\n", final_string);
                            token_gene = strtok(NULL, delimiters_gene);
                        }
                        //while tokenizing ends
                    }
                    //check ends if value is not a dot
                    else
                    {
                        //we couldn't get any gene
                        break;
                    }
                }
                //last column check ends

                count_split++;
                token = strtok(NULL, delimiters);
            }
            //tokenizing ends until the 7th column
        }
        //check of line counter ends

        check_line_counter++;
    }
    //while loop ends
    fclose(fp_annovar);
    return 0;
}
