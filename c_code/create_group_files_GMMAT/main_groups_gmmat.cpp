#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <fstream>

#include "functions_groups_gmmat.hpp"

static struct option long_opts[] = {
    /* ... */
    {"annov_file", required_argument, NULL, 'a'}, //use this for annotation input file
    {"out_file", required_argument, NULL, 'o'},   //use this for output file
    {NULL, 0, NULL, 0}};
static char short_option[] = "a:o:h"; //use this for parsing // this could have  *short_option
// //////////////////////////////////////
int main(int argc, char *argv[])
{
    /**
     * //Date 06 15 2021
     * we will read file from the user and output it in a tab delimited manner
     * 
    */
    char *annov_file = NULL;  // arg for annovar file .
    char *output_file = NULL; // arg for output file .
    assign_variables(argc, argv, short_option, long_opts, &annov_file, &output_file);

    if (annov_file == NULL)
    {
        printf("We do not have annovar file. Exiting\n");
        exit(0);
    }

    if (output_file == NULL)
    {
        printf("We do not have output file. Exiting\n");
        exit(0);
    }

    std::ofstream gmmat_output_stream(output_file);

    FILE *fp_annovar = fopen(annov_file, "r"); //FILE *fp_annovar = fopen(argv[1], "r"); //FILE *fp_annovar = fopen("input_annotation.txt", "r");

    if (fp_annovar == NULL)
    {
        fprintf(stderr, "Input file for annotation isn't correct. Exiting\n");
        exit(EXIT_FAILURE);
    }
    char *line = NULL;
    size_t len = 0;
    int check_line_counter = 0;

    const char *delimiters = "\t ";    //use this for spliting lines
    const char *delimiters_gene = ","; //use this for spliting gene column

    while ((getline(&line, &len, fp_annovar)) != -1)
    {
        remove_trailingspaces(line); //remove trailing new line character
        std::string annot_string = "chr";
        std::string gene_name="";

        if (check_line_counter > 1)
        {
            char *token = NULL; // we use this for splitting and storing
            int count_split = 1; //use this for counting split. as columns are fixed

            token = strtok(line, delimiters);

            while (token != NULL && count_split <= 7)
            {
                if (count_split == 1)
                {
                    annot_string = annot_string + token + '\t'; //make chr1
                }
                //first column check ends

                if (count_split == 2)
                {
                    annot_string = annot_string + token + '\t'; // CHR1 123(position)
                }
                //second column check ends

                if (count_split == 4)
                {
                    if (strlen(token) == 1 && strcmp(token, ".") != 0)
                    {
                        //Ref allele
                        annot_string = annot_string + token + '\t';//mystrcat(string_annotation, "\t");    mystrcat(string_annotation, token);
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
                        //alt allele
                        annot_string = annot_string + token + '\t' + "1";
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
                    //iterate over the column with gene information
                    if (strcmp(token, ".") != 0)
                    {
                        char *token_gene = NULL;
                        token_gene = strtok(token, delimiters_gene);

                        while (token_gene != NULL)
                        {
                            gene_name = token_gene;
                            gene_name = gene_name + '\t';   
                            gmmat_output_stream << gene_name + annot_string << '\n';
                            gene_name=""; //assign no string to gene
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
    //while loop ends of file reading
    if (line)
    {
        free(line);
    }

    //while loop ends

    fclose(fp_annovar);
    gmmat_output_stream.close(); //outfile.close();
    return 0;
}
