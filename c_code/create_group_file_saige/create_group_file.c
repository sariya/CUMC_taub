#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions_group.h"
//Date 05 17 2021

void replace_semi_colon(char *temp_snp_string)
{
    /**
     * //make chr20:54159075:G:A to chr20:54159075_G/A
     * 
    */
    size_t count_colon = 0, index_itr = 0;

    while (*(temp_snp_string + index_itr))
    {  
        if (*(temp_snp_string + index_itr) == ':')
        {
            count_colon++; //make chr20:54159075:G:A to chr20:54159075_G/A
            if (count_colon == 2)
            {
                *(temp_snp_string + index_itr) = '_';
            }
            if (count_colon == 3)
            {
                *(temp_snp_string + index_itr) = '/';
            }
        }
        index_itr++;
    }
}
// //////////////////////////////////////
int main(int argc, char *argv[])
{
    /**
     * 
     * gcc -Wpedantic -Wextra -Wall main_parse_gtf.c -o gtf_parse
     * 
     * ./create_group  input_gene_annot.txt  checkoutput
     * // New output in seedNum_126820	chr1:32300_A/C	chr1:32301_A/C	chr1:32302_A/C
     * 
    */
    if (argc != 3)
    {
        printf("we have incorrect number of input params\n");
        return -1;
    }

    //add check if file doesn't exists
    //FILE *fp_gtf = fopen("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/ADSP20K/ANNOTATION/files_fixed/VEP_HIGH/VEP_ADSP20K_CHR22_high", "r");
    FILE *fp = fopen(argv[1], "r");

    if (fp == NULL)
    {
        fprintf(stderr, "Input file for annotation isn't correct. Exiting\n");
        exit(EXIT_FAILURE);
    }

    gene_store *start_gene = NULL; //store all genes in this

    char *line = NULL;
    size_t len = 0;

    const char *delimiters = "\t "; //use this for spliting lines
    while ((getline(&line, &len, fp)) != -1)
    {
        char *token = NULL;
        remove_trailingspaces(line); //remove trailing new line character

        int count_split = 1; //use this for counting split. as columns are fixed
        token = strtok(line, delimiters);
        char snp_name_parsed[MAX_LEN_SNP];
        char gene_name[MAX_LEN_GENE];

        while (token != NULL)
        {
            snp_name_parsed[0] = '\0';
            if (count_split == 1)
            {
                copy_strings(token, gene_name, MAX_LEN_GENE);
            }
            if (count_split == 2)
            {

                mystrcat(snp_name_parsed, "chr");     
                mystrcat(snp_name_parsed, token);
                replace_semi_colon(snp_name_parsed);
     
                char *snp_info = NULL;
                snp_info = strtok(token, ":"); //
                int count_snp_split = 1;
                while (snp_info)
                {
                    if (count_snp_split == 2)
                    {
                        if (search_gene(&start_gene, gene_name) == 0)
                        {
                            gene_add_node(&start_gene, gene_name, snp_name_parsed, (unsigned int)strtoull(snp_info, NULL, 0));
                            fprintf(stderr, "No gene is not added %s\n", gene_name);
                        }
                        else
                        {
                            add_SNP_to_exiting_gene(&start_gene, gene_name, snp_name_parsed, (unsigned int)strtoull(snp_info, NULL, 0));
                        }
                    }
                    snp_info = strtok(NULL, ":");
                    count_snp_split++;
                }
            }
            count_split++;
            token = strtok(NULL, delimiters);
        }
    }
    // New output in seedNum_126820	chr1:32300_A/C	chr1:32301_A/C	chr1:32302_A/C
    printf("total number of nodes in gene are %u\n", get_length_gene_nodes(start_gene));

    //display_gene_list(start_gene);

    print_group_files(start_gene, argv[2]); //print onto the file

    fclose(fp);
    delete_linked_list_gene(&start_gene); //delete gene linked list
    return 0;
}