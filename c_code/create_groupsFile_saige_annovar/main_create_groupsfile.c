#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions_create_saige.h"
//Date 07 26 2021

int main(int argc, char *argv[])
{
    /**
     * 
     * // New output in seedNum_126820	chr1:32300_A/C	chr1:32301_A/C	chr1:32302_A/C
     * 
    */
    if (argc != 3)
    {
        printf("we have incorrect number of input params\n");
        return -1;
    }

    FILE *fp = fopen(argv[1], "r"); //add check if file doesn't exists

    if (fp == NULL)
    {
        fprintf(stderr, "Input file for annotation isn't correct. Exiting\n");
        exit(EXIT_FAILURE);
    }

    gene_store *start_gene = NULL; //store all genes in this

    char *line = NULL;
    size_t len = 0;

    const char *delimiters = "\t "; //use this for spliting lines
    int check_line_counter = 0;
    const char *delimit_gene = ";,"; //use this for the seventh column  SMIM12,DLGAP3 or  SMIM1;CCDC27

    while ((getline(&line, &len, fp)) != -1)
    {
        //printf("%s\n ", line);
        char *token = NULL;
        remove_trailingspaces(line); //remove trailing new line character

        int count_split = 1; //use this for counting split. as columns are fixed
        token = strtok(line, delimiters);
        char snp_name_parsed[MAX_LEN_SNP]; //chr1:33879_A/C
        snp_name_parsed[0] = '\0';
        //char gene_name[MAX_LEN_GENE];
        char *SNP_POS = NULL;
        if (check_line_counter >= 1)
        {
            while (token != NULL)
            {
                if (count_split == 1)
                {
                    mystrcat(snp_name_parsed, "chr"); //make chr
                    mystrcat(snp_name_parsed, token); // make chr22
                    mystrcat(snp_name_parsed, ":");   // make chr22:
                }
                if (count_split == 2)
                {
                    //position
                    SNP_POS = token;                  //use this to add linked list and perform juggling
                    mystrcat(snp_name_parsed, token); // make chr1:32300
                    mystrcat(snp_name_parsed, "_");
                }

                if (count_split == 4)
                {
                    //reference allele
                    if (strlen(token) != 1)
                    {
                        snp_name_parsed[0] = '\0';
                        break;
                    }

                    mystrcat(snp_name_parsed, token); // make chr1:32300_A
                    mystrcat(snp_name_parsed, "/");   // make chr1:32300_A/
                }

                if (count_split == 5)
                {
                    //alternate allele
                    if (strlen(token) != 1)
                    {
                        snp_name_parsed[0] = '\0';
                        break;
                    }
                    mystrcat(snp_name_parsed, token); // make chr1:32300_T/A
                }

                if (count_split == 7)
                {
                    //7th column is where we have gene information
                    //printf("geneREF have %s\n", token);
                    //printf("checking %s\n",snp_name_parsed); ////chr1:33879_

                    char *gene_token;
                    gene_token = strtok(token, delimit_gene);

                    while (gene_token)
                    {
                        //printf("we have gene %s %s\n", snp_name_parsed, gene_token);

                        if (search_gene(&start_gene, gene_token) == 0)
                        {
                            fprintf(stderr, "No gene is not added %s\n", gene_token);
                            gene_add_node(&start_gene, gene_token, snp_name_parsed, (unsigned int)strtoull(SNP_POS, NULL, 0));
                            //printf("we have added SNP %s to gene %s. New list \n", snp_name_parsed, gene_token);
                        }
                        else
                        {
                            add_SNP_to_exiting_gene(&start_gene, gene_token, snp_name_parsed, (unsigned int)strtoull(SNP_POS, NULL, 0));
                            //  printf("we have added SNP %s to gene %s. Existing list \n", snp_name_parsed, gene_token);
                        }

                        gene_token = strtok(NULL, delimit_gene);
                    }
                }
                if (count_split > 7)
                {
                    //no more splitting after column 7
                    break;
                }
                count_split++;
                token = strtok(NULL, delimiters);
            }
        }

        check_line_counter++;
    }
    fclose(fp);

    // New output in seedNum_126820	chr1:32300_A/C	chr1:32301_A/C	chr1:32302_A/C
    printf("total number of nodes in gene are %u\n", get_length_gene_nodes(start_gene));

    //display_gene_list(start_gene);
    get_SNP_counts_within_gene(start_gene); //this will get counts of SNPs within each gene
    printf("all done. now printing......\n");
    print_group_files(start_gene, argv[2]); //print onto the file

    delete_linked_list_gene(&start_gene); //delete gene linked list
    printf("all done. exiting after deleting SNPs and genes list......\n");
    return 0;
}