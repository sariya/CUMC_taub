#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include <list>

#include <map>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <string>
////
#include "parse_function.hpp"

int main(int argc, char *argv[])
{
    // gcc -Wpedantic -Wextra -Wall *.c -o ./parse_annot

    //gene_store *start_gene = NULL;
    //
    //FILE *fp = fopen("/mnt/mfs/hgrcgrid/shared/MHAS/data/GENOTYPED/ANALYSES/gene_basedanalysis/annotation/INFO40/CHR22_MHAS_INFO40_plink.hg38_multianno.txt", "r");
    FILE *fp = fopen(argv[1], "r");
    printf("we have %s %s %d\n", argv[1], argv[2], argc);

    if (argc != 3)
    {
        printf("we have incorrect number of input params\n");
        return -1;
    }

    //FILE *fp = fopen("test_commas.txt", "r");
    if (fp == NULL){
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    size_t len = 0;
    int line_number = 0;
    const char *gene_delimiters = ",;"; //use this for spliting gene names
    while ((getline(&line, &len, fp)) != -1)
    {
        remove_trailingspaces(line);
        char *token = NULL;
        //char SNP[400];
        token = strtok(line, "\t ");
        int count_split = 1;

        if (line_number >= 1) //if line counter is more than 0
        {
            while (token != NULL && ++count_split <= 8)
            {
                if (count_split == 2)
                {
                    //chr_concat(token, SNP); //make 21 as chr21:
                }
                if (count_split == 3)
                {
                    //join_strings(token, ':', SNP); //make  chr21 as chr21:1234:
                }
                if (count_split == 5)
                {
                    //join_strings(token, ':', SNP); //make   chr21:1234: as chr21:1234:A:
                }
                if (count_split == 6)
                {
                    //join_strings(token, '\0', SNP); //make  chr21:1234:A as chr21:1234:A:T
                }

                if (count_split == 8)
                {
                    char *gene_token = NULL;
                    gene_token = strtok(token, gene_delimiters); //we can have ; or , or none in the column for genes

                    while (gene_token)
                    {
                        // if (search_gene(&start_gene, gene_token) == 0)
                        // {
                        //     gene_add_node(&start_gene, gene_token, SNP);
                        // }
                        // else
                        // {
                        //     /**
                        //      * If gene exists in the linked list. simply add SNP name
                        //     */
                        //     add_SNP_to_exiting_gene(&start_gene, gene_token, SNP);
                        // }
                        gene_token = strtok(NULL, gene_delimiters);
                    }
                    //while loop ends of token for gene split
                }
                //check for column 8th - that is where all genes are stored.

                token = strtok(NULL, "\t ");
            }
            /// token and split column 8 check ends
            //SNP[0] = '\0';
        }

        //line_number check ends

        line_number++; //parse from second line
    }
    //while loop ends of file reading
    fclose(fp);

    if (line)
    {
        free(line);
    }

    printf("we'll print into file now\n");
    //display_gene_list(start_gene);

    //printf("The number of nodes are %d\n", get_length_gene_nodes(start_gene));
    //print_annotations(start_gene, argv[2]);
    printf("SNP-gene annotations have been printed in the output file provided\n");
    //delete_linked_list_gene(&start_gene);
    printf("Exiting code\n");
    return 0;
}
//////////////////////////////