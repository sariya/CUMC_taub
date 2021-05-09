#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

////
#include "parse_function.h"

//////////////////

int main(int argc, char *argv[])
{
    // gcc -Wpedantic -Wextra -Wall *.c -o ./parse_annot

    gene_store *start_gene = NULL;
    //
    //   FILE *fp = fopen("/mnt/mfs/hgrcgrid/shared/MHAS/data/GENOTYPED/ANALYSES/gene_basedanalysis/annotation/INFO40/CHR22_MHAS_INFO40_plink.hg38_multianno.txt", "r");

    FILE *fp = fopen("test_commas.txt", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    char *line = NULL;
    size_t len = 0;
    int line_number = 0;
    while ((getline(&line, &len, fp)) != -1)
    {
        remove_trailingspaces(line);
        char *token = NULL;
        char SNP[400];
        token = strtok(line, "\t");
        int count_split = 1;

        if (line_number > 1)
        {
            while (token != NULL && ++count_split <= 8)
            {

                if (count_split == 2)
                {
                    chr_concat(token, SNP); //make 21 as chr21:
                }
                if (count_split == 3)
                {
                    join_strings(token, ':', SNP); //make  chr21 as chr21:1234:
                }
                if (count_split == 5)
                {
                    join_strings(token, ':', SNP); //make   chr21:1234: as chr21:1234:A:
                }
                if (count_split == 6)
                {
                    join_strings(token, '\0', SNP); //make  chr21:1234:A as chr21:1234:A:T
                }

                if (count_split == 8)
                {

                    char *gene_token = NULL;
                    const char *delimiters = ",;";
                    //gene_token = strtok(token, ",;");
                    gene_token = strtok(token, delimiters);

                    while (gene_token)
                    {
                        printf("we have gene name as %s\n", gene_token);
                        if (search_gene(&start_gene, gene_token) == 0)
                        {
                            gene_add_node(&start_gene, gene_token, SNP);
                        }
                        else
                        {
                            add_SNP_to_exiting_gene(&start_gene, gene_token, SNP);
                        }
                        //gene_token = strtok(NULL, ",;");

                        gene_token = strtok(NULL, delimiters);
                    }
                }
                //check for column 8th - that is where all genes are stored.

                token = strtok(NULL, "\t");
            }
            SNP[0] = '\0';
            // using printf() in all tests for consistency
            //printf("line is %s\n", line);
        }
        // printf("we have line number as %d\n", line_number);
        line_number++;
    }
    fclose(fp);

    if (line)
    {
        free(line);
    }
    //    display_gene_list(start_gene);
    delete_linked_list_gene(&start_gene);

    return 0;
}
