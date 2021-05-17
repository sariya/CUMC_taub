
//Date 05 17 2021
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions_group.h"

// void copy_strings(char *temp_string, char *dest_string, size_t str_max_length)
// {

//     dest_string[0] = '\0';
//     if (temp_string == NULL)
//     {
//         printf("Null string in the function\n");
//         return;
//     }

//     size_t len_string = strlen(temp_string);
//     size_t itr_string = 0;

//     if (len_string > str_max_length)
//     {
//         fprintf(stderr, "We have error in string length %s\n", temp_string);
//         return;
//     }

//     while (itr_string < len_string && temp_string[itr_string] != '\0')
//     {
//         dest_string[itr_string] = temp_string[itr_string];
//         itr_string++;
//     }
//     dest_string[itr_string] = '\0';
// }
// //////////////////////////////////////
int main(int argc, char *argv[])
{

    // gcc -Wpedantic -Wextra -Wall main_parse_gtf.c -o gtf_parse
    //FILE *fp_gtf = fopen("test.gtf", "r");
    //  cat test.gtf | grep gene | grep -v exon | grep -v transcript

    //FILE *fp_gtf = fopen("error_snp.txt", "r"); // /mnt/mfs/ctcn/resources/GRCh38/v1/Homo_sapiens.GRCh38.93.filtered.gtf

    FILE *fp_gtf = fopen("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/ADSP20K/ANNOTATION/files_fixed/VEP_HIGH/VEP_ADSP20K_CHR22_high", "r"); // /mnt/mfs/ctcn/resources/GRCh38/v1/Homo_sapiens.GRCh38.93.filtered.gtf //input_gene_annot_pos.txt
    // FILE *fp_gtf = fopen("input_gene_annot.txt", "r"); // /mnt/mfs/ctcn/resources/GRCh38/v1/Homo_sapiens.GRCh38.93.filtered.gtf //input_gene_annot_pos.txt

    if (fp_gtf == NULL)
    {
        exit(EXIT_FAILURE);
    }

    gene_store *start_gene = NULL;

    char *line = NULL;
    size_t len = 0;

    const char *delimiters = "\t "; //use this for spliting lines
    while ((getline(&line, &len, fp_gtf)) != -1)
    {
        char *token = NULL;
        remove_trailingspaces(line); //remove trailing new line character

        int count_split = 1; //use this for counting split. as columns are fixed
        token = strtok(line, delimiters);
        char snp_name_parsed[MAX_LEN_SNP];
        char gene_name[MAX_LEN_GENE];

        while (token != NULL)
        {
            if (count_split == 1)
            {
                copy_strings(token, gene_name, MAX_LEN_GENE);
            }
            if (count_split == 2)
            {
                copy_strings(token, snp_name_parsed, MAX_LEN_SNP);

                char *snp_info = NULL;
                snp_info = strtok(token, ":"); //
                int count_snp_split = 1;
                while (snp_info)
                {
                    if (count_snp_split == 2)
                    {
                        // printf("SNP position is %s\n", snp_info);
                        // printf("returned we have gene as %s\n", gene_name);
                        // printf("returned SNP as %s\n", snp_name_parsed);                               //(unsigned int *)
                        // printf("SNP position casted %u\n", (unsigned int)strtoull(snp_info, NULL, 0)); //https://stackoverflow.com/a/13025496/2740831

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
    //printf("total number of nodes in gene are %x\n", get_length_gene_nodes(start_gene));
    //get_length_gene_nodes(start_gene);

    display_gene_list(start_gene);
    return 0;
}