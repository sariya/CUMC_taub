#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions_parse_gtf.h"

//  gcc -Wpedantic -Wextra -Wall main_parse_gtf.c -o gtf_parse

//  make clean ; make ; ./gtf_parse  > genes
int check_chromosome(char *chr_string)
int main(int argc, char *argv[])
{
    // /mnt/mfs/ctcn/resources/GRCh37/v1/Homo_sapiens.GRCh37.75.gtf
    //FILE *fp_gtf = fopen("test.gtf", "r");
    //  cat test.gtf | grep gene | grep -v exon | grep -v transcript

    FILE *fp_gtf = fopen("/mnt/mfs/ctcn/resources/GRCh37/v1/Homo_sapiens.GRCh37.75.gtf", "r"); // /mnt/mfs/ctcn/resources/GRCh38/v1/Homo_sapiens.GRCh38.93.filtered.gtf

    if (fp_gtf == NULL)
        exit(EXIT_FAILURE);

    char *line = NULL;
    size_t len = 0;

    const char *delimiters = "\t"; //use this for spliting lines
    while ((getline(&line, &len, fp_gtf)) != -1)
    {
        if (line[0] != '#')
        {
            char *token = NULL;
            remove_trailingspaces(line); //remove trailing new line character

            token = strtok(line, delimiters);
            int count_split = 1;        //use this for counting split. as columns are fixed
            int found_gene = 0;         //if column 3rd is gene
            char gene_info_print[5000]; //use this to print CHR1 start end ENS_ID GENE_NAME origin and biotype
            gene_info_print[0] = '\0';

            while (token != NULL)
            {
                if (count_split == 1)
                {
                    //printf("we have chr as %s\n", token);

                    if (check_chromosome(token) == 0)
                    {
                        break;
                    }
                    else
                    {
                        mystrcat(gene_info_print, "chr");
                        mystrcat(gene_info_print, token);
                    }

                    //printf("what to print %s\n", gene_info_print);
                }
                if (count_split == 3 && (strcmp(token, "gene") == 0))
                {
                    found_gene = 1;
                }
                if (count_split == 4 && found_gene == 1)
                {
                    /**
                     * join with start position
                    */
                    //printf("we have start as %s\n", token);
                    mystrcat(gene_info_print, "\t");
                    mystrcat(gene_info_print, token);
                }

                if (count_split == 5 && found_gene == 1)
                {
                    // printf("we have end as %s\n", token);
                    /**
                     * join with end position
                    */
                    mystrcat(gene_info_print, "\t");
                    mystrcat(gene_info_print, token);
                }

                if (count_split == 3 && (strcmp(token, "gene") != 0))
                {
                    break;
                }
                if (count_split == 9 && found_gene == 1)
                {
                    char *gene_info = NULL;
                    gene_info = strtok(token, ";"); // gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";

                    while (gene_info)
                    {
                        if (gene_info[5] == 'i') // gene_id
                        {
                            remove_double_quotes(&gene_info); //gene_id "ENSG00000223972"; as ENSG00000223972
                            // printf("the name of ided is %s\n", gene_info);
                            mystrcat(gene_info_print, "\t");
                            mystrcat(gene_info_print, gene_info);
                        }

                        if (gene_info[6] == 'n') //gene_name
                        {
                            remove_double_quotes(&gene_info); //make gene_name "DDX11L1"  as DDX11L1
                            //   printf("the name of genename ided is %s\n", gene_info);
                            mystrcat(gene_info_print, "\t");
                            mystrcat(gene_info_print, gene_info);
                        }

                        if (gene_info[6] == 's') //  gene_
                        {
                            remove_double_quotes(&gene_info); //make gene_source "ensembl_havana" as  ensembl_havana
                            //printf("the name of source ided is %s\n", gene_info);
                            mystrcat(gene_info_print, "\t");
                            mystrcat(gene_info_print, gene_info);
                        }

                        if (gene_info[6] == 'b') //gene_biotype
                        {
                            remove_double_quotes(&gene_info); //make gene_biotype "pseudogene" as  pseudogene
                            // printf("the name of biotype ided is %s\n", gene_info);
                            mystrcat(gene_info_print, "\t");
                            mystrcat(gene_info_print, gene_info);
                        }

                        gene_info = strtok(NULL, ";");
                    }
                    printf("%s\n", gene_info_print); //print it on the screen
                    break;                           //move out of the while loop of line parsing.
                }
                ///parsing of column with gene info ends

                token = strtok(NULL, delimiters);
                count_split++;
            }
            //while loop ends for token split of delimiters
        }
        //if line doesn't start with an #

    } //while loop of file ends

    //while loop ends of file reading
    fclose(fp_gtf);

    if (line)
    {
        free(line);
    }

    return 0;
}
//main function ends