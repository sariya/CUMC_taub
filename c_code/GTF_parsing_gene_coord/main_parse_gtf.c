#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions_parse_gtf.h"

void again_remove_double_quotes(char **temp_gene_info)
{
    if (*temp_gene_info == NULL)
    {
        printf("String recived is NULL\n");
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

//  gcc -Wpedantic -Wextra -Wall main_parse_gtf.c -o gtf_parse

int main(int argc, char *argv[])
{
    // /mnt/mfs/ctcn/resources/GRCh37/v1/Homo_sapiens.GRCh37.75.gtf
    //FILE *fp_gtf = fopen("test.gtf", "r");
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
            remove_trailingspaces(line);

            token = strtok(line, delimiters);
            int count_split = 1;
            int found_gene = 0; //if column 3rd is gene

            while (token != NULL)
            {
                if (count_split == 1)
                {
                    //printf("we have chr as %s\n", token);
                    ;
                }
                if (count_split == 3 && (strcmp(token, "gene") == 0))
                {
                    found_gene = 1;
                }
                if (count_split == 4 && found_gene == 1)
                {
                    printf("we have start as %s\n", token);
                }

                if (count_split == 5 && found_gene == 1)
                {
                    printf("we have end as %s\n", token);
                }

                if (count_split == 3 && (strcmp(token, "gene") != 0))
                {
                    //printf("garbage as %s\n", token);
                    ;
                }
                if (count_split == 9 && found_gene == 1)
                {
                    char *gene_info = NULL;
                    char *gene_id = NULL;
                    char *gene_source = NULL;
                    char *gene_biotype = NULL;

                    gene_info = strtok(token, ";"); // gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";
                    
                    while (gene_info)
                    {
                        if (gene_info[5] == 'i') // gene_id  
                        {
                            again_remove_double_quotes(&gene_info);
                            printf("the name of ided is %s\n", gene_info);
                        }

                        if (gene_info[6] == 'n') //gene_name 
                        {
                            again_remove_double_quotes(&gene_info);
                            printf("the name of genename ided is %s\n", gene_info);
                        }

                        if (gene_info[6] == 's') //  gene_
                        {
                            again_remove_double_quotes(&gene_info);
                            printf("the name of source ided is %s\n", gene_info);
                        }

                        if (gene_info[6] == 'b') //gene_biotype 
                        {
                            again_remove_double_quotes(&gene_info);
                            printf("the name of biotype ided is %s\n", gene_info);
                        }
                        //printf("Gene name is about to print %s\n", gene_name);
                        gene_info = strtok(NULL, ";");
                    }
                    break; //move out of the while loop of line parsing.
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