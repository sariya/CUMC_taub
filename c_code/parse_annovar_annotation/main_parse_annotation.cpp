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

static struct option long_opts[] = {
    /* ... */
    {"annov_file", required_argument, NULL, 'a'}, //use this for annotation input file
    {"out_file", required_argument, NULL, 'o'},   //use this for output file
    {NULL, 0, NULL, 0}};
static char short_option[] = "a:o:h"; //use this for parsing // this could have  *short_option
int main(int argc, char *argv[])
{   
    char *annov_file = NULL;  // arg for annovar file .
    char *output_file = NULL; // arg for output file .
    assign_variables(argc, argv, short_option, long_opts, &annov_file, &output_file); //func to assign variables
    FILE *fp = fopen(annov_file, "r");

    if (fp == NULL)
    {
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    size_t len = 0;
    int line_number = 0;
    const char *gene_delimiters = ",;"; //use this for spliting gene names
    std::map<std::string, std::list<std::string>> gene_snp_store; //store list of snp-gene in a map

    while ((getline(&line, &len, fp)) != -1)
    {
        std::string gene_name = "";
        std::string snp_name = "chr"; //make chr12:3456:A:T
        remove_trailingspaces(line);
        char *token = NULL;
        
        token = strtok(line, "\t ");
        int count_split = 1;

        if (line_number >= 1) //if line counter is more than 0
        {
            while (token != NULL && ++count_split <= 8)
            {
                if (count_split == 2)
                {
                    snp_name = snp_name + token + ":"; //chr_concat(token, SNP); //make 21 as chr21:
                }
                if (count_split == 3)
                {
                    snp_name = snp_name + token + ":"; //join_strings(token, ':', SNP); //make  chr21 as chr21:1234:
                }
                if (count_split == 5)
                {   
                    if (strlen(token) == 1)
                    {
                        snp_name = snp_name + token + ":"; //join_strings(token, ':', SNP); //make   chr21:1234: as chr21:1234:A:
                    }
                    else
                    {
                        break;
                    }
                }
                if (count_split == 6)
                {
                    if (strlen(token) == 1)
                    {
                        snp_name = snp_name + token; //join_strings(token, '\0', SNP); //make  chr21:1234:A as chr21:1234:A:T
                    }
                    else
                    {
                        break;
                    }
                }

                if (count_split == 8)
                {
                    char *gene_token = NULL;
                    gene_token = strtok(token, gene_delimiters); //we can have ; or , or none in the column for genes

                    while (gene_token)
                    {
                        gene_name = gene_token;
                        if (gene_name.compare("NONE") != 0)
                        {
                            gene_snp_store[gene_name].push_back(snp_name);
                        }

                        gene_name = "";
                        gene_token = strtok(NULL, gene_delimiters);
                    }
                    //while loop ends of token for gene split
                }
                //check for column 8th - that is where all genes are stored.

                token = strtok(NULL, "\t ");
            }
            /// token and split column 8 check ends
            snp_name = "";
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

    std::ofstream annovar_output_stream(output_file); //use this to print to an output file
    
    std::map<std::string, std::list<std::string>>::iterator it_map; //use this to iterate over genes and SNPs within

    for (it_map = gene_snp_store.begin(); it_map != gene_snp_store.end(); ++it_map)
    {
        std::string print_annot = "";
        print_annot = it_map->first + '\n';

        std::list<std::string>::iterator it_snp_list;

        //iterate over SNPs within the gene
        for (it_snp_list = (it_map->second).begin(); it_snp_list != (it_map->second).end(); ++it_snp_list)
        {
            print_annot = print_annot + it_snp_list->c_str() + '\n';
        }
        //for loop of snp list ends
        print_annot=print_annot+"END\n"; //add final END to it

        annovar_output_stream <<print_annot <<'\n' ;
        
    }
    //for loop ends
    annovar_output_stream.close(); 

    printf("SNP-gene annotations have been printed in the output file provided\n");
    printf("Exiting code\n");
    return 0;
}
//////////////////////////////