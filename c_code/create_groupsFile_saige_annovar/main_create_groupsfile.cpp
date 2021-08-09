#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions_create_saige.hpp"
//Date 07 26 2021
#include <list>

#include <map>
#include <algorithm>

#include <iostream>
#include <fstream>

static struct option long_opts[] = {
    /* ... */
    {"annov_file", required_argument, NULL, 'a'}, //use this for annotation input file
    {"out_file", required_argument, NULL, 'o'},   //use this for output file
    {NULL, 0, NULL, 0}};
static char short_option[] = "a:o:h"; //use this for parsing // this could have  *short_option

int main(int argc, char *argv[])
{
    /**
     * 
     * // New output in seedNum_126820	chr1:32300_A/C	chr1:32301_A/C	chr1:32302_A/C
     * 
    */
    char *annov_file = NULL;  // arg for annovar file .
    char *output_file = NULL; // arg for output file .
    assign_variables(argc, argv, short_option, long_opts, &annov_file, &output_file);

    if (annov_file == NULL)
    {
        printf("We do not have annovar file. Exiting\n");
        print_usage();
        exit(0);
    }

    if (output_file == NULL)
    {
        printf("We do not have output file. Exiting\n");
        print_usage();
        exit(0);
    }

    FILE *fp = fopen(annov_file, "r"); //add check if file doesn't exists

    if (fp == NULL)
    {
        fprintf(stderr, "Input file for annotation isn't correct. Exiting\n");
        exit(EXIT_FAILURE);
    }

    //gene_store *start_gene = NULL; //store all genes in this

    char *line = NULL;
    size_t len = 0;

    const char *delimiters = "\t "; //use this for spliting lines
    int check_line_counter = 0;
    const char *delimit_gene = ";,";                            //use this for the seventh column  SMIM12,DLGAP3 or  SMIM1;CCDC27
    std::map<std::string, std::list<std::string>> gene_snp_map; //use this to store gene and SNPs within it

    while ((getline(&line, &len, fp)) != -1)
    {
        //printf("%s\n ", line);
        char *token = NULL;
        remove_trailingspaces(line); //remove trailing new line character

        int count_split = 1; //use this for counting split. as columns are fixed
        token = strtok(line, delimiters);
        //char snp_name_parsed[MAX_LEN_SNP]; //chr1:33879_A/C
        //snp_name_parsed[0] = '\0';
        //char gene_name[MAX_LEN_GENE];
        //char *SNP_POS = NULL;
        if (check_line_counter >= 1)
        {
            std::string snp_string = "chr";
            std::string gene_name = "";
            while (token != NULL)
            {
                if (count_split == 1)
                {
                    snp_string = snp_string + token + ":"; //make chr1:
                    //mystrcat(snp_name_parsed, "chr"); //make chr  //mystrcat(snp_name_parsed, token); // make chr22  //mystrcat(snp_name_parsed, ":");   // make chr22:
                }
                if (count_split == 2)
                {
                    //position
                    //SNP_POS = token;  //use this to add linked list and perform juggling
                    snp_string = snp_string + token + "_";
                    //mystrcat(snp_name_parsed, token); // make chr1:32300 //mystrcat(snp_name_parsed, "_");
                }

                if (count_split == 4)
                {
                    //reference allele
                    if (strlen(token) != 1)
                    {
                        //snp_name_parsed[0] = '\0';
                        break;
                    }

                    //mystrcat(snp_name_parsed, token); // make chr1:32300_A
                    snp_string = snp_string + token + "/"; //mystrcat(snp_name_parsed, "/"); // make chr1:32300_A/
                }

                if (count_split == 5)
                {
                    //alternate allele
                    if (strlen(token) != 1)
                    {
                        //snp_name_parsed[0] = '\0';
                        break;
                    }

                    snp_string = snp_string + token; //mystrcat(snp_name_parsed, token); // make chr1:32300_T/A
                }

                if (count_split == 7)
                {
                    //7th column is where we have gene information
                    //printf("geneREF have %s\n", token); //printf("checking %s\n",snp_name_parsed); ////chr1:33879_

                    char *gene_token;
                    gene_token = strtok(token, delimit_gene);

                    while (gene_token)
                    {
                        //printf("we have gene %s %s\n", snp_name_parsed, gene_token);

                        gene_name = gene_token;
                        if (gene_name.compare("NONE") != 0)
                        {
                            gene_snp_map[gene_token].push_back(snp_string);
                        }
                        //gene_add_node(&start_gene, gene_token, snp_name_parsed, (unsigned int)strtoull(SNP_POS, NULL, 0));
                        //printf("we have added SNP %s to gene %s. New list \n", snp_name_parsed, gene_token);
                        //}
                        //else
                        //{
                        //add_SNP_to_exiting_gene(&start_gene, gene_token, snp_name_parsed, (unsigned int)strtoull(SNP_POS, NULL, 0));
                        //  printf("we have added SNP %s to gene %s. Existing list \n", snp_name_parsed, gene_token);
                        //}

                        gene_token = strtok(NULL, delimit_gene);
                    }
                }
                if (count_split > 7)
                { //no more splitting after column 7
                    break;
                }
                count_split++;
                token = strtok(NULL, delimiters);
            }
        }
        check_line_counter++;
    }
    fclose(fp);
    //while loop ends of file reading

    if (line)
    {
        free(line);
    }
    // New output in seedNum_126820	chr1:32300_A/C	chr1:32301_A/C	chr1:32302_A/C  //printf("total number of nodes in gene are %u\n", get_length_gene_nodes(start_gene));

    //display_gene_list(start_gene);     //get_SNP_counts_within_gene(start_gene); //this will get counts of SNPs within each gene
    std::ofstream saige_output_stream(output_file);

    std::map<std::string, std::list<std::string>>::iterator it_map; //use this to iterate over genes and SNPs within

    for (it_map = gene_snp_map.begin(); it_map != gene_snp_map.end(); ++it_map)
    {
        std::string print_annot = "";
        print_annot = it_map->first + "\t";

        std::list<std::string>::iterator it_snp_list;

        //iterate over SNPs within the gene
        for (it_snp_list = (it_map->second).begin(); it_snp_list != (it_map->second).end(); ++it_snp_list)
        {
            print_annot = print_annot + it_snp_list->c_str() + "\t";
        }
        saige_output_stream << print_annot << '\n';
        ///for loop ends
    }
    saige_output_stream.close(); //outfile.close();
    /// for loop ends

    printf("all done. exiting after printing SNPs and genes list......\n");
    return 0;
}