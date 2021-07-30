#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <map>
#include <iostream>
#include <algorithm>


#include <string>
#include "functions_regenie.hpp"

// g++ -Wall -Wextra -O3   functions_regenie.cpp  create_regenie_files_main.cpp 
// compile with full optimization https://stackoverflow.com/a/17468320/2740831

int main(int argc, char *argv[])
{
    char *line = NULL;
    size_t len = 0;

    FILE *fp = fopen("small_input", "r"); //add check if file doesn't exists

    if (fp == NULL)
    {
        fprintf(stderr, "Input file for annotation isn't correct. Exiting\n");
        exit(EXIT_FAILURE);
    }

    const char *delimiters = "\t "; //use this for spliting lines
    int check_line_counter = 0;

    std::map<std::string, std::list<snp_store> > gene_snp_map;

    const char *delimit_gene = ";,"; //use this for the seventh column  SMIM12,DLGAP3 or  SMIM1;CCDC27

    while ((getline(&line, &len, fp)) != -1)
    {

        if (check_line_counter >= 1)
        {
            remove_trailingspaces(line); //remove trailing new line character

            char *token = NULL;

            int count_split = 1; //use this for counting split. as columns are fixed

            std::string temp_snp_name = "chr"; //we will use this for later uses : chr19:12456:A:T
            std::string func_ann = "";         // use this with SNP struct

            token = strtok(line, delimiters);

            while (token != NULL)
            {
                if (count_split == 1)
                {
                    temp_snp_name = temp_snp_name + token + ":"; // make chr1:
                }
                if (count_split == 2)
                {
                    //position
                    temp_snp_name = temp_snp_name + token + ":"; // make chr1:123:
                }
                if (count_split == 4)
                {
                    //reference allele
                    temp_snp_name = temp_snp_name + token + ":"; // make chr1:123:A
                }
                if (count_split == 5)
                {
                    //alternate allele
                    temp_snp_name = temp_snp_name + token; // make chr1:123:A:T
                }
                if (count_split == 6)
                {
                    //func information
                    func_ann = token; 
                }
                if (count_split == 7)
                {
                    char *gene_token; //use this to parse gene information  //gene_store gene_info;
                    
                    snp_store snp_temp_info = snp_store(); //make a new struct of SNP

                    snp_temp_info.func_annotation = func_ann; //store annotation
                    snp_temp_info.snp_name = temp_snp_name;   //store name

                    gene_token = strtok(token, delimit_gene); //use gene token as key in the map storage

                    while (gene_token)
                    {
                        //keep parsing that information until nothing is left
                        gene_snp_map[gene_token].push_back(snp_temp_info);
                        gene_token = strtok(NULL, delimit_gene);
                    }

                    //gene_info.gene_name = gene_token;                     //gene_names_list.push_back(gene_info);
                }
                if (count_split > 7)
                {
                    //no more splitting after column 7
                    break;
                }

                count_split++;
                token = strtok(NULL, delimiters);
            }
            token = strtok(line, delimiters); // spluit endlessly
        }

        check_line_counter++;
    }
    fclose(fp);

    std::map<std::string, std::list<snp_store> >::iterator it_map;

    for (it_map = gene_snp_map.begin(); it_map != gene_snp_map.end(); ++it_map)
    {
        ///do not use endl  https://stackoverflow.com/a/17468264/2740831

        //this is key value pair
        std::cout << "gene is " << it_map->first << " with length " << (it_map->second).size()  << '\n';
        //std::list<snp_store>::iterator it_snp_list;

        //iterate over SNPs within the gene
        // for (it_snp_list = (it_map->second).begin(); it_snp_list != (it_map->second).end(); ++it_snp_list)
        // {
        //     std::cout << "SNP value is " << it_snp_list->snp_name << "\t" << it_snp_list->func_annotation << std::endl;
        // }
    }

    return 0;
}

    //     std::list<Student>::iterator it;
    // std::list<gene_store>::iterator it_gene_list;
