#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <map>
#include <iostream>
#include <algorithm>
#include <getopt.h>
#include <string>
#include "functions_regenie.hpp"

// g++ -Wall -Wextra -O3   functions_regenie.cpp  create_regenie_files_main.cpp
// compile with full optimization https://stackoverflow.com/a/17468320/2740831

void print_usage()
{
    printf("Usage: regenie_files_annovar -a annovar_file -o output_file \n");
}
/////

static struct option long_opts[] = {
    /* ... */
    {"ann_file", required_argument, NULL, 'a'}, //use this for annotation input file
    {"out_file", required_argument, NULL, 'o'}, //use this for output file
    {NULL, 0, NULL, 0}};

static char short_option[] = "a:o:h"; //use this for parsing // this could have  *short_option

int main(int argc, char *argv[])
{
    char *line = NULL;         //read line
    size_t len = 0;            //this is read while reading file
    char *annovar_file = NULL; // arg for annovar file .
    char *output_file = NULL;  // arg for output  file .
    int opt;

    while ((opt = getopt_long(argc, argv, short_option, long_opts, NULL)) != -1)
    {
        switch (opt)
        {
        case 'a':
            annovar_file = optarg;
            break;

        case 'o':
            output_file = optarg;
            break;
        case 'h':
            print_usage();
            exit(0);
            //break;
        }
        //switch case ends
    }
    //while ends

    if (argc != 5)
    {
        printf("we do not have sufficient arguments\n");
        print_usage();
        exit(0);
    }

    printf("Input annovar file is %s\n", annovar_file);
    printf("output file is %s\n", output_file);

    if (annovar_file == NULL)
    {
        printf("No annovar file. Exiting\n");
    }

    if (output_file == NULL)
    {
        printf("No output file. Exiting\n");
    }

    FILE *fp_annovar = fopen(annovar_file, "r"); //add check if file doesn't exists     // FILE *fp = fopen("small_input", "r"); //add check if file doesn't exists

    if (fp_annovar== NULL)
    {
        fprintf(stderr, "Input file for annotation isn't correct. Exiting\n");
        exit(EXIT_FAILURE);
    }

    const char *delimiters = "\t "; //use this for spliting lines
    int check_line_counter = 0; //we do not want to read headers of the annovar output

    std::map<std::string, std::list<snp_store> > gene_snp_map; //use this to store gene and SNPs within it

    const char *delimit_gene = ";,"; //use this for the seventh column  SMIM12,DLGAP3 or  SMIM1;CCDC27

    while ((getline(&line, &len, fp_annovar)) != -1)
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

                    while (gene_token){
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
    //while loop ends for reading annovar input file
    fclose(fp_annovar);
    
    std::map<std::string, std::list<snp_store> >::iterator it_map;

    for (it_map = gene_snp_map.begin(); it_map != gene_snp_map.end(); ++it_map)
    {
        ///do not use endl  https://stackoverflow.com/a/17468264/2740831

        //this is key value pair
        std::cout << "gene is " << it_map->first << " with length " << (it_map->second).size() << '\n';
        std::list<snp_store>::iterator it_snp_list;

        //iterate over SNPs within the gene
         for (it_snp_list = (it_map->second).begin(); it_snp_list != (it_map->second).end(); ++it_snp_list)
         {
             std::cout << "SNP value is " << it_snp_list->snp_name << "\t" << it_snp_list->func_annotation << std::endl;
         }
         ///for loop ends
    }

    return 0;
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

//     std::list<Student>::iterator it;
// std::list<gene_store>::iterator it_gene_list;
    // while (1)
    // {
    //     int option_index = 0;

    //     switch (opt)
    //     {
    //     // case 0: /* for long option */
    //     //     printf("option %s", long_options[option_index].name);
    //     //     if (optarg)
    //     //         printf(" with arg '%s'", optarg);
    //     //     printf("\n");
    //     //     break;

    //     case '?':
    //     default:
    //         printf("unknown return code:%o\n", opt);
    //     }
    // }
    ///