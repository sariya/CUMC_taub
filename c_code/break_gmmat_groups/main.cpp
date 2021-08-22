#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "helper_functions.hpp"
//Date 08 20 2021
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
    char *line = NULL; //read line
    size_t len = 0;    //this is read while reading file

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
    FILE *fp_annovar = fopen(annov_file, "r"); //add check if file doesn't exists

    if (fp_annovar == NULL)
    {
        fprintf(stderr, "Input file for annotation isn't correct. Exiting\n");
        exit(EXIT_FAILURE);
    }
    if (argc != 5)
    {
        printf("we do not have sufficient arguments\n");
        print_usage();
        exit(0);
    }

    const char *delimiters = "\t "; //use this for spliting lines
    int check_line_counter = 0;     //we do not want to read headers of the annovar output

    std::map<std::string, std::list<snp_store> > gene_snp_map; //use this to store gene
    const char *delimit_gene = ";,";                           //use this for the seventh column  SMIM12,DLGAP3 or  SMIM1;CCDC27

    int count_split; //use this for counting split. as columns are fixed
    while ((getline(&line, &len, fp_annovar)) != -1)
    {
        if (check_line_counter >= 1)
        {
            count_split = 1;                   //use this for counting split. as columns are fixed
            std::string temp_snp_name = "chr"; //we will use this for later uses : chr19:12456:A:T
            std::string func_ann = "";         // use this with SNP struct

            remove_trailingspaces(line); //remove trailing new line character

            char *token = NULL;
            //int count_split = 1; //use this for counting split. as columns are fixed
            token = strtok(line, delimiters);
            snp_store snp_temp_info = snp_store(); //make a new struct of SNP
            while (token != NULL)
            {
                if (count_split == 1)
                {
                    //temp_snp_name = temp_snp_name + token + ":"; // make chr1:
                    snp_temp_info.chromosome = token;
                }
                if (count_split == 2)
                {
                    //position
                    //temp_snp_name = temp_snp_name + token + ":"; // make chr1:123:
                    snp_temp_info.position = token;
                }
                if (count_split == 4)
                {
                    if (strlen(token) == 1)
                    {
                        //reference allele
                        //temp_snp_name = temp_snp_name + token + ":"; // make chr1:123:A
                        snp_temp_info.ref_allele = token[0];
                    }
                    else
                    {
                        break;
                    }
                }
                if (count_split == 5)
                {
                    //alternate allele
                    if (strlen(token) == 1)
                    {
                        //temp_snp_name = temp_snp_name + token; // make chr1:123:A:T
                        snp_temp_info.alt_allele = token[0];
                    }
                    else
                    {
                        break;
                    }
                }
                if (count_split == 6)
                {
                    //func information
                    func_ann = token;
                }
                if (count_split == 7)
                {
                    char *gene_token; //use this to parse gene information  //gene_store gene_info;
                    gene_token = strtok(token, delimit_gene); //use gene token as key in the map storage
                    while (gene_token)
                    {
                        //keep parsing that information until nothing is left                         //snp_temp_info.weight = '1';
                        gene_snp_map[gene_token].push_back(snp_temp_info);
                        gene_token = strtok(NULL, delimit_gene);
                    }
                }
                if (count_split > 7)
                {
                    break; //no more splitting after column 7
                }
                count_split++;
                token = strtok(NULL, delimiters);
            }
        }
        check_line_counter++;
    }
    if (line)
    {
        free(line);
    }
    fclose(fp_annovar);

    /**
     * we will print genes to each sep output file
    */

    std::map<std::string, std::list<snp_store> >::iterator it_map; //use it to iterate gene-snpp information
    std::cout << "Unique genes as " << gene_snp_map.size() << '\n';
    std::string new_filename = output_file;     //use this for each gene
    std::list<snp_store>::iterator it_snp_list; //use it later in the for loop

    for (it_map = gene_snp_map.begin(); it_map != gene_snp_map.end(); ++it_map)
    {
        //do not use endl  https://stackoverflow.com/a/17468264/2740831
        //this is key value pair
        new_filename = output_file;
        
        new_filename = new_filename + it_map->first + ".txt";
        //std::cout << "gene is " << it_map->first << " with length " << (it_map->second).size() << '\t' << "file name " << new_filename << '\n';
        if ((it_map->second).size() >= 2)
        {
            ///only if genes have two SNPs at least
            std::ofstream gene_output_stream(new_filename); //use this to output genes
            std::string gene_pos_alleles = ""; //make a string for each gene-snp-ref-alt allele

            //only if we have more than 2 SNPs in it
            for (it_snp_list = (it_map->second).begin(); it_snp_list != (it_map->second).end(); ++it_snp_list)
            {
                //std::cout << it_map->first << '\t' << it_snp_list->chromosome << '\t'  << it_snp_list->position << '\t' << it_snp_list->ref_allele << '\t' << it_snp_list->alt_allele << '\t' << it_snp_list->weight << '\n';
                gene_pos_alleles = it_map->first+'\t'+it_snp_list->chromosome+'\t'+it_snp_list->position+'\t'+it_snp_list->ref_allele+'\t'+it_snp_list->alt_allele+'\t'+'1'+'\n';
                gene_output_stream << gene_pos_alleles;
                gene_pos_alleles = "";
            }
            gene_output_stream.close(); //outfile.close();
        }
        //     //iterate over SNPs within the gene

    }
    //for loop ends of the iterator
    printf("Check current working directory for gene files with prefix of output provided. \n");
    return 0;
}