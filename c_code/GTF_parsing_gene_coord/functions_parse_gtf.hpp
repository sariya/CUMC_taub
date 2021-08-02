#ifndef FUNC_parse_H
#define FUNC_parse_H

#define MAX_LEN 500
#define MAX_LEN_GENE 200

using namespace std;
#include <string>
#include <getopt.h>

//added on 05 11 2021

struct gene
{
    std::string gene_id;        //ENS
    std::string gene_name;      // abca7
    std::string chromosome;     // 9
    std::string start_position; // 123
    std::string end_position; // 321 
    std::string GENE_SOURCE; // havana 
    std::string GENE_BIOTYPE; //anti sense, 3', overlallping
};

void remove_trailingspaces(char *newline);
void remove_double_quotes(char **temp_gene_info); //remove double quotes and parse gene info
int check_chromosome(char *chr_string);
void print_usage() ;
void assign_variables(int count_argc, char **argv_array, const char *short_array, const struct option long_opt_temp[], char **temp_gtf_file, char **temp_out_file);

#endif

//char *mystrcat(char *dest, char *src);