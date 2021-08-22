#include <getopt.h>
#include <string>

//date 08 20 2021
typedef struct snp
{
    char ref_allele;
    char alt_allele;
    std::string position;
    std::string chromosome;
} snp_store;

void print_usage();
void assign_variables(int count_argc, char **argv_array,
                        const char *short_array, const struct option long_opt_temp[], 
                        char **temp_annovar_file, char **temp_out_file);
void remove_trailingspaces(char *newline); //remove 

