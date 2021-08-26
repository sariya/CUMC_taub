#ifndef FUNC_parse_HPP
#define FUNC_parse_HPP

// #define MAX_LEN 500
// #define MAX_LEN_GENE 200
#include <getopt.h>

void remove_trailingspaces(char *newline) ;
void print_usage();
void assign_variables(int count_argc, char **argv_array, const char *short_array, const struct option long_opt_temp[], char **temp_annovar_file, char **temp_out_file);

#endif


