#ifndef FUNC_create_group_H
#define FUNC_create_group_H
using namespace std;
#include <string>
#include <getopt.h>
void print_usage();
void assign_variables(int count_argc, char **argv_array, const char *short_array, const struct option long_opt_temp[], char **temp_annovar_file, char **temp_out_file);

void remove_trailingspaces(char *newline);
//char *mystrcat(char *dest, char *src) ;
#endif

