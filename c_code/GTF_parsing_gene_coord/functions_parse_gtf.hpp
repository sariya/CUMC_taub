#ifndef FUNC_parse_H
#define FUNC_parse_H

#define MAX_LEN 500
#define MAX_LEN_GENE 200
 
//added on 05 11 2021

void remove_trailingspaces(char *newline) ;
void remove_double_quotes(char **temp_gene_info) ; //remove double quotes and parse gene info
char *mystrcat(char *dest, char *src) ;
int check_chromosome(char *chr_string);
#endif
