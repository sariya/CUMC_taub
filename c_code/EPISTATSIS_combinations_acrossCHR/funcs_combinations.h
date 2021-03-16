#ifndef FUNC_COMBINATIONS_H
#define FUNC_COMBINATIONS_H

#include "common_include_file.h"

/*
data structure of gene. Single linked list
*/
typedef struct genes
{
    char name[MAX_LEN];
    struct genes *next;
} data_store; 


data_store *return_node() ; //create node 

int check_length(data_store *headnode) ; // return length of linked list
char *gene_combination(char *st1, char *st2, char *temp_string_with_asterix, char glue) ;

void add_node(data_store **headnode, char *temp); // add node to the data linned list
void create_combination( data_store *node_temp_gene, data_store *start_second_list, char *, char*, char *file_name) ; //GENE1 and GENE2 returned GENE1*GENE2
void remove_trailingspaces(char *newlines) ; //remove last line chracter from file read
void get_gene_name(char temp_line[]) ; // extract gene name from read line from file
void display_list(data_store *headnode); //print linked list 
void delete_linked_list(); // delete nodes from the linked list

void concatenate_str(char *, char *) ; //use it for CHR making
#endif