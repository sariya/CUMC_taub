#ifndef FUNC_create_group_H
#define FUNC_create_group_H

//Date 05 17 2021

#define MAX_LEN_SNP 500U
#define MAX_LEN_GENE 200U

typedef struct snp
{
    char snp_name[MAX_LEN_SNP];
    unsigned int snp_pos;
    struct snp *next_snp_node, *previous_snp_node;
} snp_store; 

typedef struct genes
{
    char name[MAX_LEN_GENE];
    struct genes *next;
    snp_store *first_snp;
} gene_store; 

snp_store *return_snp_node();
gene_store *return_node();

void remove_trailingspaces(char *newline);
void gene_add_node(gene_store **headnode, char *temp_gene, char *temp_snp, unsigned int temp_bp) ; //
//void add_SNP_to_exiting_gene(gene_store **headnode, char *temp_gene_name, char *snp_name_temp) ;
void add_SNP_to_exiting_gene(gene_store **headnode, char *temp_gene_name, char *snp_name_temp, unsigned int temp_bp) ;
int search_gene(gene_store **headnode, char *temp_gene_name);

void copy_strings(char *temp_string, char *dest_string, size_t str_max_length);
unsigned int get_length_gene_nodes(gene_store *headnode) ;

void display_gene_list(gene_store *headnode) ;

char *mystrcat(char *dest, char *src) ;


void print_group_files(gene_store *headnode, char *outfile);
void delete_linked_list_gene(gene_store **headnode); 

#endif

