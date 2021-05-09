#ifndef FUNC_parse_H
#define FUNC_parse_H

#define MAX_LEN 500
#define MAX_LEN_GENE 200

typedef struct snp
{
    char snp_name[MAX_LEN];
    struct snp *next;
} snp_store; 

typedef struct genes
{
    char name[MAX_LEN_GENE];
    struct genes *next;
    snp_store *first_snp;
} gene_store; 

snp_store *return_snp_node(); //create node 

gene_store *return_node() ; //create node 

void delete_linked_list_gene();  // delete nodes from the linked list
void display_gene_list(gene_store *headnode); //print linked list of the gene
//void gene_add_node( gene_store **headnode, char *temp);  //add gene 

void gene_add_node(gene_store **headnode, char *temp_gene, char *temp_snp) ;

int search_gene(gene_store **headnode, char *temp); 

void chr_concat(char *tempCHR, char *temp_SNP);

void join_strings(char *string_tojoin_toSNP, char glue_chr, char *snp_temp) ; ////

int check_comma(char *temp_column_genes) ;

void remove_trailingspaces(char *newline) ;
#endif