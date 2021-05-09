
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "parse_function.h"

// 05 09 2021

///////////////////////////////////////////
///////////////////////////////////////
//////////////////////////////

//////////////////
void remove_trailingspaces(char *newline)
{
    /**
    * remove trailing new line character from the line 
    */
    size_t line_length = strlen(newline);

    if (line_length== 0)
    {
        printf("remove_trailingspaces issues with line trimming\n");
        return;
    }
    else
    {
        newline[line_length - 1] = '\0';
    }
}
//////////////////////////////////////////////////
int check_comma(char *temp_column_genes)
{
    /**
     * if comma found in the string return 1, else 0
    */

    int found_comma = 0;
    size_t len_column_genes = strlen(temp_column_genes);

    if (len_column_genes == 0)
    {
        printf("we have issues in column of genes\n");
    }
    else
    {
        size_t itr = 0;
        while (itr < len_column_genes)
        {
            if (temp_column_genes[itr] == ',')
            {
                found_comma = 1;
                break;
            }
            itr++;
        }
        //while loop ends
    }
    //if length is more than 0

    return found_comma;
}
//////////////////
void join_strings(char *string_tojoin_toSNP, char glue_chr, char *snp_temp)
{
    /**
    * 
    */
    if (strlen(snp_temp) == 0 || strlen(string_tojoin_toSNP) == 0)
    {
        printf("we have issues in join_strings function \n");
        return;
    }

    else
    {
        size_t length_SNP = strlen(snp_temp); //SNP length
        size_t itr = 0; //iterate over the string
        size_t length_to_join = strlen(string_tojoin_toSNP);

        while (itr < length_to_join)
        {
            snp_temp[length_SNP + itr] = string_tojoin_toSNP[itr];
            itr++;
        }

        if (glue_chr != '\0')
        {
            snp_temp[length_SNP + itr] = glue_chr;
            snp_temp[length_SNP + itr + 1] = '\0';
        }
        if (glue_chr == '\0')
        {
            snp_temp[length_SNP + itr] = '\0';
        }
    }
}
////////
void chr_concat(char *tempCHR, char *temp_SNP)
{

    /**
     *   //tempCHR is 21
     * 
    */
    if (strlen(tempCHR) > 0)
    {
        temp_SNP[0] = 'c';
        temp_SNP[1] = 'h';
        temp_SNP[2] = 'r';

        size_t length_chr = strlen(tempCHR);
        size_t itr_len = 0;

        while (itr_len < length_chr)
        {
            temp_SNP[3 + itr_len] = tempCHR[itr_len];
            itr_len++;
        }

        temp_SNP[3 + itr_len] = ':';
        temp_SNP[3 + itr_len + 1] = '\0';
    }
    else
    {
        printf("we have issues in CHR number\n");
    }
}
//////////////////////////////
void display_gene_list(gene_store *headnode)
{
    /**
     * Function to display linkedlist 
     * 
    */
    if (headnode == NULL)
    {
        printf("Nothing to display in display_list function\n");
        return;
    }
    while ((headnode) != NULL)
    {
        printf("we have so and so name %s\n", (headnode)->name);
        (headnode) = (headnode)->next;
    }
}

gene_store *return_node()
{
    /**
     * 
     * Return malloc node address
     */
    gene_store *tempaddress = malloc(sizeof(gene_store));
    return (tempaddress);
}

/////////////////////////////////////////////////
void delete_linked_list_gene(gene_store **headnode)
{
    /**
     * 
     * This function is to delete linked list
    */

    gene_store *next_node;           //use it for iteration
    gene_store *current = *headnode; //store current - begin from the start

    if (*headnode == NULL)
    {
        printf("No linked list exists\n");
        return;
    }
    else
    {
        next_node = current;
        while (next_node)
        {
            next_node = current->next; //store next node
            free(current);             //delete node
            current = next_node;       //
        }
        *headnode = NULL;
    }
}
//////////////////////////////////////////////////////
void gene_add_node(gene_store **headnode, char *temp_gene, char *temp_snp)
{
    /**
     * 
     * Add node to the linked list
     */
    gene_store *temp_node_pointer = *headnode;
    gene_store *temp_node = return_node();

    strcpy(temp_node->name, temp_gene);
    temp_node->next = NULL;

    if (temp_node_pointer == NULL)
    {
        temp_node_pointer = temp_node;
        *headnode = temp_node_pointer;
    }
    else
    {
        while ((temp_node_pointer->next) != NULL)
        {
            temp_node_pointer = temp_node_pointer->next;
        }
        temp_node_pointer->next = temp_node;
    }
}
//////////////////////////////////

int search_gene(gene_store **headnode, char *temp_gene_name)
{
    int flag_found_gene = 0;
    printf("we have value as %s\n", temp_gene_name);
    if (*headnode == NULL)
    {
        printf("No linked list exists. Exiting search_gene\n");
    }
    
    else if(strlen(temp_gene_name) ==0){
        printf("in function search_gene gene length is 0\n");
    }

    else
    {
        gene_store *temp_node_pointer = *headnode;
        while ((temp_node_pointer) != NULL)
        {
            if(strcmp(temp_gene_name,temp_node_pointer->name) ==0){
                flag_found_gene=1;
                break;
            }
            temp_node_pointer = temp_node_pointer->next;
        }
    }

    return flag_found_gene;
}