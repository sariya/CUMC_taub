
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "parse_function.h"

// 05 09 2021

///////////////////////////////////////////
int get_length_gene_nodes(gene_store *headnode)
{
    /**
     * get total number of nodes in the gene linked list
    */
    int number_of_nodes = 0;
    if (headnode == NULL)
    {
        printf("Nothing to count in get_length_gene_nodes function\n");
        return number_of_nodes;
    }
    else
    {
        while ((headnode) != NULL)
        {
            ++number_of_nodes;
            headnode = (headnode)->next;
        }
    }

    return number_of_nodes;
}
/////////////////// get_length_gene_nodes function ends

void remove_trailingspaces(char *newline)
{
    /**
    * remove trailing new line character from the line 
    */
    size_t line_length = strlen(newline);

    if (line_length == 0)
    {
        printf("remove_trailingspaces issues with line trimming\n");
        return;
    }
    else
    {
        newline[line_length - 1] = '\0';
    }
}
///////////////remove_trailingspaces function ends

void print_annotations(gene_store *headnode, char *outfile)
{
    /**
     * Function to display linkedlist 
     * 
    */
    FILE *write_ptr = fopen(outfile, "a"); //use this pointer to print gene and SNP

    if (write_ptr == NULL)
    {
        printf("In print_annotations we have error in file opening for print\n");
        return;
    }

    if (headnode == NULL)
    {
        printf("Nothing to display in print_annotations function\n");
        return;
    }
    while ((headnode) != NULL)
    {
        //printf("%s\n", (headnode)->name);
        fprintf(write_ptr, "%s\n", (headnode)->name);

        snp_store *temp_address_to_print = headnode->first_snp;

        while (temp_address_to_print != NULL)
        {
            //printf("%s\n", temp_address_to_print->snp_name);
            fprintf(write_ptr, "%s\n", temp_address_to_print->snp_name);
            temp_address_to_print = temp_address_to_print->next_snp_node;
        }

        //printf("END\n\n");
        fprintf(write_ptr, "%s\n", "END\n");
        (headnode) = (headnode)->next;
    }
    fclose(write_ptr);
}
// //////////////print_annotations function ends

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
        size_t itr = 0;                       //iterate over the string
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
////////join_strings function ends
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
///////////chr_concat function ends

void display_gene_list(gene_store *headnode)
{
    /**
     * Function to display linkedlist of genes with SNPs within it
     * 
    */
    if (headnode == NULL)
    {
        printf("Nothing to display in display_list function\n");
        return;
    }
    while ((headnode) != NULL)
    {
        printf("we have gene name %s \n", (headnode)->name);
        snp_store *temp_address_to_print = headnode->first_snp;

        while (temp_address_to_print != NULL)
        {
            printf("we have so and so name SNP %s \n", temp_address_to_print->snp_name);
            temp_address_to_print = temp_address_to_print->next_snp_node;
        }
        (headnode) = (headnode)->next;
    }
}

//////////////display_gene_list function ends
snp_store *return_snp_node()
{
    /**
     * Return malloc node address
     */
    snp_store *tempaddress_snp = malloc(sizeof(snp_store));
    return (tempaddress_snp);
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

////////// gene_store fumnntions ends
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
        snp_store *next_snp_node, *current_snp_node;

        while (next_node)
        {
            current_snp_node = current->first_snp;
            //  printf("we are doing gene as %s\n", next_node->name);

            next_snp_node = current_snp_node;

            while (next_snp_node)
            {
                next_snp_node = current_snp_node->next_snp_node;
                // printf("we will delete %s\n", current_snp_node->snp_name);
                free(current_snp_node);
                current_snp_node = next_snp_node;
            }

            next_node = current->next; //store next node
            free(current);             //delete node
            current = next_node;       //
        }
        *headnode = NULL;
    }
}
////////////////// delete_linked_list_gene function ends
void gene_add_node(gene_store **headnode, char *temp_gene, char *temp_snp)
{
    /**
     * Add node to the linked list
     */
    gene_store *temp_node_pointer = *headnode;
    gene_store *temp_node = return_node();

    snp_store *temp_snp_node = return_snp_node();

    temp_snp_node->next_snp_node = NULL;
    strcpy(temp_snp_node->snp_name, temp_snp);

    strcpy(temp_node->name, temp_gene);
    temp_node->next = NULL;

    if (temp_node_pointer == NULL)
    {
        temp_node_pointer = temp_node;
        temp_node_pointer->first_snp = temp_snp_node;
        *headnode = temp_node_pointer;
    }
    else
    {
        while ((temp_node_pointer->next) != NULL)
        {
            temp_node_pointer = temp_node_pointer->next;
        }

        temp_node->first_snp = temp_snp_node;
        temp_node_pointer->next = temp_node;
    }
}
////////////////gene_add_node function ends

void add_SNP_to_exiting_gene(gene_store **headnode, char *temp_gene_name, char *snp_name_temp)
{
    /**
     * If there is node of gene. then simply add SNP to the linked list
    */

    if (*headnode == NULL)
    {
        printf("No linked list exists. Exiting search_gene\n");
    }

    else if (strlen(temp_gene_name) == 0)
    {
        printf("in function add_SNP_to_exiting_gene gene length is 0\n");
    }
    else if (strlen(snp_name_temp) == 0)
    {
        printf("in function add_SNP_to_exiting_gene SNP length is 0\n");
    }
    else
    {
        gene_store *temp_node_pointer = *headnode;
        snp_store *temp_snp_node = return_snp_node();

        temp_snp_node->next_snp_node = NULL;
        strcpy(temp_snp_node->snp_name, snp_name_temp);

        while ((temp_node_pointer) != NULL)
        {
            if (strcmp(temp_gene_name, temp_node_pointer->name) == 0)
            {
                snp_store *store_to_swap = temp_node_pointer->first_snp;
                temp_node_pointer->first_snp = temp_snp_node;
                temp_node_pointer->first_snp->next_snp_node = store_to_swap;
            }
            temp_node_pointer = temp_node_pointer->next;
        }
        //iteration over the linked list ends
    }
}
///////// add_SNP_to_exiting_gene Function ends

int search_gene(gene_store **headnode, char *temp_gene_name)
{
    int flag_found_gene = 0;
    //  printf("we have value as %s\n", temp_gene_name);
    if (*headnode == NULL)
    {
        printf("No linked list exists. Exiting search_gene\n");
    }

    else if (strlen(temp_gene_name) == 0)
    {
        printf("in function search_gene gene length is 0\n");
    }

    else
    {
        gene_store *temp_node_pointer = *headnode;
        while ((temp_node_pointer) != NULL)
        {
            if (strcmp(temp_gene_name, temp_node_pointer->name) == 0)
            {
                flag_found_gene = 1;
                break;
            }
            temp_node_pointer = temp_node_pointer->next;
        }
    }

    return flag_found_gene;
}
///////// search_gene Function ends