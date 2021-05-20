
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "functions_group.h"

//Date 05 17 2021

/////////////
unsigned int nodes_snp_list(snp_store *temp_head_node)
{
    /**
     * return count of nodes in the snp list within a gene
    */

    unsigned int number_of_nodes = 0;
    if (temp_head_node == NULL)
    {
        printf("No SNP nodes in the gene nodes_snp_list function\n");
    }
    else
    {

        while ((temp_head_node) != NULL)
        {
            ++number_of_nodes;
            temp_head_node = (temp_head_node)->next_snp_node;
        }
    }
    return number_of_nodes;
}
//////Function ends

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
        char string_annotation[50000];
        string_annotation[0] = '\0';

        //printf("we have gene name %s \n", (headnode)->name);
        snp_store *temp_address_to_print = headnode->first_snp;

        unsigned int temp_count_snps_gene_list = nodes_snp_list(temp_address_to_print);
        mystrcat(string_annotation, (headnode)->name);

        if (temp_count_snps_gene_list == 1)
        {
            mystrcat(string_annotation, "\t"); //we alraedy added gene in the beginning
            mystrcat(string_annotation, temp_address_to_print->snp_name);
            mystrcat(string_annotation, "\n");
            printf("%s", string_annotation);
        }
        else
        {
            unsigned int itr_counter = 0;

            while (temp_address_to_print != NULL)
            {
                if (itr_counter < temp_count_snps_gene_list - 1)
                {
                    //printf("we have so and so name SNP %s \n", temp_address_to_print->snp_name);
                    mystrcat(string_annotation, "\t");
                    mystrcat(string_annotation, temp_address_to_print->snp_name);
                }
                else
                {
                    /////////////////
                    mystrcat(string_annotation, "\t"); //we alraedy added gene in the beginning
                    mystrcat(string_annotation, temp_address_to_print->snp_name);
                    mystrcat(string_annotation, "\n");
                }
                // if not the last node

                itr_counter++;
                temp_address_to_print = temp_address_to_print->next_snp_node;
            }
            //while iteration ends

            printf("%s", string_annotation);
        }
        //if more than one node in SNP list

        (headnode) = (headnode)->next;
    }
}
///////////////////////////////////////////

void print_group_files(gene_store *headnode, char *outfile)
{
    /**
     * Function to display linkedlist 
     * 
    */
  // FILE *write_ptr = fopen(outfile, "a+");
    FILE *write_ptr = fopen(outfile, "a"); //use this pointer to print gene and SNP

    if (write_ptr == NULL)
    {
        printf("In print_group_files we have error in file opening for print\n");
        return;
    }

    if (headnode == NULL)
    {
        printf("Nothing to display in print_group_files function\n");
        return;
    }

    while ((headnode) != NULL)
    {
        char string_annotation[700000];
        string_annotation[0] = '\0';

        //printf("we have gene name %s \n", (headnode)->name);
        snp_store *temp_address_to_print = headnode->first_snp;

        unsigned int temp_count_snps_gene_list = nodes_snp_list(temp_address_to_print);
        mystrcat(string_annotation, (headnode)->name);

        if (temp_count_snps_gene_list == 1)
        {
            mystrcat(string_annotation, "\t"); //we already added gene in the beginning
            mystrcat(string_annotation, temp_address_to_print->snp_name);
            mystrcat(string_annotation, "\n");
            ///printf("%s", string_annotation);
            fprintf(write_ptr, "%s", string_annotation);
        }
        else
        {
            unsigned int itr_counter = 0;

            while (temp_address_to_print != NULL)
            {
                if (itr_counter < temp_count_snps_gene_list - 1)
                {
                    //printf("SNP %s \n", temp_address_to_print->snp_name);
                    mystrcat(string_annotation, "\t");
                    mystrcat(string_annotation, temp_address_to_print->snp_name);
                }
                else
                {
                    /////////////////
                    mystrcat(string_annotation, "\t"); //we alraedy added gene in the beginning
                    mystrcat(string_annotation, temp_address_to_print->snp_name);
                    mystrcat(string_annotation, "\n");
                    //printf("we have so and so name SNP %s \n", temp_address_to_print->snp_name);
                }
                // if not the last node

                itr_counter++;
                temp_address_to_print = temp_address_to_print->next_snp_node;
            }
            //while iteration ends

            //printf("%s", string_annotation);
            fprintf(write_ptr, "%s", string_annotation);
        }
        //if more than one node in SNP list

        (headnode) = (headnode)->next;
    }
    // while iteration over the linked list ends
    fclose(write_ptr);
}
// //////////////print_group function ends

/////////////////////////////////////////////
unsigned int get_length_gene_nodes(gene_store *headnode)
{
    /**
     * get total number of nodes in the gene linked list
    */
    unsigned int number_of_nodes = 0;
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
/////////////////////////////////////////

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
snp_store *return_snp_node()
{
    /**
     * Return malloc node address
     */
    snp_store *tempaddress_snp = malloc(sizeof(snp_store));
    return (tempaddress_snp);
}
/////////////////////////
gene_store *return_node()
{
    /**
     * 
     * Return malloc node address
     */
    gene_store *tempaddress = malloc(sizeof(gene_store));
    return (tempaddress);
}
/////////////////////////

char *mystrcat(char *dest, char *src)
{
    /**
     * https://stackoverflow.com/a/21881314/2740831
    */
    while (*dest)
    {
        dest++;
    }
    while (*dest++ = *src++)
    {
        ;
    }

    return --dest;
}
////////////////////////////

void copy_strings(char *temp_string, char *dest_string, size_t str_max_length)
{
    /**
     * Do not use strcpy series https://devblogs.microsoft.com/oldnewthing/20050107-00/?p=36773
    */

    dest_string[0] = '\0';
    if (temp_string == NULL)
    {
        printf("Null string in the function\n");
        return;
    }

    size_t len_string = strlen(temp_string);
    size_t itr_string = 0;

    if (len_string > str_max_length)
    {
        fprintf(stderr, "We have error in string length %s\n", temp_string);
        return;
    }

    while (itr_string < len_string && temp_string[itr_string] != '\0')
    {
        dest_string[itr_string] = temp_string[itr_string];
        itr_string++;
    }
    dest_string[itr_string] = '\0';
}
//////////////////////////////////////

void gene_add_node(gene_store **headnode, char *temp_gene, char *temp_snp, unsigned int temp_snp_pos)
{
    /**
     * Add node to the linked list
     */
    gene_store *temp_node_pointer = *headnode;
    gene_store *temp_node = return_node();

    snp_store *temp_snp_node = return_snp_node();

    temp_snp_node->snp_pos = temp_snp_pos;
    temp_snp_node->next_snp_node = NULL;
    temp_snp_node->previous_snp_node = NULL;

    //strcpy(temp_snp_node->snp_name, temp_snp);
    copy_strings(temp_snp, temp_snp_node->snp_name, MAX_LEN_SNP);

    //strcpy(temp_node->name, temp_gene);
    copy_strings(temp_gene, temp_node->name, MAX_LEN_GENE);

    temp_node->next = NULL; //gene node next null

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
    // check for if head pointer is not NULL ends
}
////////////////gene_add_node function ends

void add_SNP_to_exiting_gene(gene_store **headnode, char *temp_gene_name, char *snp_name_temp, unsigned int temp_snp_pos)
{
    /**
     * If there is node of gene. then simply add SNP to the linked list
    */

    if (*headnode == NULL)
    {
        fprintf(stderr, "No linked list exists. Exiting search_gene %s\n", temp_gene_name);
    }

    else if (strlen(temp_gene_name) == 0)
    {
        fprintf(stderr, "in function add_SNP_to_exiting_gene gene length is 0\n");
    }
    else if (strlen(snp_name_temp) == 0)
    {
        fprintf(stderr, "in function add_SNP_to_exiting_gene SNP length is 0\n");
    }
    else
    {
        gene_store *temp_node_pointer = *headnode;
        snp_store *temp_snp_node = return_snp_node();

        temp_snp_node->previous_snp_node = NULL;
        temp_snp_node->next_snp_node = NULL;
        temp_snp_node->snp_pos = temp_snp_pos;

        copy_strings(snp_name_temp, temp_snp_node->snp_name, MAX_LEN_SNP);

        while ((temp_node_pointer) != NULL)
        {
            if (strcmp(temp_gene_name, temp_node_pointer->name) == 0)
            {
                snp_store *iterate_SNP_list = temp_node_pointer->first_snp;

                if (iterate_SNP_list->snp_pos >= temp_snp_pos)
                {
                    temp_snp_node->next_snp_node = temp_node_pointer->first_snp;
                    temp_snp_node->next_snp_node->previous_snp_node = temp_snp_node;
                    temp_node_pointer->first_snp = temp_snp_node;
                }
                else
                {

                    while (iterate_SNP_list->next_snp_node && iterate_SNP_list->snp_pos < temp_snp_pos)
                    {
                        iterate_SNP_list = iterate_SNP_list->next_snp_node;
                    }
                    temp_snp_node->next_snp_node = iterate_SNP_list->next_snp_node;

                    if (temp_snp_node->next_snp_node != NULL)
                    {
                        temp_snp_node->next_snp_node->previous_snp_node = temp_snp_node;
                    }
                    iterate_SNP_list->next_snp_node = temp_snp_node;
                    temp_snp_node->previous_snp_node = iterate_SNP_list;
                }

                //* /
            }
            temp_node_pointer = temp_node_pointer->next;
        }
        //iteration over the linked list ends
    }
}
///////// add_SNP_to_exiting_gene Function ends

int search_gene(gene_store **headnode, char *temp_gene_name)
{
    /**
     * function returns 1 if gene is already in the linked list
    */
    int flag_found_gene = 0; //return 1 if gene is found in the linked list
    
    if (*headnode == NULL)
    {
        fprintf(stderr, "No linked list exists. Exiting search_gene %s\n", temp_gene_name);
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
        // while loop ends of the iteration to find gene in the linked list
    }
    //else check ends to check if input gene length is not 0

    return flag_found_gene;
}
///////// search_gene Function ends

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
            next_snp_node = current_snp_node;

            while (next_snp_node)
            {
                next_snp_node = current_snp_node->next_snp_node;
                free(current_snp_node);
                current_snp_node = next_snp_node;
            }

            next_node = current->next; //store next node
            free(current);             //delete node
            current = next_node;       //
        }
        *headnode = NULL;
    }
    // else check ends when headnode is not NULL
}
////////////////// delete_linked_list_gene function ends

///////////