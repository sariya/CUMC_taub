
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "common_include_file.h"
#include "funcs_combinations.h"

void delete_linked_list(data_store **headnode)
{
    data_store *next_node; //use it for iteration
    data_store *current = *headnode; //store current - begin from the start

    if (*headnode == NULL)
    {
        printf("No linked list exists\n");
        return;
    }
    else
    {
        while (next_node)
        {
            next_node = current->next; //store next node
            free(current); //delete node
        }
        *headnode = NULL;
    }
}
//////////////////////////////////////////////////////
void remove_trailingspaces(char *newlines)
{
    size_t length = strlen(newlines);
    newlines[length - 1] = '\0';
}
///////////////////////////////////////////
void get_gene_name(char temp_line[])
{
    int space_count = 0; //store location of space
    size_t itr = 0;
    size_t len_line = strlen(temp_line);
    if (len_line == 0)
    {
        return;
    }

    while (itr < len_line)
    {
        if (temp_line[itr] == ' ' || temp_line[itr] == '\t')
        {
            space_count = itr;
        }
        itr++;
    }

    itr = 0;       //iterate again
    space_count++; //start from next position
    while (space_count < len_line)
    {
        temp_line[itr++] = temp_line[space_count++]; //set character
    }
    temp_line[itr] = '\0';
}

///////////////////////////////////////////
data_store *return_node()
{
    data_store *tempaddress = malloc(sizeof(data_store));
    return (tempaddress);
}
//////////////////////////////////

void add_node(data_store **headnode, char *temp)
{
    data_store *temp_node_pointer = *headnode;
    data_store *temp_node = return_node();

    strcpy(temp_node->name, temp);
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

int check_length(data_store *headnode)
{
    int length = 0;
    if (headnode == NULL)
    {
        printf("nothing to show \n");
        return -1;
    }
    else
    {
        while ((headnode->next) != NULL)
        {
            length++;
            (headnode) = (headnode)->next;
        }
        length++;
        return (length);
    }
}
///////////////////////////////////////////////////
char *gene_combination(char *st1, char *st2, char *temp_string_with_asterix, char glue)
{
    size_t dest_maxlen = MAX_LEN; //max len is from common include functions

    size_t left_len = strlen(st1);

    temp_string_with_asterix[0] = '\0';
    strncat(temp_string_with_asterix, st1, dest_maxlen);
    if (dest_maxlen > left_len)
    {
        temp_string_with_asterix[left_len] = glue;
        temp_string_with_asterix[left_len + 1] = '\0';
        strncat(temp_string_with_asterix + left_len + 1, st2, dest_maxlen - left_len - 1);
    }

    return (temp_string_with_asterix);
}

///////////////////////////////////////
void display_list(data_store *headnode)
{
    while ((headnode) != NULL)
    {
        printf("we have so and so name %s\n", (headnode)->name);
        (headnode) = (headnode)->next;
    }
}

////////////////////////////////////////////////
void create_combination(data_store *node_temp_gene, char *file_name)
{
    data_store *iterate = node_temp_gene;
    data_store *iterate2;
    char temp_string_with_asterix[MAX_LEN];

    FILE *write_ptr; //use this pointer to write combinations

    if (file_name == NULL)
    {
        printf("we have error with file name \n");
    }

    write_ptr = fopen(file_name, "a");

    if (check_length(node_temp_gene) == 1 || check_length(node_temp_gene) == 0)
    {
        printf("nothing to make combination of\n");
    }
    else
    {
        while (iterate != NULL)
        {
            iterate2 = iterate->next;
            while (iterate2 != NULL)
            {
                gene_combination(iterate->name, iterate2->name, temp_string_with_asterix, '*');

                fprintf(write_ptr, "%s\n", temp_string_with_asterix);

                iterate2 = iterate2->next;
            }
            /// while iterate2 ends
            iterate = iterate->next;
        } //while iterate ends
        fclose(write_ptr);
    }
}
////////////////////////////////////////////////////////////////////////////
