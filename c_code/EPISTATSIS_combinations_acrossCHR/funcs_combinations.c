
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "common_include_file.h"
#include "funcs_combinations.h"

void concatenate_str(char *name_str, char *chr_str)
{
    /**
     * 
     * Function to create CHR to CHR9, CHR10 and such.
    */
    int itr = 0;
    int second_itr = 0;
    if (strlen(name_str) == 0)
    {
        printf("Length 0 in concatenate_str function. Exiting\n");
        return;
    }
    if (strlen(chr_str) == 0)
    {
        printf("Length 0 in concatenate_str function, second argument. Exiting\n");
        return;
    }
    while (itr < strlen(name_str))
    {
        itr++;
    }

    while (second_itr < strlen(chr_str))
    {
        *(name_str + itr) = *(chr_str + second_itr);
        second_itr++;
        itr++;
    }
    *(name_str + CHR_LEN-1 ) = '\0'; //chr_LEN is defined in the include headers
}
/////////////////////////////////////////////////
void delete_linked_list(data_store **headnode)
{
    /**
     * 
     * This function is to delete linked list
    */

    data_store *next_node;           //use it for iteration
    data_store *current = *headnode; //store current - begin from the start

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
void remove_trailingspaces(char *newlines)
{
    /**
     * 
     * This function is to remove trailing space from the line read. file opened
    */

    size_t length = strlen(newlines);
    if (length == 0)
    {
        printf("we have zero length in remove_trailingspaces function\n");
    }

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
    /**
     * Function to display linkedlist 
     * 
    */
   if(headnode ==NULL){
       printf("Nothing to display\n");
       return;

   }
    while ((headnode) != NULL)
    {
        printf("we have so and so name %s\n", (headnode)->name);
        (headnode) = (headnode)->next;
    }
}

////////////////////////////////////////////////
void create_combination(data_store *node_temp_gene, data_store *start_second_list, char *chr_one, char *chr_second, char *file_name)
{
    data_store *iterate = node_temp_gene;
    data_store *iterate2 = start_second_list;
    char temp_string_with_asterix[MAX_LEN];

    // FILE *write_ptr; //use this pointer to write combinations

    if (file_name == NULL)
    {
        printf("we have error with file name \n");
    }

    if (iterate2 == NULL)
    {
        printf("Nothing in the second linked list\n");
        return;
    }
    //write_ptr = fopen(file_name, "a");

    if (check_length(node_temp_gene) == 1 || check_length(node_temp_gene) == 0)
    {
        printf("nothing to make combination of\n");
        return;
    }
    else
    {
        while (iterate != NULL)
        {
            iterate2 = start_second_list;
            while (iterate2 != NULL)
            {
                char gene_one[MAX_LEN]; //make CHR10_gene1
                char gene_two[MAX_LEN]; //make  CHR1_gene2
                gene_one[0] = '\0';
                gene_two[0] = '\0';

                gene_combination(chr_one, iterate->name, gene_one, '_');
                gene_combination(chr_second, iterate2->name, gene_two, '_');

                gene_combination(gene_one, gene_two, temp_string_with_asterix, '*'); //make asterikx

                printf("%s\n", temp_string_with_asterix);
                
                //fprintf(write_ptr, "%s\n", temp_string_with_asterix);

                iterate2 = iterate2->next;
            }
            /// while iterate2 ends
            iterate = iterate->next;
        } //while iterate ends
          // fclose(write_ptr);
    }
}
////////////////////////////////////////////////////////////////////////////
