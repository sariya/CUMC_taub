
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions_groups_gmmat.h"

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
///
char *mystrcat(char *dest, char *src)
{
    /** function to concatenate two strings. 
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
