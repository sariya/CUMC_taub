#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "functions_parse_gtf.h"

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
