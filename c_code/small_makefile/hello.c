#include <stdio.h>
#include <stdlib.h>

//  gcc -Wpedantic -Wextra -Wall hello.c -o print 
int main(int argc, char *argv[])
{
    char array[] = "3 hello";
    printf("we will print it %s\n", array);
    return 0;
}