#include "spkmeans.h"

/*
gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans_main.c -o spkmeans_main -lm
run using:
make && ./spkmeans wam tests/test1.txt
*/

int parseFile(char *fileName)
/*
reads data from file with directory fileName.
return : 0 if run was successful, 1 otherwise.
*/
{
    if (fileName){}
    return 0;
}

int main(int argc, char *argv[])
{
    char goal, fileName;
    if(argc && argv[0] && goal && fileName){}
    /*
    goal = argv[GOAL_IDX];
    fileName = argv[FILENAME_IDX];
    res = parseFile(data);
    */
    printf("here\n");
    return 0;
}
