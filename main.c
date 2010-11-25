/*  main.c
    
    See STABIL.c for an explanation.

    - Keshav Kini <kini@member.ams.org>, 2010-10-25
*/

#define ALLOC(X, Y) (X = malloc(Y * sizeof *(X)))

#include <stdlib.h>
#include <stdio.h>
#include "STABIL.h"

int main (int argc, char* argv[]) {
    unsigned long n, d;                                                         /* rank and dimension */
    unsigned long* matrix;                                                      /* for storing the matrix */
    unsigned long i, j, ij;                                                     /* counters */
    int result;                                                                 /* for checking the return value of STABIL() */
    FILE* f;                                                                    /* input file */
    
    if (argc < 2)
        f = stdin;
    else {
        f = fopen(argv[1], "r");
        if (!f) {
            fprintf(stderr, "HALT: File %s not found!\n", argv[1]);
            fprintf(stderr, "      Usage: `STABIL filename` (read from file) or `STABIL` (read from stdin)\n");
            return EXIT_BAD_INPUT;
        }
    }

    fscanf(f, "%lu%lu", &d, &n);                                                /* get rank and dimension from the first two lines of the input */
    
    ALLOC(matrix, n*n);
    for (ij = 0, i = 0; i < n; ++i)
        for (j = 0; j < n; ++j, ++ij)
            fscanf(f, "%lu", matrix + ij);                                      /* get matrix entries from input */
    fclose(f);

    if (n > 0xffffUL) {
        fprintf(stderr, "HALT: Predicted integer overflow!\n");
        fprintf(stderr, "      Please do not use matrices larger than 65535x65535.\n");
        return EXIT_OVERFLOW;
    }

    result = STABIL(matrix, n, &d);

    switch(result) {
    case EXIT_BAD_INPUT:
        fprintf(stderr, "HALT: Malformed input data!\n");
        fprintf(stderr, "      Enter first a dimension d, then a rank n, then the entries of a matrix\n");
        fprintf(stderr, "      of dimension n x n whose entries are, as a set, {0, 1, ..., d-1} and whose\n");
        fprintf(stderr, "      diagonal and off-diagonal entries are disjoint.\n");
        break;
    case EXIT_ALLOC_ERROR:
        fprintf(stderr, "HALT: Could not allocate sufficient memory!\n");
        break;
    case EXIT_SUCCESS:
        printf("%lu \n%lu \n", d, n);
        for (ij = 0, i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j, ++ij)
                printf("%lu ", matrix[ij]);
            printf("\n");
        }
    }
    
    free(matrix);
    return result;
}
