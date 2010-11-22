/*  STABIL-tests.h
    
    This file contains testing functions for STABIL() in STABIL.c .
    
    - Keshav Kini <kini@member.ams.org>, 2010-11-02
*/

#ifndef STABIL_TESTS_H
#define STABIL_TESTS_H

#include "STABIL.h"
#include <stdio.h>

int print_matrix(unsigned long* matrix, unsigned long n, unsigned long d)
{
    unsigned long i, j, ij;
    printf("Current configuration is of dimension %lu and degree %lu. Matrix:\n", d, n);
    for (ij = 0, i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j, ++ij)
            printf("%lu ", matrix[ij]);
        printf("\n");
    }
    printf("\n");
    
    return EXIT_SUCCESS;
}

int print_triples(struct triple* triples, struct theader* theaders, struct theader* hnav1, unsigned long k, unsigned long d) {
    struct theader* hnav;
    struct triple* tnav;
    unsigned long i;
    
    printf("Color %lu out of a current %lu colors has so far been refined as follows:", k, d);
    for (hnav = theaders; hnav < hnav1; hnav++) {
        printf("\n[%lu: color %lu, %lu triples, LDR = (", hnav - theaders, hnav->color, hnav->len);
        if(hnav->left) printf("%lu,", hnav->left - theaders); else printf("_,");
        if(hnav->down) printf("%lu,", hnav->down - theaders); else printf("_,");
        if(hnav->right) printf("%lu)]", hnav->right - theaders); else printf("_)]");
        for (i = 0, tnav = hnav->data; i < hnav->len; ++i, ++tnav)
            printf(" (%lu,%lu)x%lu;", tnav->uw, tnav->wv, tnav->coeff);
    }
    printf("\n\n");
    
    return EXIT_SUCCESS;
}

#endif
