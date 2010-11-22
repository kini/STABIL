/*  STABIL.h
    
    See STABIL.c for an explanation.
    
    - Keshav Kini <kini@member.ams.org>, 2010-10-14
*/

#ifndef STABIL_H
#define STABIL_H

#define EXIT_SUCCESS 0
#define EXIT_BAD_INPUT 1
#define EXIT_ALLOC_ERROR 2
#define EXIT_OVERFLOW 3

struct edge {                                                                   /* linkable edge struct */
    unsigned long row;
    unsigned long col;
    struct edge* next;
};
struct triple {                                                                 /* coeff indicates proto- p_{i,j}^k for current edge (u,v) during step 1 of STABIL */
    unsigned long uw;                                                           /* color of (u,w) */
    unsigned long wv;                                                           /* color of (w,v) */
    unsigned long coeff;                                                        /* p_{uw,wv}^k as predicted by current edge (u,v) of color k */
};
struct triple2 {                                                                /* linkable version of struct triple, pre-sorted by .uw */
    unsigned long wv;
    unsigned long coeff;
    struct triple2* next;
};
struct theader {                                                                /* a header for lists of struct triple elements, making them both doubly (left/right) and singly (down) linked */
    struct triple* data;                                                        /* where in the array of triples the list this header describes begins */
    unsigned long color;                                                        /* stores a new (refined) color for the edge that generated the coefficient list */
    unsigned long len;                                                          /* length of the subsequent list (measured in sizeof struct triple)*/
    struct theader* left;                                                       /* pointer to the header to the "left" (header of immediately shorter list) */
    struct theader* right;                                                      /* pointer to the header to the "right" (header of immediately longer list) */
    struct theader* down;                                                       /* pointer to the header "below" (header of next list of same length) */
};

int STABIL(unsigned long* matrix, unsigned long n, unsigned long* d);

#endif
