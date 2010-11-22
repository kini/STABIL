/*  STABIL.c
    
    This file contains an adaptation of the program "STABIL" written by Luitpold Babel and Dmitrii Pasechnik and
    described in the paper "Program Implementation of the Weisfeiler-Leman Algorithm", arXiv:1002.1921v1 . The paper is
    available at <http://arxiv.org/abs/1002.1921v1>, and the original code is made available by Dmitrii Pasechnik at
    <http://bit.ly/aVF0BH>.
    
    The adaptation consists of much rephrasing / renaming, fleshed-out code comments, and a Cython wrapper to
    facilitate usage of the program from Sage or other Python-based environments. Preprocessor constants for hard
    limiting of memory usage have also been done away with in favor of dynamic allocation. The storage requirement has
    been increased (though not beyond the previous complexity order) in the interest of type safety. For the
    theoretical details of the algorithm and a discussion of its history, purpose, and applications, please see the
    paper linked above, which is free-access.
    
    As the original implementation was released under the GNU General Public License v3, so too is this file and the
    associated helper files.
    
    - Keshav Kini <kini@member.ams.org>, 2010-10-05
*/

#define ALLOC(X, Y) (X = malloc(Y * sizeof *(X)))
#define CALLOC(X, Y) (X = calloc(Y, sizeof *(X)))

#ifdef DEBUG
    #include "STABIL-tests.h"
    #define DEBUG_PRINT_MATRIX() print_matrix(matrix, n, *d);
    #define DEBUG_PRINT_TRIPLES() print_triples(triples, theaders, hnav1, k, *d)
#else
    #define DEBUG_PRINT_MATRIX()
    #define DEBUG_PRINT_TRIPLES()
#endif

#include <stdlib.h>
#include "STABIL.h"

/*  STABIL()
    
    STABIL, an implementation of the Weisfeiler-Leman refinement algorithm for coherent configurations. The function
    accepts a coherent configuration in the form of an integer matrix (called a generalized graph) "matrix", which is
    equal to $\sum_{i=0}^{d-1} i A_i$, where $A_0,\dots,A_{d-1}$ are the elements of the coherent configuration as
    $n\times n$ matrices. The function arguments n and d are as above. Internal variable naming convention is from
    visualization of "matrix" as a colored digraph with cycles.
*/
int STABIL(unsigned long* matrix, unsigned long n, unsigned long* d)
{
    /* constants */
    const unsigned long triples_len = n*n;                                      /* naively this should be of size n^3, but a favorable tradeoff is made by allocating only n^2 blocks and doing the
                                                                                    computation more than once if n^2 blocks proves insufficient (i.e. too much disagreement in predicted structure
                                                                                    coefficients) */
    const unsigned long theaders_len = n*n;                                     /* naively, up to n^2 new colors could be discovered per old color per iteration, though this is actually impossible */
    
    /* (pointers to) large-scale structures */
    struct edge** color_classes;                                                /* a list of root pointers of linked lists of edge structs, indexed by color */
    struct edge* edges;                                                         /* memory for the linked lists in struct edge** color_classes to live in */
    char* color_type;                                                           /* to keep track of vertex / edge colors */
    struct triple2** uw_classes;                                                /* a list of root pointers of linked lists of triple2 structs, indexed by the phantom .uw */
    struct triple2* triple2s;                                                   /* memory for the linked lists in struct triple2** uw_classes to live in */
    struct triple* triples;                                                     /* memory for storing predicted structure coefficients p_{i,j}^k for fixed k */
    struct theader* theaders;                                                   /* memory for storing metadata about the predicted structure coefficient lists */
    
    /* values */
    char stable, overflow;                                                      /* flags */
    unsigned long a, b, i, j, k;                                                /* counters - a,b are vertices, i,j,k are colors corresponding to the indices of p_{i,j}^k, the structure coefficients */
    unsigned long ab, ua, av;                                                   /* 2-dimensional iterators */
    unsigned long d_;                                                           /* next available color for new color classes within the refinement loop on a fixed color k; current total number of
                                                                                    color classes */
    
    /* navigational pointers */
    struct edge* uv;                                                            /* pointer to current edge (u,v) */
    struct triple2* newt2;                                                      /* pointer to next free position in block "triple2s" */
    struct triple2* t2nav1;                                                     /* pointers for navigation in linked list "uw_classes[i]" */
    struct triple2* t2nav2;
    struct theader* hnav1;                                                      /* pointers for navigation in the block "theaders" */
    struct theader* hnav2;
    struct triple* tnav1;                                                       /* pointers for navigation in the block "triples" */
    struct triple* tnav2;
    
    /* copies */
    struct edge* uv_;
    
    /* check matrix for shenanigans */
    if (
        !CALLOC(color_type, *d)
    )
        return EXIT_ALLOC_ERROR;
    for (ab = 0, a = 0; a < n; ++a)
        for (b = 0; b < n; ++b, ++ab) {
            if (matrix[ab] < 0 || matrix[ab] >= *d)
                return EXIT_BAD_INPUT;                                          /* die if out-of-range color found */
            if (color_type[matrix[ab]] == 0)
                color_type[matrix[ab]] = 1 + (a == b);                          /* type 1 is off-diagonal (edge), type 2 is diagonal (vertex) */
            else if (color_type[matrix[ab]] != 1 + (a == b))
                return EXIT_BAD_INPUT;                                          /* die if diagonal/off-diagonal conflict found */
        }
    for (i = 0; i < *d; ++i)
        if (!color_type[i])
            return EXIT_BAD_INPUT;                                              /* die if any color in range is not found */
    free(color_type);                                                           /* we don't care about this anymore */
    DEBUG_PRINT_MATRIX();
    
    
/*  STEP 0
    populate struct edge** color_classes
*/
    if (
        !CALLOC(color_classes, n*n) ||                                          /* we may eventually have as many as n^2 colors */
        !ALLOC(edges, n*n)
    )
        return EXIT_ALLOC_ERROR;
    for (ab = 0, a = 0; a < n; ++a)
        for (b = 0; b < n; ++b, ++ab) {
            edges[ab].row = a;
            edges[ab].col = b;
            edges[ab].next = color_classes[matrix[ab]];                         /* prepend this edge to the correct linked list... */
            color_classes[matrix[ab]] = edges + ab;                             /* ... and rebase the linked list */
        }
    
    
/*  STEP 1
    do computation
*/
    if (
        !ALLOC(triples, triples_len) ||
        !ALLOC(theaders, theaders_len) ||
        !ALLOC(uw_classes, n*n) ||                                              /* we could eventually have up to n^2 colors in the final refined matrix */
        !ALLOC(triple2s, n)                                                     /* only n triangles can exist on a given (u,v) */
    )
        return EXIT_ALLOC_ERROR;
    do {
        stable = 1;
        d_ = *d;                                                                /* original colors run from 0 to d-1, so d is the first available new color */
        for (k = 0; k < *d; ++k) {                                              /* for each color k... */
            overflow = 0;
            hnav1 = theaders;                                                   /* theaders is overwritten for each new color k */
            tnav1 = triples;                                                    /* triples is overwritten for each new color k */
            uv = color_classes[k];                                              /* get the first edge of color k - this must always exist */
            uv_ = uv;                                                           /* old (u,v) */
            if (uv->next == NULL)
                continue;                                                       /* there's only one edge of color k, and thus one (true) prediction for each p_{i,j}^k, nothing to do */
            
            do {                                                                /* for each edge (u,v) with color k... */
                if (overflow) {                                                 /* no more space for new coeff lists; just set to the last color */
                    uv_->next = uv->next;                                       /* when (u,v) is the first edge of its color, overflow will not have occurred, so here uv_ is edge previous to uv */
                    uv->next = color_classes[d_ - 1];                           /* prepend (u,v) */
                    color_classes[d_ - 1] = uv;                                 /* rebase linked list */
                    uv = uv_->next;                                             /* move on to next (u,v) of original color k */
                    continue;
                }
                
                
/*  STEP 1(i)
    compute the predicted structure coefficients (in struct triple2 linked lists)
*/
                newt2 = triple2s;                                               /* triple2s is overwritten for each new edge (u,v) */
                for (i = 0; i < *d; ++i)
                    uw_classes[i] = NULL;                                       /* zero the pointers before building the linked lists of predicted coefficient data */
                for (
                    a = 0, ua = uv->row * n, av = uv->col;                      /* loop over possible vertices w (= a) with which to build a triangle on uv */
                    a < n;
                    ++a, ++ua, av += n
                ) {
                    i = matrix[ua];                                             /* color of (u,w) */
                    j = matrix[av];                                             /* color of (w,v) */
                    
                    if (!uw_classes[i]) {                                       /* begin a linked list for .uw = i; this will be stored in increasing order of j */
                        *newt2 = (struct triple2){j, 1, NULL};
                        uw_classes[i] = newt2++;
                    } else if (j < uw_classes[i]->wv) {                         /* we can prepend new triple2, as its j is minimal */
                        *newt2 = (struct triple2){j, 1, uw_classes[i]};
                        uw_classes[i] = newt2++;
                    } else if (j == uw_classes[i]->wv) {                        /* this triangle is already at the root - we can increment its count */
                        ++uw_classes[i]->coeff;
                    } else {                                                    /* j > uw_classes[i]->wv, so new triple2 will not be minimal - we must insert it appropriately */
                        t2nav1 = uw_classes[i];
                        while (1) {
                            if (!(t2nav2 = t2nav1->next)) {                     /* move forward; if cannot move forward, append a new triple */
                                *newt2 = (struct triple2){j, 1, NULL};
                                t2nav1->next = newt2++;
                            } else if (j == t2nav2->wv) {                       /* this triangle coloring was already found, so increment its count */
                                t2nav2->coeff++;
                            } else if (j < t2nav2->wv) {                        /* the next is too far, and the current is too soon, so insert new triple2 between them */
                                *newt2 = (struct triple2){j, 1, t2nav2};
                                t2nav1->next = newt2++;
                            } else {
                                t2nav1 = t2nav2;                                /* move along the linked list */
                                continue;
                            }
                            break;
                        }
                    }
                }
                
                
/*  STEP 1(ii)
    collect the nonzero predicted structure coefficients into the block "triples"
*/
                hnav1->data = tnav1;                                            /* position of hnav1's list */
                hnav1->len = 0;                                                 /* number of color-different triangles found on (u,v) so far */
                for (i = 0; i < *d; ++i) {                                      /* for each possible color i of (u,w), in increasing order ... */
                    t2nav1 = uw_classes[i];
                    while(t2nav1) {                                             /* ... record all color-distinct triangles with that color on their (u,w), in increasing order of j */
                        ++hnav1->len;
                        *tnav1++ = (struct triple){i, t2nav1->wv, t2nav1->coeff};
                        t2nav1 = t2nav1->next;
                    }                                                           /* note: this is a while rather than a do because this linked list may be empty */
                }
                
                
/*  STEP 1(iii)
    search the block MEMORY, here called "struct triple* triples", for structure coefficient prediction lists matching
    the one produced by (u,v); if such is found, set (u,v) to that color, else create a new color for (u,v)
*/
                if (hnav1 == theaders)                                          /* we're still on the first edge of color k, so we need to start building our structure of lists in block "triples" */
                    *hnav1 = (struct theader){hnav1->data, k, hnav1->len, NULL, NULL, NULL};
                else {                                                          /* coefficient lists already exist, so now we need to search through them */
/*
    The coefficient lists, one for each class of edges of color k, live in the memory block "triples", and are arranged
    in several linked lists, each containing coefficient lists of a certain length. These linked lists are in turn
    linked together into a doubly linked list, which is ordered 'left-to-right' in increasing order of the 'certain
    length' belonging to its elements. This structural information, and indeed the demarcations of the coefficient
    lists themselves, are maintained by a series of headers in the memory block "theaders".The following procedure will
    start at the beginning of "theaders" (where the chronologically first coefficient list's data was written) and
    search through this structure looking for a match with the current (u,v)'s coefficient list, or alternatively a
    place to insert it into the structure if it represents a new color class.
*/
                    hnav2 = theaders;                                           /* start at the header of the first edge of color k, wherever it may be in the structure */
                    while (hnav2->len > hnav1->len && hnav2->left)              /* move left until *hnav2 <= *hnav1 in length or reach left end */
                        hnav2 = hnav2->left;
                    while (hnav2->len < hnav1->len && hnav2->right)             /* move right until *hnav2 >= *hnav1 in length or reach right end */
                        hnav2 = hnav2->right;
                    if (hnav2->len != hnav1->len) {                             /* there was no linked list with the right size, so we need to create and correctly place a new one. */
                        if (hnav2->len < hnav1->len) {                          /* rightward search undershot, so we must have reached the right end, so append */
                            hnav1->left = hnav2;
                            hnav1->right = NULL;
                            hnav2->right = hnav1;
                        } else if (hnav2->len > hnav1->len) {                   /* rightward search overshot, so the correct position must be directly to the left */
                            if (hnav2->left) {                                  /* there exists stuff to the left, so insert */
                                hnav1->left = hnav2->left;
                                hnav1->right = hnav2;
                                hnav2->left->right = hnav1;
                                hnav2->left = hnav1;
                            } else {                                            /* this is the left end, so prepend */
                                hnav1->left = NULL;
                                hnav1->right = hnav2;
                                hnav2->left = hnav1;
                            }
                        }
                        hnav1->down = NULL;                                     /* we are creating a new linked list for this size of coeff lists, so initialize hnav1 as the end of the list */
                        hnav1->color = d_;                                      /* as this is a new coeff list, it must belong to a new color class */
                    } else {                                                    /* a linked list with the right size was found, so we need to check for matches more closely */
                        do {                                                    /* for each coeff list with the same size as (u,v)'s coeff list, compare element-by-element */
                            tnav1 = hnav1->data;
                            tnav2 = hnav2->data;
                            for (i = 0; i < hnav1->len; ++i, ++tnav1, ++tnav2) {
                                if (
                                    tnav1->uw != tnav2->uw ||
                                    tnav1->wv != tnav2->wv ||
                                    tnav1->coeff != tnav2->coeff
                                ) break;                                        /* found a mismatch, so get out of here */
                            }
                            if (i == hnav1->len) {                              /* didn't find any mismatches, so we have found the correct class */
                                hnav1->color = hnav2->color;
                                break;                                          /* no need to search further down the linked list, so get out of here */
                            } else
                                tnav1 = hnav1->data + hnav1->len;               /* reset tnav1 to point to the beginning of free memory in struct triple* triples */
                            if (!hnav2->down) {                                 /* there are no more classes in this linked list, so we need to create a new one and append it */
                                hnav1->down = NULL;
                                hnav1->color = d_;
                                hnav2->down = hnav1;
                                break;                                          /* we're done here */
                            }
                        } while ((hnav2 = hnav2->down));
                    }
                }
                
                if (hnav1->color != k) {                                        /* this won't happen on the first (u,v), so inside we can assume uv_ is uv's predecessor */
                    uv_->next = uv->next;                                       /* remove uv from color class k, which we are iterating through */
                    uv->next = color_classes[hnav1->color];                     /* prepend to uv to newly refined color's linked list */
                    color_classes[hnav1->color] = uv;                           /* rebase list */
                } else
                    uv_ = uv;                                                   /* uv has not been recolored, so move along color class k normally */
                uv = uv_->next;                                                 /* get next uv */
                
                if (hnav1->color == d_) {                                       /* if a new color was added, go to next free color and next free space for a header in block "theaders" */
                    ++d_;
                    ++hnav1;
                } else if (hnav1 == theaders)                                   /* if color k was "added" (i.e. if this was the first edge), d_ doesn't need to increase, but the header and data still
                                                                                    need to be preserved */
                    ++hnav1;
                else
                    tnav1 = hnav1->data;                                        /* tnav1 should by this point naturally be after hnav1's data block, and needs to be moved back */
                
                if (theaders + theaders_len - hnav1 < 1 || triples + triples_len - tnav1 < n) {
                    overflow = 1;                                               /* if no more space for coeff lists and headers, mark an overflow */
                    ++d_;                                                       /* increment d_ (again, maybe) in advance to provide a color class of edges (d-1) found after overflow */
                }
            } while (uv);                                                       /* move on to the next color when this color is exhausted */
            
            if (overflow && !color_classes[d_ - 1])
                --d_;                                                           /* undo provisional incrementation of d_ performed above if no new edges were added after we went into overflow mode */
            
            for (i = *d; i < d_; ++i) {                                         /* save color changes to matrix; no need to check i < *d as old color classes are only shrunk */
                uv = color_classes[i];
                do {
                    matrix[uv->row*n+uv->col] = i;
                } while ((uv = uv->next));
            }
            
            if (d_ > *d) {                                                      /* more color classes have been added during this iteration */
                *d = d_;
                stable = 0;
            }
            
            DEBUG_PRINT_TRIPLES();
        }                                                                       /* next color */
        
        DEBUG_PRINT_MATRIX();
    } while(!stable);                                                           /* continue until the process stabilizes */
    
    free(color_classes);
    free(edges);
    free(uw_classes);
    free(triple2s);
    free(theaders);
    free(triples);
    
    return EXIT_SUCCESS;
}
