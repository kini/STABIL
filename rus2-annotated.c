/* An implementation of the Weisfeiler-Leman graph stabilization procedure 
    see http://arxiv.org/abs/1002.1921 and references therein for more info.
    Complete sources can be found at http://www.ntu.edu.sg/home/dima/software.htm

    Copyright (C) 1989-2010 Luitpold Babel, Dmitrii V. Pasechnik, and others 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
    Authors can be contacted at
    e-mail: Luitpold Babel  <luitpold.babel@unibw.de>
    e-mail: Dmitrii Pasechnik <dimpase@gmail.com>                          */

#include <stdlib.h>
#include <stdio.h>
#define MAXN 200
#define MAXD MAXN*MAXN
#define MEMLENGTH 3*MAXN*MAXN
#define mlong int                                                               /* to fix an MSDOS artefact */
struct triple {
    int col;
    int val;
    struct triple *ptr;
};
struct edge {
    int row;
    int col;
    struct edge *ptr;
};

int main(int narg, char *arg[10])
{
    struct edge *color[MAXD];                                                   /* list of linked lists of edges indexed by color */
    int graph[MAXN][MAXN];                                                      /* colored adjacency matrix of the graph */
    int rank;                                                                   /* number of colors (rank of the configuration) */
    int vert;                                                                   /* number of vertices (dimension of the configuration) */
    int i, j;
    long time, start_time, end_time;
    FILE *f;

    printf("--- STABIL2 implementation of the Weisfeiler-Leman algorithm ---\n");

    f = fopen("input1", "r");                                                       /* char */
    if (f == NULL) {
        printf("HALT: File not found.\n");
        exit(0);
    }

    fscanf(f, "%d%d", &rank, &vert);                                            /* first two lines of the file are the rank
                                                                                   and dimension of the configuration */
    for (i = 0; i < vert; i++)
        for (j = 0; j < vert; j++)
            fscanf(f, "%d", &graph[i][j]);                                      /* subsequent lines are the adjacency matrix */
    i = edgepack(graph, rank, vert, color);
    if (i == 0) {
        printf("HALT: Data could not be interpreted.\n");
        exit(0);
    }
    printf("Number of colors (rank): ");

    start_time = clock() / 1000;
    stabil(&rank, vert, graph, color);
    end_time = clock() / 1000 - start_time;

    printf("\b\b\b\b\b\b%6d", rank);
    for (i = 0; i < rank; i++)
        color[i]->row = -1;
    j = 0;
    for (i = 0; i < vert; ++i) {
        if (color[graph[i][i]]->row < 0)
            color[graph[i][i]]->row = j++;
    }
    printf("\n\n number of cells: %6d", j);
    for (i = 0; i < rank; ++i)
        if (color[i]->row < 0)
            color[i]->row = j++;
    printf("\n\n adjacency matrix of the cellular algebra:\n\n");
    for (i = 0; i < vert; ++i) {
        for (j = 0; j < vert; ++j) {
            graph[i][j] = color[graph[i][j]]->row;
            if (graph[i][j] < 10)
                printf("  %d ", graph[i][j]);
            else {
                if (graph[i][j] < 100)
                    printf(" %d ", graph[i][j]);
                else
                    printf("%d ", graph[i][j]);
            }
        }
        printf("\n");
    }
    /* printf("\n\n%ld msec \n\n",end_time); */
}

/* check for bad input, populate the linked lists of edges */
int edgepack(int graph[MAXN][MAXN], int rank, int vert, struct edge *color[MAXD])
{
    static struct edge space[MAXD];
    int k, i, j;
    struct edge *free;

    for (i = 0; i < vert; ++i)
        for (j = 0; j < vert; ++j) {
            if (graph[i][j] < 0 || graph[i][j] > rank - 1)
                return (0);
        }
    for (free = space; free < space + rank; free++) /* clears the first d edges in space */
        free->row = 0;
    for (i = 0; i < vert; ++i) /* set .row to 1 for each color represented in the matrix diagonal */
        space[graph[i][i]].row += (space[graph[i][i]].row) ? 0 : 1;
    for (i = 0; i < vert; ++i) /* set .row to 2 for each color represented off the matrix diagonal */
        for (j = i + 1; j < vert; ++j) {
            if (space[graph[i][j]].row == 1)
                return (0); /* die if off-diagonal entry contains diagonal color */
            else
                space[graph[i][j]].row = 2;
            if (space[graph[j][i]].row == 1)
                return (0); /* die if off-diagonal entry contains diagonal color */
            else
                space[graph[j][i]].row = 2;
        }
    for (free = space; free < space + rank; free++)
        if (free->row == 0)
            return (0); /* die if any color 0 to d-1 is not found in the matrix */

    free = &space[0];
    for (k = 0; k < MAXD; k++) /* kill everything in struct edge** color */
        color[k] = NULL;
    for (i = 0; i < vert; i++) /* populate color[] */
        for (j = 0; j < vert; j++) {
            free->row = i;
            free->col = j;
            if (color[graph[i][j]] == NULL)
                free->ptr = NULL;
            else
                free->ptr = color[graph[i][j]];
            color[graph[i][j]] = free++;
        }
    return (1);
}

/* main algorithm */
int stabil(int* arank, int vert, int graph[MAXN][MAXN], struct edge* color[MAXD])
{
    int k, p, i, j, rank, klass, c, s, t, truth, overfl, q, oldq;
    int newrank, oldnrank, oldp;
    int *gamma;
    int memory[MEMLENGTH];
    struct edge *free, *w, *o, *oo;
    gamma = memory;
    rank = *arank;
    printf("%6d", rank);
    fflush(stdout);
    do {                                                                            /*! until new colors would not appear */
        truth = 0;                                                                  /*! new colors were not appear */
        newrank = rank;
        for (k = 0; k < rank; k++) {                                                /*! cycle on colors */
            overfl = 0;
            klass = 0;                                                              /*! number of new colors */
            p = 0;                                                                  /*! the begin of newgamma */
            *gamma = k;
            w = color[k];                                                       /* start with first edge of color k */
            o = w;                                                                  /*! the previous edge of color k */
            if (w->ptr == NULL)
                continue;                                                       /* no more edges of color k, move on to next k */

            do {                                                                /* for each edge (u,v) of color k */
                triangl(graph, w->row, w->col, gamma + p, rank, vert);          /* find all triangles with base (u,v), store their info */
                oldnrank = newrank;                                             /* in case we run out of memory */
                oldp = p;
                search(k, gamma, &p, &c, &klass, &s, &newrank, &truth, &q, &oldq);      /* do something */
                if (p >= (MEMLENGTH - (vert * 3 + 5)) || overfl == 1) {         /* forget everything gained from this edge (!?) because there's no space for further edges */
                    p = oldp;
                    newrank = oldnrank;
                    overfl = 1;
                    if (q == -1)
                        *(gamma + oldq + 4) = -1;
                    else if (oldq != -1) {
                        if (*(gamma + p + 3) != -1)
                            *(gamma + *(gamma + p + 3) + 2) = *(gamma + p + 2);
                        if (*(gamma + p + 2) != -1)
                            *(gamma + *(gamma + p + 2) + 3) = *(gamma + p + 3);
                    }
                }
                if (oldnrank != newrank) {
                    printf("\b\b\b\b\b\b%6d", newrank);
                    fflush(stdout);
                }
                if (c != k) {
                    o->ptr = w->ptr;
                    w->ptr = color[c];
                    color[c] = w;
                } else
                    o = w;
                w = o->ptr;
            } while (w != NULL);                                                /*! the last edge of color #k */
            
            if (overfl == 1)
                newrank++;
            for (i = rank; i < newrank; i++) {
                w = color[i];
                while (w != NULL) {
                    graph[w->row][w->col] = i;
                    w = w->ptr;
                }
            }
            rank = newrank;
        }                                                                       /*! next color */
        if (truth == 0)
            break;
        *arank = rank;
    } while (1);
}

/* search all paths of length 2 between vertices i and j, classify by color */
/* triangl(graph = graph, i = w->row, j = w->col, newgamma = gamma + p, rank = rank, vert = vert); */
int triangl(int graph[MAXN][MAXN], int i, int j, mlong* newgamma, int rank, int vert)
{
    static struct triple *lines[MAXD]; /* an array of linked lists, one for each possible color */
    struct triple *w;
    int s, t, p, numval, q;
    struct triple cnst[MAXN], *freemem;
    mlong *nnn;
    numval = 0;                                                                 /*! the number of nonzero const */
    freemem = cnst;
    for (p = 0; p < rank; p++)
        lines[p] = NULL;
    for (p = 0; p < vert; p++) {
        s = graph[i][p];    /* get the colors of the path from i to j through p */
        t = graph[p][j];
        pack(&lines[s], &freemem, t);
    }
    nnn = newgamma + 5;
    for (s = 0; s < rank; s++) {
        w = lines[s];
        while (w != NULL) {
            numval++;
            *(nnn++) = s;
            *(nnn++) = w->col;
            *(nnn++) = w->val;
            w = w->ptr;
        }
    }
    *(newgamma + 1) = numval;
}

/* write nonzero predicted structure coefficients into linked lists indexed by color(u,w) */
/* pack(line = lines+s, free = &freemem, t = t) */
int pack(struct triple **line, struct triple **free, int t) /* line is a pointer to a pointer, free is a pointer to an array */
{
    struct triple *w, *o;
    if (*line == NULL) {                                                        /* t is the first column in the line s */
        (*free)->col = t;
        (*free)->val = 1;
        (*free)->ptr = NULL;
        *line = *free;
        (*free)++;
    } else {
        if ((*line)->col > t) {                                                 /* the first column in the line > t */
            (*free)->col = t;
            (*free)->val = 1;
            (*free)->ptr = *line;
            *line = *free;
            (*free)++;
        } else {                                                                /* the first column in the line <=t  */
            w = *line;
            while (w != NULL) {
                if (w->col == t) {
                    w->val++;
                    break;
                }
                if (w->col > t) {
                    (*free)->col = t;
                    (*free)->val = 1;
                    (*free)->ptr = w;
                    o->ptr = *free;
                    (*free)++;
                    break;
                }
                o = w;
                w = w->ptr;
            }
            if (w == NULL) {
                o->ptr = *free;
                (*free)->col = t;
                (*free)->val = 1;
                (*free)->ptr = NULL;
                (*free)++;
            }
        }
    }
}

/* search(k = k, gamma = gamma, ap = &p, ac = &c, aklass = &klass, as = &s, anewrank = &newrank, atruth = &truth, aq = &q, aoldq = &oldq); */
int search(int k, mlong* gamma, int* ap, int* ac, int* aklass, int* as, int* anewrank, int* atruth, int* aq, int* aoldq)
{
    int q, oldp, oldq, nexte, dl, t, i;
    int c, p, klass, newrank, truth;
    p = *ap;
    klass = *aklass; /* number of new colors (?) */
    truth = *atruth; /* a flag */
    newrank = *anewrank; /* d_ */
    oldp = p;                                                                   /*! the begin of newgamma */
    q = 0;                                                                      /*! the begin of searching in gamma */
    if (klass) { /* color k has yielded divisions (?) */
        while (*(gamma + q + 1) != *(gamma + p + 1)) {
            if (*(gamma + q + 1) > *(gamma + p + 1)) {
                if (*(gamma + q + 2) == -1) {                                   /*!   prev==-1  */
                    *(gamma + q + 2) = p;
                    *(gamma + p + 2) = -1;
                    *(gamma + p) = newrank;
                    c = newrank;
                    newrank++;
                    klass++;
                    *(gamma + p + 3) = q;
                    *(gamma + p + 4) = -1;
                    truth = 1;
                    p = oldp + *(gamma + oldp + 1) * 3 + 5;
                    goto IR1;
                } else if (*(gamma + p + 1) > *(gamma + *(gamma + q + 2) + 1)) {    /*! between */
                    *(gamma + p + 2) = *(gamma + q + 2);
                    *(gamma + p + 3) = q;
                    *(gamma + q + 2) = p;
                    /* neu */ *(gamma + *(gamma + p + 2) + 3) = p;
                    *(gamma + p + 4) = -1;
                    *(gamma + p) = newrank;
                    c = newrank;
                    newrank++;
                    klass++;
                    truth = 1;
                    p = p + *(gamma + p + 1) * 3 + 5;
                    goto IR1;
                } else
                    q = *(gamma + q + 2);
            } else if (*(gamma + q + 3) == -1) {                                    /*! next==-1  */
                *(gamma + q + 3) = p;
                *(gamma + p + 3) = -1;
                *(gamma + p) = newrank;
                c = newrank;
                newrank++;
                klass++;
                *(gamma + p + 2) = q;
                *(gamma + p + 4) = -1;
                truth = 1;
                p = oldp + *(gamma + oldp + 1) * 3 + 5;
                goto IR1;
            } else if (*(gamma + p + 1) < *(gamma + *(gamma + q + 3) + 1)) {    /*! between */
                *(gamma + p + 3) = *(gamma + q + 3);
                *(gamma + p + 2) = q;
                *(gamma + q + 3) = p;
                /* neu */ *(gamma + *(gamma + p + 3) + 2) = p;
                *(gamma + p + 4) = -1;
                *(gamma + p) = newrank;
                c = newrank;
                newrank++;
                klass++;
                truth = 1;
                p = p + *(gamma + p + 1) * 3 + 5;
                goto IR1;
            } else
                q = *(gamma + q + 3);
        }
        
        do {
            oldq = q;
            dl = *(gamma + q + 1);
            nexte = *(gamma + q + 4);
            q += 5;
            p += 5;
            for (t = 1; t < dl * 3; t++, q++, p++) {
                if (*(gamma + q) != *(gamma + p))
                    break;
            }
            if (t == dl * 3) {                                                  /*! old class  */
                c = *(gamma + oldq);                                            /*! colour */
                p = oldp;
                oldq = -1;
                break;
            }
            if (nexte == -1) {                                                  /*! create a new class */
                klass++;
                *(gamma + oldp) = newrank;
                c = newrank;
                *(gamma + oldp + 4) = -1;
                *(gamma + oldq + 4) = oldp;
                *(gamma + oldp + 2) = *(gamma + oldq + 2);
                *(gamma + oldp + 3) = *(gamma + oldq + 3);
                newrank++;
                truth = 1;
                p = oldp + *(gamma + oldp + 1) * 3 + 5;
                /* neu */ q = nexte;
                break;
            }
            q = nexte;
            p = oldp;
        } while (q != -1);
    } else {
        klass++;                                                                /*! the first edge of colour k */
        *(gamma + p + 4) = -1;
        *(gamma + p + 3) = -1;
        *(gamma + p + 2) = -1;
        *(gamma + p) = k;
        p += *(gamma + p + 1) * 3 + 5;
        c = k;
    }
  IR1:*ac = c;
    *ap = p;
    *anewrank = newrank;
    *atruth = truth;
    *aklass = klass;
    *aq = q;
    *aoldq = oldq;
}
