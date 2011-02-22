r"""
Weisfeiler-Leman algorithm

The Weisfeiler-Leman algorithm for finding the coarsest coherent refinement of
a partition of `\Omega \times \Omega` for some finite set Omega. Applicable
to colored digraphs as a recoloration that refines the partition of arcs and
loops by colors. In some cases may return the automorphism partition of the
graph; in all cases returns something of which the automorphism partition of
the graph is a refinement.

AUTHORS:

- Keshav Kini (2010-12-16)

EXAMPLES:

Refine the product with itself of the cyclic graph on five points::

    >>> from sage.graphs.wlrefine import GraphWL
    >>> g = graphs.CycleGraph(5)
    >>> gg = g.cartesian_product(g)
    >>> gg2 = GraphWL(gg)
    >>> gg2
    Looped digraph on 25 vertices
    >>> set(gg2.edge_labels())
    set([1, 2, 3, 4, 5])

Refine a test matrix::

    >>> from sage.graphs.wlrefine import WL
    >>> m = Matrix(8, 8, [3, 1, 2, 1, 1, 2, 2, 2, 1, 0, 1, 2, 2, 1, 2, 2, 2,
    ...   1, 3, 1, 2, 2, 1, 2, 1, 2, 1, 0, 2, 2, 2, 1, 1, 2, 2, 2, 0, 1, 2,
    ...   1, 2, 1, 2, 2, 1, 3, 1, 2, 2, 2, 1, 2, 2, 1, 0, 1, 2, 2, 2, 1, 1,
    ...   2, 1, 3])
    >>> m
    [3 1 2 1 1 2 2 2]
    [1 0 1 2 2 1 2 2]
    [2 1 3 1 2 2 1 2]
    [1 2 1 0 2 2 2 1]
    [1 2 2 2 0 1 2 1]
    [2 1 2 2 1 3 1 2]
    [2 2 1 2 2 1 0 1]
    [2 2 2 1 1 2 1 3]
    >>> m2 = WL(m)
    >>> m2
    [1 2 3 2 2 3 5 3]
    [4 0 4 6 6 4 6 7]
    [3 2 1 2 5 3 2 3]
    [4 6 4 0 6 7 6 4]
    [4 6 7 6 0 4 6 4]
    [3 2 3 5 2 1 2 3]
    [7 6 4 6 6 4 0 4]
    [3 5 3 2 2 3 2 1]

"""
# XXX Using >>>/... doctest notation until trac #10458 is resolved. Please
# change to sage:/....: later if possible.

#   wlrefine.pyx
#
#   See STABIL.c for an explanation.
#
#   - Keshav Kini <kini@member.ams.org>, 2010-10-14
#

from libc.stdlib cimport malloc, free

cdef extern from "STABIL.c":
    int STABIL(unsigned long* matrix, unsigned long n, unsigned long* d)

def WL(mat, fix_colors=True, algorithm="STABIL"):
    r"""
    Perform Weisfeiler-Leman refinement on a matrix.

    INPUT:

    - ``mat`` -- a square Sage matrix whose set of entries is the set of
      consecutive integers from 0 to some d-1 and whose diagonal entries do
      not occur outside the diagonal

    - ``fix_colors`` -- (default: true) if true, we enforce the separation of
      colors on the diagonal from colors off the diagonal. Optional because
      it may be desirable not to do so (it is not a necessary condition for
      the formation of a cellular algebra).

    - ``algorithm`` -- (default: "STABIL") choose the algorithm to use.
      Currently supported algorithms are: "STABIL" (default)

    OUTPUT:

    - The Weisfeiler-Leman refinement of mat

    EXAMPLES:

    Refine a test matrix::
        
        >>> from sage.graphs.wlrefine import WL
        >>> m = Matrix(8, 8, [3, 1, 2, 1, 1, 2, 2, 2, 1, 0, 1, 2, 2, 1, 2,
        ...   2, 2, 1, 3, 1, 2, 2, 1, 2, 1, 2, 1, 0, 2, 2, 2, 1, 1, 2, 2, 2,
        ...   0, 1, 2, 1, 2, 1, 2, 2, 1, 3, 1, 2, 2, 2, 1, 2, 2, 1, 0, 1, 2,
        ...   2, 2, 1, 1, 2, 1, 3])
        >>> m
        [3 1 2 1 1 2 2 2]
        [1 0 1 2 2 1 2 2]
        [2 1 3 1 2 2 1 2]
        [1 2 1 0 2 2 2 1]
        [1 2 2 2 0 1 2 1]
        [2 1 2 2 1 3 1 2]
        [2 2 1 2 2 1 0 1]
        [2 2 2 1 1 2 1 3]
        >>> m2 = WL(m)
        >>> m2
        [1 2 3 2 2 3 5 3]
        [4 0 4 6 6 4 6 7]
        [3 2 1 2 5 3 2 3]
        [4 6 4 0 6 7 6 4]
        [4 6 7 6 0 4 6 4]
        [3 2 3 5 2 1 2 3]
        [7 6 4 6 6 4 0 4]
        [3 5 3 2 2 3 2 1]

    WL() optionally enforces diagonal/off-diagonal separation and mandatorily
    enforces reversibility of the matrix (the condition that the transpose of
    the matrix must be a recoloring of the matrix).

    Try a non-reversible matrix without diagonal/off-diagonal separation, with
    and without "color fixing"::
        
        >>> from sage.graphs.wlrefine import WL
        >>> m = Matrix(8, 8, [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1,
        ...   1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
        ...   0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1])
        >>> m
        [1 1 1 1 0 0 0 0]
        [1 1 1 1 0 0 0 0]
        [1 1 1 1 0 0 0 0]
        [1 1 1 1 0 0 0 0]
        [0 0 0 0 1 0 0 1]
        [0 0 0 0 1 0 0 1]
        [0 0 0 0 1 0 0 1]
        [0 0 0 0 1 0 0 1]
        >>> m2 = WL(m)
        >>> m2
        [ 1  3  3  3  2  7  7  2]
        [ 3  1  3  3  2  7  7  2]
        [ 3  3  1  3  2  7  7  2]
        [ 3  3  3  1  2  7  7  2]
        [ 8  8  8  8  6  4  4 11]
        [ 9  9  9  9  5  0 10  5]
        [ 9  9  9  9  5 10  0  5]
        [ 8  8  8  8 11  4  4  6]
        >>> m3 = WL(m, fix_colors=False)
        >>> m3
        [1 1 1 1 0 4 4 0]
        [1 1 1 1 0 4 4 0]
        [1 1 1 1 0 4 4 0]
        [1 1 1 1 0 4 4 0]
        [5 5 5 5 8 2 2 8]
        [6 6 6 6 3 7 7 3]
        [6 6 6 6 3 7 7 3]
        [5 5 5 5 8 2 2 8]

    NOTES:

    Uses a reimplementation of STABIL by Keshav Kini, based on original work by
    Luitpold Babel and Dmitrii Pasechnik as described in [Bab]_.

    REFERENCES:

    .. [Bab] L. Babel, I. V. Chuvaeva, M. Klin, D. V. Pasechnik. Program
       Implementation of the Weisfeiler-Leman Algorithm. arXiv preprint
       1002.1921v1.

    AUTHORS:

    - Keshav Kini (2010-12-10)

    """
    from sage.matrix.constructor import Matrix

    cdef unsigned long* c_matrix
    cdef unsigned long c_d
    if (mat.nrows() != mat.ncols()):
        raise ValueError, "Malformed input data! Please provide a square matrix."
    n = mat.nrows()
    mat = [x for y in mat for x in y]

    # fix colors, or not
    if (fix_colors):
        diag_map = dict([(y,x) for (x,y) in enumerate(set([mat[i*n+i] for i in range(n)]))])
        c_d = len(diag_map)
        offdiag_map = dict([(y,x+c_d) for (x,y) in enumerate(set([mat[i*n+j] for i in range(n) for j in range(n) if i != j]))])
        c_d += len(offdiag_map)
        for i in range(n):
            for j in range(n):
                if (i == j):
                    mat[i*n + j] = diag_map[mat[i*n + j]]
                else:
                    mat[i*n + j] = offdiag_map[mat[i*n + j]]
    else:
        c_d = max(mat) + 1

    # prepare C matrix for passing to STABIL()
    c_matrix = <unsigned long*>malloc(n*n*sizeof(unsigned long))
    for i in range(n*n):
        c_matrix[i] = mat[i]

    # run STABIL() and interpret the results
    try:
        result = STABIL(c_matrix, n, &c_d)
        if (result == 1):
            raise ValueError, "Malformed input data! Entries of matrix must consist, as a set, of consecutive integers from 0 to some d-1, and diagonal and non-diagonal entries must be disjoint."
        elif (result == 2):
            raise MemoryError, "Could not allocate enough memory!"
        elif (result == 3):
            raise OverflowError, "Predicted overflow! Please do not use matrices larger than 65535x65535."
        result = Matrix(n, n, [c_matrix[x] for x in range(n*n)])
    finally:
        free(c_matrix)

    return result

def GraphWL(g, digraph=True, ignore_weights=False):
    r"""
    Perform Weisfeiler-Leman refinement on a graph.

    INPUT:

    - ``g`` -- a Sage looped graph with edge weights representing colors

    - ``digraph`` -- (default: True) if false, return an undirected graph
      with colors merged by Sage's coercion (?), which may be desirable in
      some situations

    - ``ignore_weights`` -- (default: False) if true, ignore edge weights

    OUTPUT:

    - g after performing Weisfeiler-Leman refinement on its coloring

    EXAMPLES:

    Refine the product with itself of the cyclic graph on five points::
        
        >>> from sage.graphs.wlrefine import GraphWL
        >>> g = graphs.CycleGraph(5)
        >>> gg = g.cartesian_product(g)
        >>> gg2 = GraphWL(gg)
        >>> gg2
        Looped digraph on 25 vertices
        >>> set(gg2.edge_labels())
        set([1, 2, 3, 4, 5])

    NOTES:

    This is a wrapper for the WL() function.

    AUTHORS:

    - Keshav Kini (2010-12-10)

    """
    from sage.graphs.graph import Graph, DiGraph

    if (not g.weighted() or ignore_weights) and not g.has_loops():
        mat = WL(g.adjacency_matrix() + 2, fix_colors=False)                    # g.adjacency_matrix() is a 0-1 matrix so setting the diagonal to 2 satisfies the color conditions for the algorithm
    else:
        mat = WL(g.weighted_adjacency_matrix())

    if digraph:
        return DiGraph(mat, pos=g.get_pos(), format='weighted_adjacency_matrix')
    else:
        return Graph(mat, pos=g.get_pos(), format='weighted_adjacency_matrix')
