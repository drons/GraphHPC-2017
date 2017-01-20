#ifndef __GRAPH_HPC_DEFS_H
#define __GRAPH_HPC_DEFS_H

#include <cstdio>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <ctime>
#include <sys/types.h>
#include <iostream>
#include <queue>
#include <cstring>
#include <cfloat>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <utility>

#define AVG_VERTEX_DEGREE 16
#define FILENAME_LEN 256
#define eps 1e-6

typedef unsigned vertex_id_t;
typedef unsigned long long edge_id_t;

/* The graph data structure*/
typedef struct
{
    /***
     The minimal graph repesentation consists of:
     n        -- the number of vertices
     m        -- the number of edges
     endV     -- an array of size m that stores the
                 destination ID of an edge <src->dest>.
     rowsIndices -- an array of size (n + 1) that stores pointers to the endV array (CRS format).
                 The degree of vertex i is given by rowsIndices[i + 1] - rowsIndices[i], and the
                 edges out of i are stored in the contiguous block
                 endV[rowsIndices[i] .. rowsIndices[i + 1] - 1].
     Vertices are numbered from 0 in our internal representation.
     ***/
    vertex_id_t n;
    edge_id_t m;
    edge_id_t *rowsIndices;
    vertex_id_t *endV;
    
    /* other graph parameters */
    int scale; /* log2 of vertices number */
    int avg_vertex_degree; /* relation m / n */
    
    /* RMAT graph parameters */
    double a, b, c;
    bool permute_vertices;
    
    /* filename for storing of the graph */
    char filename[FILENAME_LEN];
} graph_t;

/* write graph to file */
void writeGraph(graph_t *G, char *filename);

/* read graph from file */
void readGraph(graph_t *G, char *filename);

/* print text graph in std out */
void printGraph(graph_t *G);

/* free graph */
void freeGraph(graph_t *G);

/* algorithm */
void run(graph_t *G, double *result);

#endif
