#include "defs.h"

using namespace std;

/* read graph from the file */
void readGraph(graph_t *G, char *filename)
{
    unsigned char align;
    FILE *f = fopen(filename, "rb");
    assert(f != NULL);
    assert(fread(&G->n, sizeof(vertex_id_t), 1, f) == 1);
    G->scale = log(G->n) / log(2);
    assert(fread(&G->m, sizeof(edge_id_t), 1, f) == 1);
    assert(fread(&align, sizeof(unsigned char), 1, f) == 1);
    G->rowsIndices = new edge_id_t[G->n + 1];
    assert(G->rowsIndices != NULL);
    assert(fread(G->rowsIndices, sizeof(edge_id_t), G->n + 1, f) == (G->n + 1));
    G->endV = new vertex_id_t[G->m];
    assert(G->endV != NULL);
    assert(fread(G->endV, sizeof(vertex_id_t), G->m, f) == G->m);    
    fclose(f);
}

/* write graph to the file */
void writeGraph(graph_t *G, char *filename)
{
    FILE *f = fopen(filename, "wb");
    assert(f != NULL);
    assert(fwrite(&G->n, sizeof(vertex_id_t), 1, f) == 1);
    assert(fwrite(&G->m, sizeof(edge_id_t), 1, f) == 1);
    unsigned char align = 0;
    assert(fwrite(&align, sizeof(unsigned char), 1, f) == 1);
    assert(fwrite(G->rowsIndices, sizeof(edge_id_t), G->n + 1, f) == G->n + 1);
    assert(fwrite(G->endV, sizeof(vertex_id_t), G->m, f) == G->m);
    fclose(f);
}

/* print graph (text format) */
void printGraph(graph_t *G)
{
    cout << "vertex number = " << G->n << ", edge number = " << G->m << " (each edge taken into account twice)" << endl;
    for (vertex_id_t i = 0; i < G->n; i++) {
        for (edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i + 1]; j++) {
            cout << i << " " << G->endV[j] << endl;
        }
    }
}

/* free graph */
void freeGraph(graph_t *G)
{
    delete[] G->rowsIndices;
    delete[] G->endV;
}
