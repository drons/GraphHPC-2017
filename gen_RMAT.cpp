#include "defs.h"

#define FILENAME_LEN 256

char outFilename[FILENAME_LEN];

using namespace std;

/* helper */
void usage(int argc, char **argv) 
{
    printf("RMAT graph generator\n");
    printf("Usage:\n");
    printf("%s -s <scale> [options]\n", argv[0]);
    printf("Options:\n");
    printf("   -s <scale>, number of vertices is 2^<scale>\n");
    printf("   -k <half the average vertex degree>, default value is 16\n");
    exit(1);
}

/* initialization */
void init(int argc, char **argv, graph_t *G)
{
    G->scale = -1;
    G->a = 0.45;
    G->b = 0.25;
    G->c = 0.15;
    G->permute_vertices = true;
    /* default value */
    G->avg_vertex_degree = AVG_VERTEX_DEGREE;
    if (argc == 1) {
        usage(argc, argv);
    }
    
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-s")) {
            G->scale = (int)atoi(argv[++i]);
        }
        
        if (!strcmp(argv[i], "-k")) {
            G->avg_vertex_degree = (int)atoi(argv[++i]);
        }
    }
    
    sprintf(outFilename, "rmat-%d", G->scale);
    
    if (G->scale == -1) {
        usage(argc, argv);
    }
    
    G->n = (vertex_id_t)1 << G->scale;
    G->m = G->n * G->avg_vertex_degree;
}

/* RMAT graph generator */
void gen_RMAT_graph(graph_t *G)
{
    /* init */
    vertex_id_t n;
    edge_id_t m;
    edge_id_t offset;
    double a, b, c, d;
    double av, bv, cv, dv, S, p;
    int SCALE;
    double var;
    vertex_id_t step;
    bool permute_vertices;
    vertex_id_t *permV, tmpVal;
    vertex_id_t u, v;
    vertex_id_t *src;
    vertex_id_t *dest;
    unsigned *degree;
    
    a = G->a;
    b = G->b;
    c = G->c;
    assert(a + b + c < 1);
    d = 1 - (a + b + c);
    
    permute_vertices = G->permute_vertices;
    
    n = G->n;
    m = G->m;
    
    src = new vertex_id_t[m];
    assert(src != NULL);
    dest = new vertex_id_t[m];
    assert(dest != NULL);
    degree = new unsigned[n];
    assert(degree != NULL);
    memset(degree, 0, n * sizeof(unsigned));
    
    /* fixed graph */
    srand48(2387);
    SCALE = G->scale;
    
    /* generate edges */
    edge_id_t edges_cnt = 0;
    while (edges_cnt < m) {
        u = 1;
        v = 1;
        step = n / 2;
        
        av = a;
        bv = b;
        cv = c;
        dv = d;
        
        p = drand48();
        if (p < av) {
            /* Do nothing */
        } else if ((p >= av) && (p < av + bv)) {
            v += step;
        } else if ((p >= av + bv) && (p < av + bv + cv)) {
            u += step;
        } else {
            u += step;
            v += step;
        }
        
        for (int j = 1; j < SCALE; j++) {
            step = step / 2;
            
            /* Vary a, b, c, d by up to 10% */
            var = 0.1;
            av *= 0.95 + var * drand48();
            bv *= 0.95 + var * drand48();
            cv *= 0.95 + var * drand48();
            dv *= 0.95 + var * drand48();
            
            S = av + bv + cv + dv;
            av = av / S;
            bv = bv / S;
            cv = cv / S;
            dv = dv / S;
            
            /* Choose partition */
            p = drand48();
            if (p < av) {
                /* Do nothing */
            } else if ((p >= av) && (p < av + bv)) {
                v += step;
            } else if ((p >= av + bv) && (p < av + bv + cv)) {
                u += step;
            } else {
                u += step;
                v += step;
            }
        }
        
        /* graph without self-loops */
        if (u != v) {
            src[edges_cnt] = u - 1;
            dest[edges_cnt++] = v - 1;
        }
    }
    
    /* reshuffle */
    if (permute_vertices) {
        /* fixed graph */
        srand48(4791);
        
        permV = new vertex_id_t[n];
        assert(permV != NULL);
        
        for (vertex_id_t i = 0; i < n; i++) {
            permV[i] = i;
        }
        
        for (vertex_id_t i = 0; i < n; i++) {
            vertex_id_t j = n * drand48();
            tmpVal = permV[i];
            permV[i] = permV[j];
            permV[j] = tmpVal;
        }
        
        for (edge_id_t i = 0; i < m; i++) {
            src[i] = permV[src[i]];
            dest[i] = permV[dest[i]];
        }
        
        delete[] permV;
    }
    
    /* Update graph data structure */
    for (edge_id_t i = 0; i < m; i++) {
        degree[src[i]]++;
        degree[dest[i]]++;
    }
    
    G->endV = new vertex_id_t[2 * m];
    assert(G->endV != NULL);
    
    G->rowsIndices = new edge_id_t[n + 1];
    assert(G->rowsIndices != NULL);
    
    G->n = n;
    /* undirected graph, each edge is stored twice; if edge is (u, v), then it's
     *stored at the vertex u and at the vertex v */
    G->m = 2 * m;
    
    G->rowsIndices[0] = 0;
    for (vertex_id_t i = 1; i <= G->n; i++) {
        G->rowsIndices[i] = G->rowsIndices[i - 1] + degree[i - 1];
    }
    
    for (edge_id_t i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];
        offset = degree[u]--;
        G->endV[G->rowsIndices[u] + offset - 1] = v;
        offset = degree[v]--;
        G->endV[G->rowsIndices[v] + offset - 1] = u;
    }
    
    delete[] src;
    delete[] dest;
    delete[] degree;
}

int main(int argc, char **argv)
{
    graph_t g;
    init(argc, argv, &g);
    gen_RMAT_graph(&g);
    writeGraph(&g, outFilename);
    
    return 0;
}
