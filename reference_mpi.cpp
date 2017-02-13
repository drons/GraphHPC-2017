#include "defs.h"

using namespace std;

struct msg_t
{
    vertex_id_t src_vertex, dest_vertex;
    unsigned paths_cnt;
};

/* distributed bfs */
/* bfs from the vertex start_v */
void bfs(graph_t *G, vertex_id_t start_v, unsigned **shortest_path, unsigned **shortest_paths_cnt)
{
    /* init */
    unsigned *d = new unsigned[G->local_n], *cnt = new unsigned[G->local_n];
    memset(d, 0xFF, sizeof(unsigned) * G->local_n);
    memset(cnt, 0, sizeof(unsigned) * G->local_n);
    
    MPI_Datatype MSG_T;
    MPI_Type_contiguous(3, MPI_UNSIGNED, &MSG_T);
    MPI_Type_commit(&MSG_T);
    
    vector<msg_t> send_msgs[G->nproc], recv_msgs[G->nproc];
    unsigned *sizes = new unsigned[G->nproc];
    MPI_Request *requests = new MPI_Request[G->nproc];
    MPI_Status *statuses = new MPI_Status[G->nproc];
    
    /* number of the shortest paths to start_v from start_v is 1 */
    if (G->rank == VERTEX_OWNER(start_v, G->n, G->nproc)) {
        vertex_id_t local_v = VERTEX_LOCAL(start_v, G->n, G->nproc, G->rank);
        d[local_v] = 0;
        cnt[local_v] = 1;
    }
    
    /* queue for the bfs */
    queue<vertex_id_t> q[2];
    if (G->rank == VERTEX_OWNER(start_v, G->n, G->nproc)) {
        q[0].push(start_v);
    }
    
    /* main part */
    unsigned level = 0;
    while (1) {
        /* if all queues are empty then break */
        int q_size = q[level % 2].size();
        int sum_q_size;
        MPI_Allreduce(&q_size, &sum_q_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (!sum_q_size) {
            break;
        }
        
        /* view the current queue (boundary vertices) */
        while (!q[level % 2].empty()) {
            /* u - boundary vertex */
            vertex_id_t u = q[level % 2].front();
            vertex_id_t local_u = VERTEX_LOCAL(u, G->n, G->nproc, G->rank);
            q[level % 2].pop();
            /* view neighbors of u and send them a message */
            for (edge_id_t j = G->rowsIndices[local_u]; j < G->rowsIndices[local_u + 1]; j++) {
                msg_t msg;
                msg.src_vertex = u;
                msg.dest_vertex = G->endV[j];
                msg.paths_cnt = cnt[local_u];
                send_msgs[VERTEX_OWNER(G->endV[j], G->n, G->nproc)].push_back(msg);
            }
        }
        
        /* send aggregated messages */
        for (int i = 0; i < G->nproc; i++) {
            MPI_Irecv(&sizes[i], 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD, &requests[i]);
        }
        
        for (int i = 0; i < G->nproc; i++) {
            unsigned msg_size = send_msgs[i].size();
            MPI_Send(&msg_size, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD);
        }
        
        MPI_Waitall(G->nproc, requests, statuses);
        
        for (int i = 0; i < G->nproc; i++) {
            recv_msgs[i].resize(sizes[i]);
            MPI_Irecv(&recv_msgs[i][0], sizes[i], MSG_T, i, 0, MPI_COMM_WORLD, &requests[i]);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        for (int i = 0; i < G->nproc; i++) {
            MPI_Send(&send_msgs[i][0], send_msgs[i].size(), MSG_T, i, 0, MPI_COMM_WORLD);
        }
        
        MPI_Waitall(G->nproc, requests, statuses);
        
        /* processing of received messages */
        for (int i = 0; i < G->nproc; i++) {
            for (int j = 0; j < (int)sizes[i]; j++) {
                msg_t msg = recv_msgs[i][j];
                vertex_id_t local_v = VERTEX_LOCAL(msg.dest_vertex, G->n, G->nproc, G->rank);
                if (level + 1 < d[local_v]) {
                    d[local_v] = level + 1;
                    cnt[local_v] = msg.paths_cnt;
                    q[(level + 1) % 2].push(msg.dest_vertex);
                } else if (level + 1 == d[local_v]) {
                    cnt[local_v] += msg.paths_cnt;
                }
            }
        }
        
        for (int i = 0; i < G->nproc; i++) {
            send_msgs[i].clear();
            recv_msgs[i].clear();
        }
        
        level++;
    }
    
    /* fill input arrays */
    for (vertex_id_t i = 0; i < G->local_n; i++) {
        shortest_path[i][start_v] = d[i];
        shortest_paths_cnt[i][start_v] = cnt[i];
    }
    
    /* free */
    delete[] d;
    delete[] cnt;
    delete[] sizes;
    delete[] requests;
    delete[] statuses;
}

/* algorithm */
void run(graph_t *G, double *result)
{
    /* init */
    unsigned **shortest_path = new unsigned *[G->local_n], **shortest_paths_cnt = new unsigned *[G->local_n];
    assert(shortest_path != NULL);
    assert(shortest_paths_cnt != NULL);
    for (vertex_id_t i = 0; i < G->local_n; i++) {
        shortest_path[i] = new unsigned[G->n];
        assert(shortest_path[i] != NULL);
        shortest_paths_cnt[i] = new unsigned[G->n];
        assert(shortest_paths_cnt[i] != NULL);
    }
    
    /* bfs for calculating shortest_path and shortest_paths_cnt for all pairs of vertices */
    for (vertex_id_t i = 0; i < G->n; i++) {
        bfs(G, i, shortest_path, shortest_paths_cnt);
    }
    
    /* main part */
    for (vertex_id_t i = 0; i < G->n; i++) {
        /* first send info about vertex i, i.e. shortest_path and shortest_paths_cnt from i to all other vertices */
        vertex_id_t local_i = VERTEX_LOCAL(i, G->n, G->nproc, G->rank);
        int owner_i = VERTEX_OWNER(i, G->n, G->nproc);
        unsigned *d, *cnt;
        if (G->rank == owner_i) {
            d = shortest_path[local_i];
            cnt = shortest_paths_cnt[local_i];
            for (int j = 0; j < G->nproc; j++) {
                if (j != G->rank) {
                    MPI_Send(d, G->n, MPI_UNSIGNED, j, 0, MPI_COMM_WORLD);
                    MPI_Send(cnt, G->n, MPI_UNSIGNED, j, 1, MPI_COMM_WORLD);
                }
            }
        } else {
            d = new unsigned[G->n];
            cnt = new unsigned[G->n];
            MPI_Status status;
            MPI_Recv(d, G->n, MPI_UNSIGNED, owner_i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(cnt, G->n, MPI_UNSIGNED, owner_i, 1, MPI_COMM_WORLD, &status);
        }
        
        /* update result */
        for (vertex_id_t j = i + 1; j < G->n; j++) {
            for (vertex_id_t u = 0; u < G->local_n; u++) {
                vertex_id_t global_u = VERTEX_TO_GLOBAL(u, G->n, G->nproc, G->rank);
                if (global_u == i || global_u == j) {
                    continue;
                }
                
                /* if one of the shortest paths between i and j contains vertex u */
                if (shortest_path[u][i] + shortest_path[u][j] == d[j]) {
                    result[u] += (double)shortest_paths_cnt[u][i] * shortest_paths_cnt[u][j] / cnt[j];
                }
            }
        }
        
        if (G->rank != owner_i) {
            delete[] d;
            delete[] cnt;
        }
    }
    
    /* free */
    for (vertex_id_t i = 0; i < G->local_n; i++) {
        delete[] shortest_path[i];
        delete[] shortest_paths_cnt[i];
    }
    
    delete[] shortest_path;
    delete[] shortest_paths_cnt;
}
