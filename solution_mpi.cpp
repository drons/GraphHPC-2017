#include "solution.h"

void run_mpi( graph_t* g_local, double* result )
{
    graph_t     storG;
    graph_t*    G = NULL;
    int         mpi_size = -1;
    int         mpi_rank = -1;

    MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

    G = &storG;
    G->n = g_local->n;
    G->m = g_local->m;
    G->nproc = g_local->nproc;
    G->rank = g_local->rank;
    G->endV = new vertex_id_t[ G->m ];
    G->rowsIndices = new edge_id_t[ G->n + 1 ];

    vertex_id_t edge_end = 0;
    for( vertex_id_t v = 0; v < G->n; ++v )
    {
        int             mpi_owner = VERTEX_OWNER( v, G->n, G->nproc );
        vertex_id_t     v_local = VERTEX_LOCAL( v, G->n, G->nproc, mpi_owner );
        vertex_id_t     edge_cnt;

        if( mpi_rank == mpi_owner )
        {
            edge_cnt = g_local->rowsIndices[v_local + 1] - g_local->rowsIndices[v_local];
        }

        MPI_Bcast( &edge_cnt, 1, MPI::UNSIGNED, mpi_owner, MPI_COMM_WORLD );

        G->rowsIndices[v] = edge_end;

        if( edge_cnt == 0 )
        {
            continue;
        }

        if( mpi_rank == mpi_owner )
        {
            memcpy( G->endV + edge_end, g_local->endV + g_local->rowsIndices[v_local], sizeof( vertex_id_t )*edge_cnt );
        }

        MPI_Bcast( G->endV + edge_end, edge_cnt, MPI::UNSIGNED, mpi_owner, MPI_COMM_WORLD );

        edge_end += edge_cnt;
    }

    G->rowsIndices[G->n] = edge_end;

    /*
    if( mpi_rank == 0 )
    {
        graph_t Gt;
        std::cout << "read start...";
        readGraph( &Gt, "/home/sudorgin/g/rmat-4" );
        std::cout << "ok\n";
        if( 0 != memcmp( G->rowsIndices, Gt.rowsIndices, sizeof( edge_id_t )*(G->n + 1) ) ||
            0 != memcmp( G->endV, Gt.endV, sizeof( vertex_id_t )*(G->m) ) )
        {
            //exit(-1);
            std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
        }
    }*/

    vertex_id_t                     n = G->n;
    std::vector<compute_buffer_t>   buffers;
    std::vector<uint32_t>           rows_indices32;
    size_t                          max_work_threads = omp_get_max_threads();

    rows_indices32.resize( G->n + 1 );
    for( size_t n = 0; n < G->n + 1; ++n )
    {
        rows_indices32[n] = G->rowsIndices[n];
    }

    buffers.resize( max_work_threads );

    #pragma omp parallel for
    for( size_t t = 0; t < max_work_threads; ++t )
    {
        compute_buffer_t&   b( buffers[ t ] );
        b.resize( G );
    }

    vertex_id_t part_size = get_local_n( g_local );
    vertex_id_t part_begin = VERTEX_TO_GLOBAL( 0, G->n, G->nproc, G->rank );
    vertex_id_t part_end = std::min( part_begin + part_size, G->n );
/*
    if( mpi_rank == 0 )
    {
        std::cout << std::endl;
    }
{
        std::vector<vertex_id_t>    x;
        x.resize( mpi_size );

        MPI_Allgather( &part_size, 1, MPI::UNSIGNED, x.data(), 1, MPI::UNSIGNED, MPI_COMM_WORLD );

        if( mpi_rank == 0 )
        {
            std::cout << "part_size = ";
            for( int r = 0; r != mpi_size; ++r )
            {
                std::cout << x[r] << " ";
            }
            std::cout << std::endl;
        }
        MPI_Allgather( &part_begin, 1, MPI::UNSIGNED, x.data(), 1, MPI::UNSIGNED, MPI_COMM_WORLD );

        if( mpi_rank == 0 )
        {
            std::cout << "part_begin = ";
            for( int r = 0; r != mpi_size; ++r )
            {
                std::cout << x[r] << " ";
            }
            std::cout << std::endl;
        }
        MPI_Allgather( &part_end, 1, MPI::UNSIGNED, x.data(), 1, MPI::UNSIGNED, MPI_COMM_WORLD );

        if( mpi_rank == 0 )
        {
            std::cout << "part_end = ";
            for( int r = 0; r != mpi_size; ++r )
            {
                std::cout << x[r] << " ";
            }
            std::cout << std::endl;
        }
}
    return;
*/
    #pragma omp parallel for
    for( vertex_id_t s = part_begin; s < part_end; ++s )
    {
        compute_buffer_t&   b( buffers[ omp_get_thread_num() ] );
        DIST_TYPE           max_distance = 0;

        std::fill( b.vertex_on_level_count, b.vertex_on_level_count + b.max_distance, 0 );
        bfs( G, rows_indices32.data(), s, b.distance, b.shortest_count, b.q, b.qnext, b.vertex_on_level_count, b.global_vertex_on_level_count, b.global_unmarked_vertex_count, max_distance );
        betweenness_centrality( G, rows_indices32.data(), s, b.distance, b.shortest_count, b.vertex_on_level_count, max_distance, b.delta, b.partial_result );

        vertex_id_t unmarked = G->n;
        for( size_t distance = 0; distance != max_distance; ++distance )
        {
            unmarked -= b.vertex_on_level_count[distance];
            b.global_vertex_on_level_count[distance] += b.vertex_on_level_count[distance];
            b.global_unmarked_vertex_count[distance] += unmarked;
        }
    }

    std::vector<double> local_result;

    local_result.resize( G->n );
    std::fill( local_result.begin(), local_result.end(), 0 );

    #pragma omp parallel for
    for( vertex_id_t s = 0; s < n; ++s )
    {
        double  r = 0;
        for( size_t t = 0; t < max_work_threads; ++t )
        {
            compute_buffer_t&   b( buffers[t] );
            r += b.partial_result[s];
        }
        local_result[s] = r*0.5;
    }

    for( size_t t = 0; t < max_work_threads; ++t )
    {
        compute_buffer_t&   b( buffers[ t ] );
        b.release();
    }

    std::vector<int>    num_vert;
    int local_n = g_local->local_n;

    num_vert.resize( mpi_size );

    MPI_Allgather( &local_n, 1, MPI::INT, num_vert.data(), 1, MPI::INT, MPI_COMM_WORLD );
    MPI_Reduce_scatter( local_result.data(), result, num_vert.data(), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );

    freeGraph( G );
}
