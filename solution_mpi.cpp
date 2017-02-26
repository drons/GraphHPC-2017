#include "solution.h"

vertex_id_t get_local_vert_count( vertex_id_t TotVertices, unsigned size, unsigned rank )
{
    vertex_id_t mod_size = MOD_SIZE(TotVertices);
    vertex_id_t div_size = DIV_SIZE(TotVertices);
    if (!mod_size) {
        return div_size;
    } else {
        if (rank < mod_size) {
            return div_size + 1;
        } else {
            return div_size;
        }
    }
}

vertex_id_t get_local_part_offset( const vertex_id_t TotVertices, const int size, const int rank)
{
    if (MOD_SIZE(TotVertices) > (unsigned int)rank) {
        return ((DIV_SIZE(TotVertices) + 1) * rank );
    } else {
        return (MOD_SIZE(TotVertices) * (DIV_SIZE(TotVertices) + 1) + DIV_SIZE(TotVertices) * (rank - MOD_SIZE(TotVertices)) );
    }
}

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
//    MPI_Barrier(MPI_COMM_WORLD);
//    double wt( omp_get_wtime() );

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

//    MPI_Barrier(MPI_COMM_WORLD);
//    if( mpi_rank == 0 )
//    {
//        std::cout << "Join time" << omp_get_wtime() - wt << std::endl;
//    }

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
            std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
        }
    }*/

    vertex_id_t                     n = G->n;
    compute_buffer_t*               buffers;
    std::vector<uint32_t>           rows_indices32;
    int                             max_work_threads = omp_get_max_threads();
    std::vector< vertex_id_t>       map( sort_graph( G, G, 0 ) );

    rows_indices32.resize( G->n + 1 );
    for( size_t n = 0; n < G->n + 1; ++n )
    {
        rows_indices32[n] = G->rowsIndices[n];
    }

    buffers = new compute_buffer_t[ max_work_threads ];

    #pragma omp parallel for
    for( int t = 0; t < max_work_threads; ++t )
    {
        compute_buffer_t&   b( buffers[ t ] );
        b.resize( G );
    }

    vertex_id_t     part_size = get_local_n( g_local );
    vertex_id_t     part_begin = get_local_part_offset( G->n, G->nproc, G->rank );
    vertex_id_t     part_end = std::min( part_begin + part_size, G->n );

    //try to balance
    {
        std::vector<int>    num_threads;
        int                 total_threads = 0;
        num_threads.resize( mpi_size );

        MPI_Allgather( &max_work_threads, 1, MPI::INT, num_threads.data(), 1, MPI::INT, MPI_COMM_WORLD );
        for( size_t i = 0; i != num_threads.size(); ++i )
            total_threads += num_threads[i];

        if( mpi_rank == 0 )
        {
            std::cout << "t = { ";
            for( size_t i = 0; i != num_threads.size(); ++i )
            {
                std::cout << num_threads[i] << " ";
            }
            std::cout << "}" << std::endl;
        }
        std::vector<int>    tpart_sizes;
        std::vector<int>    tpart_offsets;
        std::vector<int>    part_sizes;
        std::vector<int>    part_offsets;
        tpart_sizes.resize( total_threads );
        tpart_offsets.resize( total_threads );

        part_sizes.resize( mpi_size );
        part_offsets.resize( mpi_size );

        for( size_t i = 0; i != total_threads; ++i )
        {
            tpart_sizes[i] = get_local_vert_count( G->n, total_threads, i );
            tpart_offsets[i] = get_local_part_offset( G->n, G->nproc, G->rank );
        }

        if( mpi_rank == 0 )
        {
            std::cout << "total_threads = " << total_threads << std::endl;
        }

        size_t ct = 0;
        size_t off = 0;
        for( size_t r = 0; r != mpi_size; ++r )
        {
            int s = 0;
            for( size_t t = 0; t != num_threads[r]; ++t )
            {
                s += tpart_sizes[ct];
                ++ct;
            }

            if( mpi_rank == 0 )
            {
                std::cout << " rank = " << r
                          << " parts = " << s << " x " << get_local_n( g_local )
                          << " offs = " << off << " x " << get_local_part_offset( G->n, G->nproc, r )
                          << std::endl;
            }

            part_sizes[r] = s;
            part_offsets[r] = off;
            off += s;
        }
        part_size = part_sizes[mpi_rank];
        part_begin = part_offsets[mpi_rank];
        part_end = std::min( part_begin + part_size, G->n );
    }
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
        bfs( G, rows_indices32.data(), s, b.distance, b.shortest_count, b.q, b.qnext, b.vertex_on_level_count, max_distance );
        betweenness_centrality( G, rows_indices32.data(), s, b.distance, b.shortest_count, b.vertex_on_level_count, max_distance, b.delta, b.delta_precompute, b.partial_result );
    }

    std::vector<double> local_result;

    local_result.resize( G->n );
    std::fill( local_result.begin(), local_result.end(), 0 );

    #pragma omp parallel for
    for( vertex_id_t s = 0; s < n; ++s )
    {
        double  r = 0;
        for( int t = 0; t < max_work_threads; ++t )
        {
            compute_buffer_t&   b( buffers[t] );
            r += b.partial_result[s];
        }
        local_result[map[s]] = r*0.5;
    }

    for( int t = 0; t < max_work_threads; ++t )
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
    delete [] buffers;
}
