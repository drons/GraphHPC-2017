#include <mpi.h>
#include "solution.h"

void run_mpi( graph_t* G, double* result )
{
    int mpi_size = -1;
    int mpi_rank = -1;
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

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

    vertex_id_t part_size = std::max( (vertex_id_t)(n/mpi_size), (vertex_id_t)1 );
    vertex_id_t part_begin = part_size*mpi_rank;
    vertex_id_t part_end = std::min( part_begin + part_size, G->n );

//    std::cout << std::endl;
//    std::cout << " mpi_rank = " << mpi_rank
//              << " mpi_size = " << mpi_size
//              << " part_size = " << part_size
//              << " part_begin = " << part_begin
//              << " part_end = " << part_end
//              << " omp_get_max_threads " << omp_get_max_threads() << std::endl;
//    MPI_Barrier( MPI_COMM_WORLD );

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

    MPI_Reduce( local_result.data(), result, G->n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD );
}
