#include "solution.h"

void simplified_dijkstra( const graph_t* G, const uint32_t* row_indites, vertex_id_t start, DIST_TYPE* distance, SCOUNT_TYPE* shortest_count, wavefront_t& queue )
{
    static const DIST_TYPE      INVALID_DISTANCE = std::numeric_limits<DIST_TYPE>::max();

    memset( distance, INVALID_DISTANCE, sizeof(DIST_TYPE)*G->n );
    memset( shortest_count, 0, sizeof(SCOUNT_TYPE)*G->n );

    distance[ start ] = 0;
    shortest_count[ start ] = 1;

    queue.push_back( start );

    while( !queue.empty() )
    {
        vertex_id_t     v = queue.front();
        DIST_TYPE       dist_v = distance[ v ];
        SCOUNT_TYPE     shortest_count_v( shortest_count[v] );

        queue.pop_front();

        vertex_id_t* ibegin = G->endV + row_indites[ v ];
        vertex_id_t* iend = G->endV + row_indites[ v + 1 ];

        for( vertex_id_t* e = ibegin; e != iend; ++e )
        {
            vertex_id_t w( *e );
            DIST_TYPE   new_dist_v = dist_v + 1;
            DIST_TYPE&  dist_w( distance[ w ] );
            if( dist_w == INVALID_DISTANCE )
            {
                queue.push_back( w );
                dist_w = new_dist_v;
            }
            if( dist_w == new_dist_v )
            {
                shortest_count[w] += shortest_count_v;
            }
        }
    }
    return;
}

void bfs( const graph_t* G, const uint32_t* row_indites, vertex_id_t start,
          DIST_TYPE* distance, SCOUNT_TYPE* shortest_count,
          wavefront_t& q, wavefront_t& qnext,
          vertex_id_t* vertex_on_level_count,
          const double* global_vertex_on_level_count,
          const double* global_unmarked_vertex_count, DIST_TYPE& max_distance )
{
    static const DIST_TYPE      INVALID_DISTANCE = std::numeric_limits<DIST_TYPE>::max();

    memset( distance, INVALID_DISTANCE, sizeof(DIST_TYPE)*G->n );
    memset( shortest_count, 0, sizeof(SCOUNT_TYPE)*G->n );

    distance[ start ] = 0;
    shortest_count[ start ] = 1;
    vertex_on_level_count[0] = 1;

    size_t      num_verts_on_level = 1;
    DIST_TYPE   current_level = 0;
    DIST_TYPE   next_level = 1;

    q.reset();
    qnext.reset();

    q.push_back( start );

    do
    {//unroll first iteration in classic mode
        num_verts_on_level = 0;
        const vertex_id_t* rend( q.rend() );
        for( const vertex_id_t* ii = q.rbegin(); ii != rend; --ii )
        {
            vertex_id_t  v = *ii;
            vertex_id_t* ibegin = G->endV + row_indites[ v ];
            vertex_id_t* iend = G->endV + row_indites[ v + 1 ];

            for( vertex_id_t* e = ibegin; e != iend; ++e )
            {
                vertex_id_t w( *e );
                DIST_TYPE&  distance_w(distance[w]);
                if( distance_w == INVALID_DISTANCE )
                {
                    distance_w = next_level;
                    shortest_count[w] = shortest_count[v];
                    qnext.push_back( w );
                }
                else
                if( distance_w == next_level )
                {
                    shortest_count[w] += shortest_count[v];
                }
            }
        }
        ++current_level;
        ++next_level;
        num_verts_on_level = qnext.size();
        vertex_on_level_count[current_level] = num_verts_on_level;
        q.swap( qnext );
        qnext.reset();
    }while( num_verts_on_level > 0 && num_verts_on_level < 1024 );

    while( num_verts_on_level > 0 )
    {
        if( global_vertex_on_level_count[current_level] <
            global_unmarked_vertex_count[current_level] )
        {//descending
//          std::cout << "D";
            num_verts_on_level = 0;
            for( vertex_id_t v = 0; v != G->n; ++v )
            {
                const DIST_TYPE   distance_v( distance[v] );
                if( distance_v == current_level )
                {
                    const vertex_id_t* ibegin = G->endV + row_indites[ v ];
                    const vertex_id_t* iend = G->endV + row_indites[ v + 1 ];

                    for( const vertex_id_t* e = ibegin; e != iend; ++e )
                    {
                        vertex_id_t w( *e );
                        DIST_TYPE&  distance_w(distance[w]);
                        if( distance_w == INVALID_DISTANCE )
                        {
                            distance_w = next_level;
                            ++num_verts_on_level;
                            shortest_count[w] = shortest_count[v];
                        }
                        else
                        if( distance_w == next_level )
                        {
                            shortest_count[w] += shortest_count[v];
                        }
                    }
                }
            }
        }
        else
        {//ascending
//          std::cout << "A";
            num_verts_on_level = 0;
            for( vertex_id_t v = 0; v != G->n; ++v )
            {
                DIST_TYPE&   distance_v( distance[v] );
                if( distance_v == INVALID_DISTANCE )
                {
                    const vertex_id_t* ibegin = G->endV + row_indites[ v ];
                    const vertex_id_t* iend = G->endV + row_indites[ v + 1 ];

                    for( const vertex_id_t* e = ibegin; e != iend; ++e )
                    {
                        vertex_id_t     w( *e );
                        const DIST_TYPE distance_w( distance[w] );

                        if( distance_w == current_level )
                        {
                            if( distance_v == INVALID_DISTANCE )
                            {
                                distance_v = next_level;
                                shortest_count[v] = shortest_count[w];
                                ++num_verts_on_level;
                            }
                            else
                            {
                                shortest_count[v] += shortest_count[w];
                            }
                        }
                    }
                }
            }
        }
        ++current_level;
        ++next_level;
        vertex_on_level_count[current_level] = num_verts_on_level;
    }
    max_distance = current_level;
}

void betweenness_centrality( graph_t* G, const uint32_t* row_indites, vertex_id_t s,
                             const DIST_TYPE* distance,
                             const SCOUNT_TYPE* shortest_count,
                             vertex_id_t* vertex_on_level_count, DIST_TYPE max_distance,
                             DELTA_TYPE* delta,
                             PARTIAL_TYPE* result )
{
    if( max_distance <= 1 )
        return;
    memset( delta, 0, sizeof(DELTA_TYPE)*G->n );

    DIST_TYPE   next_distance = max_distance;
    --max_distance;
    for(;;)
    {
//        std::cout << "level = " << max_distance << " " << vertex_on_level_count[ max_distance + 1 ]
//                  << " " << vertex_on_level_count[ max_distance ] << std::endl;
        if( vertex_on_level_count[ next_distance ] <
            vertex_on_level_count[ max_distance ] )
        {
            for( vertex_id_t w = 0; w != G->n; ++w )
            {
                if( distance[w] != next_distance )
                {
                    continue;
                }
                const vertex_id_t*  ibegin = G->endV + row_indites[ w ];
                const vertex_id_t*  iend = G->endV + row_indites[ w + 1 ];

                for( const vertex_id_t* e = ibegin; e != iend; ++e )
                {
                    vertex_id_t v( *e );
                    if( max_distance == distance[v] )
                    {
                        const PARTIAL_TYPE sc_v( ((PARTIAL_TYPE)shortest_count[v]) );
                        delta[v] += sc_v*(1 + ((PARTIAL_TYPE)delta[w]))/((PARTIAL_TYPE)shortest_count[w]);
                    }
                }
            }
        }
        else
        {
            for( vertex_id_t v = 0; v != G->n; ++v )
            {
                if( distance[v] != max_distance )
                {
                    continue;
                }
                const vertex_id_t*  ibegin = G->endV + row_indites[ v ];
                const vertex_id_t*  iend = G->endV + row_indites[ v + 1 ];

                for( const vertex_id_t* e = ibegin; e != iend; ++e )
                {
                    vertex_id_t w( *e );

                    if( distance[w] == next_distance )
                    {
                        delta[v] += ((PARTIAL_TYPE)shortest_count[v])*(1 + ((PARTIAL_TYPE)delta[w]))/(PARTIAL_TYPE)shortest_count[w];
                    }
                }
            }
        }
        if( max_distance == 1 )
        {
            break;
        }
        --max_distance;
        --next_distance;
    }

    for( vertex_id_t w = 0; w != G->n; ++w )
    {
        if( w != s )
        {
            result[w] += delta[w];
        }
    }
}

void run( graph_t* G, double* result )
{
    vertex_id_t                     n = G->n;
    std::vector<compute_buffer_t>   buffers;
    std::vector<uint32_t>           rows_indices32;
    size_t                          max_work_threads = omp_get_max_threads();

    rows_indices32.resize( G->n + 1 );
    for( size_t n = 0; n < G->n + 1; ++n )
    {
        rows_indices32[n] = G->rowsIndices[n];
    }

    buffers.resize( omp_get_max_threads() );

    #pragma omp parallel for
    for( size_t t = 0; t < max_work_threads; ++t )
    {
        compute_buffer_t&   b( buffers[ omp_get_thread_num() ] );
        b.resize( G );
    }

//    std::cout << "omp_get_num_procs   " << omp_get_num_procs() << std::endl;
//    std::cout << "omp_get_num_threads " << omp_get_num_threads() << std::endl;
//    std::cout << "omp_get_max_threads " << omp_get_max_threads() << std::endl;
//    std::cout << "omp_get_num_teams   " << omp_get_num_teams() << std::endl;
//    std::cout << "omp_get_team_num    " << omp_get_team_num() << std::endl;

#ifdef MAXNODES
    n = MAXNODES;
#endif
    #pragma omp parallel for schedule ( dynamic )
    for( vertex_id_t s = 0; s < n; ++s )
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

        if(0)
        {
            compute_buffer_t    t;
            t.resize( G );
//            std::cout << "bfs";
//            b.dump_bfs_result( s, max_distance );
            simplified_dijkstra( G, rows_indices32.data(), s, t.distance, t.shortest_count, t.q );
//            std::cout << "simplified_dijkstra";
//            b.dump_bfs_result( s, max_distance );
            if( !t.is_equal( b ) )
            {
                std::cout << "s = " << s;
                exit( -1 );
            }
        }
    }

    #pragma omp parallel for
    for( vertex_id_t s = 0; s < n; ++s )
    {
        double  r = 0;
        for( size_t t = 0; t < max_work_threads; ++t )
        {
            compute_buffer_t&   b( buffers[ t ] );
            r += b.partial_result[s];
        }
        result[s] = r*0.5;
    }

    for( size_t t = 0; t < max_work_threads; ++t )
    {
        compute_buffer_t&   b( buffers[ t ] );
        b.release();
    }
}
