#include "solution.h"

#include <mmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

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

#define UNROLL 16
//int vector_hit_count = 0;
//int scalar_hit_count = 0;
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

    size_t      n = G->n;
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
            size_t  v = *ii;
            vertex_id_t* ibegin = G->endV + row_indites[ v ];
            vertex_id_t* iend = G->endV + row_indites[ v + 1 ];

            for( vertex_id_t* e = ibegin; e != iend; ++e )
            {
                size_t w( *e );
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
            num_verts_on_level = 0;
#ifdef UNROLL
            __m128i current_level16( _mm_set1_epi8( current_level ) );
            for( size_t vu = 0; vu < n; vu += UNROLL )
            {
                __m128i distance_v16( _mm_load_si128( (__m128i*)(distance + vu) ) );//load from mem
                __m128i comp( _mm_cmpeq_epi8( current_level16, distance_v16 ) );    //compare
                if( _mm_testz_si128( comp, comp ) )// ( (!_mm_testz_si128(a,a))  == horizontal OR )
                {
                    continue;
                }

                uint8_t b[UNROLL] __attribute__ ((aligned (16)));
                _mm_store_si128((__m128i*)b, comp );
                for( size_t u = 0; u < UNROLL; ++u )
                {
                    size_t v = vu + u;
                    if( b[u] )
#else
            {
                for( size_t v = 0; v < n; ++v )
                {
                    if( distance[v] == current_level )
#endif
                    {
                        const vertex_id_t* ibegin = G->endV + row_indites[ v ];
                        const vertex_id_t* iend = G->endV + row_indites[ v + 1 ];

                        for( const vertex_id_t* e = ibegin; e != iend; ++e )
                        {
                            size_t w( *e );
                            if( distance[w] == INVALID_DISTANCE )
                            {
                                distance[w] = next_level;
                                ++num_verts_on_level;
                                shortest_count[w] = shortest_count[v];
                            }
                            else
                            if( distance[w] == next_level )
                            {
                                shortest_count[w] += shortest_count[v];
                            }
                        }
                    }
                }
            }
        }
        else
        {//ascending
            num_verts_on_level = 0;
#ifdef UNROLL
            __m128i invalid_distance16( _mm_set1_epi8( INVALID_DISTANCE ) );
            for( size_t vu = 0; vu < n; vu += UNROLL )
            {
                __m128i distance_v16( _mm_load_si128( (__m128i*)(distance + vu) ) );//load from mem
                __m128i comp( _mm_cmpeq_epi8( invalid_distance16, distance_v16 ) );    //compare
                if( _mm_testz_si128( comp, comp ) )// ( (!_mm_testz_si128(a,a))  == horizontal OR )
                {
                    continue;
                }

                uint8_t b[UNROLL] __attribute__ ((aligned (16)));
                _mm_store_si128((__m128i*)b, comp );
                for( size_t u = 0; u < UNROLL; ++u )
                {
                    size_t v = vu + u;
                    if( b[u] )
#else
            {
                for( size_t v = 0; v < n; ++v )
                {
                    if( distance[v] == INVALID_DISTANCE )
#endif
                    {
                        const vertex_id_t* e = G->endV + row_indites[ v ];
                        const vertex_id_t* iend = G->endV + row_indites[ v + 1 ];

                        for( ; e < iend; ++e )
                        {
                            size_t     w( *e );
                            if( distance[w] == current_level )
                            {
                                distance[v] = next_level;
                                shortest_count[v] = shortest_count[w];
                                ++num_verts_on_level;
                                break;
                            }
                        }
                        ++e;
                        for( ; e < iend; ++e )
                        {
                            size_t     w( *e );
                            if( distance[w] == current_level )
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

    size_t      n = G->n;
    DIST_TYPE   next_distance = max_distance;
    --max_distance;
    for(;;)
    {
        if( vertex_on_level_count[ next_distance ] <
            vertex_on_level_count[ max_distance ] )
        {
#ifdef UNROLL
            __m128i next_distance16( _mm_set1_epi8( next_distance ) );
            for( size_t wu = 0; wu < n; wu += UNROLL )
            {
                __m128i distance_v16( _mm_load_si128( (__m128i*)(distance + wu) ) );//load from mem
                __m128i comp( _mm_cmpeq_epi8( next_distance16, distance_v16 ) );    //compare
                if( _mm_testz_si128( comp, comp ) )// ( (!_mm_testz_si128(a,a))  == horizontal OR )
                {
                    continue;
                }

                uint8_t b[UNROLL] __attribute__ ((aligned (16)));
                _mm_store_si128((__m128i*)b, comp );
                for( size_t u = 0; u < UNROLL; ++u )
                {
                    size_t w = wu + u;
                    if( b[u] )
#else
            {
                for( size_t w = 0; w < n; ++w )
                {
                    if( distance[w] == next_distance )
#endif
                    {
                        const vertex_id_t*  ibegin = G->endV + row_indites[ w ];
                        const vertex_id_t*  iend = G->endV + row_indites[ w + 1 ];

                        for( const vertex_id_t* e = ibegin; e != iend; ++e )
                        {
                            size_t v( *e );
                            if( max_distance == distance[v] )
                            {
                                const PARTIAL_TYPE sc_v( ((PARTIAL_TYPE)shortest_count[v]) );
                                delta[v] += sc_v*(1 + ((PARTIAL_TYPE)delta[w]))/((PARTIAL_TYPE)shortest_count[w]);
                            }
                        }
                    }
                }
            }
        }
        else
        {
#ifdef UNROLL
            __m128i max_distance16( _mm_set1_epi8( max_distance ) );
            for( size_t vu = 0; vu < n; vu += UNROLL )
            {
                __m128i distance_v16( _mm_load_si128( (__m128i*)(distance + vu) ) );//load from mem
                __m128i comp( _mm_cmpeq_epi8( max_distance16, distance_v16 ) );    //compare
                if( _mm_testz_si128( comp, comp ) )// ( (!_mm_testz_si128(a,a))  == horizontal OR )
                {
                    continue;
                }

                uint8_t b[UNROLL] __attribute__ ((aligned (16)));
                _mm_store_si128((__m128i*)b, comp );
                for( size_t u = 0; u < UNROLL; ++u )
                {
                    size_t v = vu + u;
                    if( b[u] )
#else
            {
                for( size_t v = 0; v < n; ++v )
                {
                    if( distance[v] == max_distance )
#endif
                    {
                        const vertex_id_t*  ibegin = G->endV + row_indites[ v ];
                        const vertex_id_t*  iend = G->endV + row_indites[ v + 1 ];

                        for( const vertex_id_t* e = ibegin; e != iend; ++e )
                        {
                            size_t w( *e );

                            if( distance[w] == next_distance )
                            {
                                delta[v] += ((PARTIAL_TYPE)shortest_count[v])*(1 + ((PARTIAL_TYPE)delta[w]))/(PARTIAL_TYPE)shortest_count[w];
                            }
                        }
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

    for( size_t w = 0; w != n; ++w )
    {
        if( w != s )
        {
            result[w] += delta[w];
        }
    }
}

struct sorter_a
{
    const graph_t* G;
    sorter_a( const graph_t* g ) : G( g ){}
    bool operator () ( vertex_id_t v, vertex_id_t w ) const
    {
        return ( G->rowsIndices[v + 1] - G->rowsIndices[v] ) <
               ( G->rowsIndices[w + 1] - G->rowsIndices[w] );
    }
};

struct sorter_d
{
    const graph_t* G;
    sorter_d( const graph_t* g ) : G( g ){}
    bool operator () ( vertex_id_t v, vertex_id_t w ) const
    {
        return ( G->rowsIndices[v + 1] - G->rowsIndices[v] ) >
               ( G->rowsIndices[w + 1] - G->rowsIndices[w] );
    }
};

std::vector< vertex_id_t> sort_graph( graph_t* G, int order )
{
    std::vector< vertex_id_t>   fmap;
    std::vector< vertex_id_t>   imap;

    fmap.resize( G->n );
    imap.resize( G->n );

    for( size_t n = 0; n != G->n; ++n )
        fmap[n] = n;

    if( order == 0 )
    {
        return fmap;
    }

    if( order > 0 )
    {
        std::sort( fmap.begin(), fmap.end(), sorter_a( G ) );
    }
    else
    {
        std::sort( fmap.begin(), fmap.end(), sorter_d( G ) );
    }

    #pragma omp parallel for
    for( size_t i = 0; i < G->n; ++i )
    {
        vertex_id_t v = fmap[i];
        imap[v] = i;
    }
/*
    std::cout << "map = { ";
    for( size_t i = 0; i != G->n; ++i )
    {
        std::cout << fmap[i] << " ";
    }
    std::cout << "}" << std::endl;

    std::cout << "ecnt = { ";
    for( size_t i = 0; i != G->n; ++i )
    {
        vertex_id_t v = fmap[i];
        std::cout << G->rowsIndices[v + 1] - G->rowsIndices[v] << " ";
    }
    std::cout << "}" << std::endl;
*/
    edge_id_t*      new_rows_indices = new edge_id_t[ G->n + 1 ];
    vertex_id_t*    new_end_v = new vertex_id_t[ G->m ];

    edge_id_t   off = 0;
    for( vertex_id_t i = 0; i != G->n; ++i )
    {
        vertex_id_t     v( fmap[i] );
        vertex_id_t*    ibegin = G->endV + G->rowsIndices[ v ];
        vertex_id_t*    iend = G->endV + G->rowsIndices[ v + 1 ];
        vertex_id_t     sz( iend - ibegin );

        memcpy( new_end_v + off, ibegin, sizeof( vertex_id_t )*sz );
        new_rows_indices[i] = off;
        off += sz;
    }
    new_rows_indices[G->n] = off;

    #pragma omp parallel for
    for( vertex_id_t e = 0; e < G->m; ++e )
    {
        new_end_v[e] = imap[ new_end_v[e] ];
    }

    delete [] G->endV;
    delete [] G->rowsIndices;

    G->endV = new_end_v;
    G->rowsIndices = new_rows_indices;

    #pragma omp parallel for schedule ( dynamic )
    for( vertex_id_t v = 0; v < G->n; ++v )
    {
        vertex_id_t* ibegin = G->endV + G->rowsIndices[ v ];
        vertex_id_t* iend = G->endV + G->rowsIndices[ v + 1 ];
        if( order > 0 )
        {
            std::sort( ibegin, iend, sorter_a( G ) );
        }
        else
        {
            std::sort( ibegin, iend, sorter_d( G ) );
        }
    }

    return fmap;
}

void run( graph_t* G, double* result )
{
    size_t                          n = G->n;
    std::vector<compute_buffer_t>   buffers;
    std::vector<uint32_t>           rows_indices32;
    size_t                          max_work_threads = omp_get_max_threads();
    std::vector< vertex_id_t>       map( sort_graph( G, +1 ) );

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
    for( size_t s = 0; s < n; ++s )
    {
        double  r = 0;
        for( size_t t = 0; t < max_work_threads; ++t )
        {
            compute_buffer_t&   b( buffers[ t ] );
            r += b.partial_result[s];
        }
        result[map[s]] = r*0.5;
    }

    for( size_t t = 0; t < max_work_threads; ++t )
    {
        compute_buffer_t&   b( buffers[ t ] );
        b.release();
    }
}
