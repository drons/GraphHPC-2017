#include "defs.h"

#include <limits>
#include <deque>
#include <map>
#include <algorithm>
#include <omp.h>

typedef uint8_t     DIST_TYPE;
typedef vertex_id_t PARENT_TYPE;
typedef uint16_t    SCOUNT_TYPE;
#define UNROLL      2

class wavefront_t
{
    vertex_id_t     m_n;
    vertex_id_t*    m_p;
    vertex_id_t*    m_front;
    vertex_id_t*    m_back;
public:
    wavefront_t()
    {
        m_n = 0;
        m_p = NULL;
        m_front = NULL;
        m_back = NULL;
    }
    ~wavefront_t()
    {
        free( m_p );
    }
    void resize( vertex_id_t n )
    {
        m_n = n;
        if( m_p != NULL )
            free( m_p );
        m_p = (vertex_id_t*)malloc( n*sizeof( vertex_id_t ) );
        m_front = m_p;
        m_back = m_p;
    }
    void push_back( vertex_id_t v )
    {
        *m_back = v;
        ++m_back;
    }
    bool empty() const
    {
        return m_front == m_back;
    }
    vertex_id_t front() const
    {
        return *m_front;
    }
    void pop_front()
    {
        ++m_front;
    }
    void reset()
    {
        m_front = m_p;
        m_back = m_p;
    }
    const vertex_id_t* rbegin() const
    {
        return m_back - 1;
    }
    const vertex_id_t* rend() const
    {
        return m_p - 1;
    }
};

void simplified_dijkstra( const graph_t* G, const uint32_t* row_indites, vertex_id_t start, DIST_TYPE* distance, SCOUNT_TYPE* shortest_count, wavefront_t& queue )
{
    static const DIST_TYPE      INVALID_DISTANCE = std::numeric_limits<DIST_TYPE>::max();

    memset( distance, INVALID_DISTANCE, sizeof(DIST_TYPE)*G->n );
    memset( shortest_count, 0, sizeof(SCOUNT_TYPE)*G->n );

    //std::cout << start << " ";
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

void bfs( const graph_t* G, const uint32_t* row_indites, vertex_id_t start, DIST_TYPE* distance, SCOUNT_TYPE* shortest_count, DIST_TYPE& max_distance )
{
    static const DIST_TYPE      INVALID_DISTANCE = std::numeric_limits<DIST_TYPE>::max();

    memset( distance, INVALID_DISTANCE, sizeof(DIST_TYPE)*G->n*UNROLL );
    memset( shortest_count, 0, sizeof(SCOUNT_TYPE)*G->n*UNROLL );

    for( size_t u = 0; u != UNROLL; ++u )
    {
        distance[ UNROLL*(start+u) + u ] = 0;
        shortest_count[ UNROLL*(start+u) + u ] = 1;
    }

    size_t      num_verts_on_level = 1;
    DIST_TYPE   current_level = 0;
    DIST_TYPE   next_level = 1;

    while( num_verts_on_level > 0 )
    {
        num_verts_on_level = 0;
        for( vertex_id_t v = 0; v != G->n; ++v )
        {
            DIST_TYPE   distance_v[UNROLL];
            SCOUNT_TYPE sc_v[UNROLL];

            bool    skip = true;
            for( size_t u = 0; u != UNROLL; ++u )
            {
                distance_v[u] = distance[ UNROLL*v + u ];
                sc_v[u] = shortest_count[ UNROLL*v + u ];
                if( distance_v[u] == current_level )
                    skip = false;
            }
            if( skip )
            {
                continue;
            }
            {
                vertex_id_t* ibegin = G->endV + row_indites[ v ];
                vertex_id_t* iend = G->endV + row_indites[ v + 1 ];

                for( vertex_id_t* e = ibegin; e != iend; ++e )
                {
                    vertex_id_t w( *e );
                    DIST_TYPE   distance_w[UNROLL];
                    DIST_TYPE   sc_w[UNROLL];

                    for( size_t u = 0; u != UNROLL; ++u )
                    {
                        distance_w[u] = distance[ UNROLL*w + u ];
                        sc_w[u] = shortest_count[ UNROLL*w + u ];
                    }

                    for( size_t u = 0; u != UNROLL; ++u )
                    {
                        if( distance_v[u] == current_level )
                        {
                            if( distance_w[u] == INVALID_DISTANCE )
                            {
                                distance_w[u] = next_level;
                                ++num_verts_on_level;
                                sc_w[u] = sc_v[u];
                            }
                            else
                            if( distance_w[u] == next_level )
                            {
                                sc_w[u] += sc_v[u];
                            }
                        }
                    }

                    if( w == v )
                    {
                        for( size_t u = 0; u != UNROLL; ++u )
                        {
                            distance_v[u] = distance_w[u];
                            sc_v[u] = sc_w[u];
                        }
                    }

                    for( size_t u = 0; u != UNROLL; ++u )
                    {
                        distance[ UNROLL*w + u ] = distance_w[u];
                        shortest_count[ UNROLL*w + u ] = sc_w[u];
                    }
                }
            }
        }
        ++current_level;
        ++next_level;
    }
    max_distance = current_level;
}

void betweenness_centrality( graph_t* G, const uint32_t* row_indites, vertex_id_t s,
                             const DIST_TYPE* distance,
                             const SCOUNT_TYPE* shortest_count,
                             DIST_TYPE max_distance,
                             double* delta,
                             double* result )
{
    memset( delta, 0, sizeof(double)*G->n*UNROLL );

    --max_distance;
    for(;;)
    {
        for( vertex_id_t w = 0; w != G->n; ++w )
        {
            DIST_TYPE    dist_w_minus_one[UNROLL];
            double       delta_w[UNROLL];
            double       sc_w[UNROLL];

            for( size_t u = 0; u != UNROLL; ++u )
            {
                dist_w_minus_one[u] = distance[ UNROLL*w + u ] - 1;
                delta_w[u] = delta[ UNROLL*w + u ];
                sc_w[u] = shortest_count[ UNROLL*w + u ];
            }
            bool    skip = true;
            for( size_t u = 0; u != UNROLL; ++u )
            {
                if( dist_w_minus_one[u] == max_distance )
                    skip = false;
            }
            if( skip )
            {
                continue;
            }
            vertex_id_t* ibegin = G->endV + row_indites[ w ];
            vertex_id_t* iend = G->endV + row_indites[ w + 1 ];

            for( vertex_id_t* e = ibegin; e != iend; ++e )
            {
                vertex_id_t     v( *e );
                DIST_TYPE       distance_v[UNROLL];
                double          sc_v[UNROLL];
                double          delta_v[UNROLL];

                for( size_t u = 0; u != UNROLL; ++u )
                {
                    distance_v[u] = distance[ UNROLL*v + u ];
                    sc_v[u] = shortest_count[ UNROLL*v + u ];
                    delta_v[u] = delta[ UNROLL*v + u ];
                }

                for( size_t u = 0; u != UNROLL; ++u )
                {
                    if( dist_w_minus_one[u] == max_distance )
                    {
                        if( dist_w_minus_one[u] == distance_v[u] )
                        {
                            delta_v[u] += (sc_v[u] + sc_v[u]*delta_w[u])/sc_w[u];
                        }
                    }
                }

                for( size_t u = 0; u != UNROLL; ++u )
                {
                    delta[ UNROLL*v + u ] = delta_v[u];
                }
            }
            for( size_t u = 0; u != UNROLL; ++u )
            {
                if( dist_w_minus_one[u] == max_distance )
                {
                    if( w != s + u )
                    {
                        result[UNROLL*w + u] += delta_w[u];
                    }
                }
            }
        }
        if( max_distance == 0 )
        {
            break;
        }
        --max_distance;
    }
}

struct compute_buffer_t
{
    std::vector<DIST_TYPE>      distance;
    std::vector<SCOUNT_TYPE>    shortest_count;
    std::vector<double>         partial_result;
    std::vector<double>         delta;
};

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

        b.distance.resize( G->n*UNROLL );
        b.shortest_count.resize( G->n*UNROLL );
        b.partial_result.resize( G->n*UNROLL );
        b.delta.resize( G->n*UNROLL );

        std::fill( b.partial_result.begin(), b.partial_result.end(), 0 );
    }

//    std::cout << "omp_get_num_procs   " << omp_get_num_procs() << std::endl;
//    std::cout << "omp_get_num_threads " << omp_get_num_threads() << std::endl;
//    std::cout << "omp_get_max_threads " << omp_get_max_threads() << std::endl;
//    std::cout << "omp_get_num_teams   " << omp_get_num_teams() << std::endl;
//    std::cout << "omp_get_team_num    " << omp_get_team_num() << std::endl;

    compute_buffer_t    t;
    wavefront_t         w;
    t.distance.resize( G->n );
    t.shortest_count.resize( G->n );
    t.partial_result.resize( G->n );
    t.delta.resize( G->n );
    w.resize( G->n );
    std::fill( t.partial_result.begin(), t.partial_result.end(), 0 );

    #pragma omp parallel for schedule ( guided )
    for( vertex_id_t s = 0; s < n/UNROLL; ++s )
    {
        compute_buffer_t&   b( buffers[ omp_get_thread_num() ] );
        DIST_TYPE           max_distance = 0;

        bfs( G, rows_indices32.data(), s*UNROLL, b.distance.data(), b.shortest_count.data(), max_distance );
        betweenness_centrality( G, rows_indices32.data(), s*UNROLL, b.distance.data(), b.shortest_count.data(), max_distance, b.delta.data(), b.partial_result.data() );
    }

    #pragma omp parallel for
    for( vertex_id_t s = 0; s < n; ++s )
    {
        double  r = 0;
        for( size_t t = 0; t < max_work_threads; ++t )
        {
            compute_buffer_t&   b( buffers[ t ] );
            for( size_t u = 0; u != UNROLL; ++u )
            {
                r += b.partial_result[ UNROLL*s + u ];
            }
        }
        result[s] = r*0.5;
    }
}
