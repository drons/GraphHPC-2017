#include "defs.h"

#include <limits>
#include <deque>
#include <map>
#include <algorithm>
#include <omp.h>

typedef vertex_id_t DIST_TYPE;
typedef vertex_id_t PARENT_TYPE;

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

void simplified_dijkstra( const graph_t* G, const uint32_t* row_indites, vertex_id_t start, DIST_TYPE* distance, vertex_id_t* shortest_count, wavefront_t& queue )
{
    static const DIST_TYPE      INVALID_DISTANCE = std::numeric_limits<DIST_TYPE>::max();

    memset( distance, INVALID_DISTANCE, sizeof(DIST_TYPE)*G->n );
    memset( shortest_count, 0, sizeof(vertex_id_t)*G->n );

    distance[ start ] = 0;
    shortest_count[ start ] = 1;

    queue.push_back( start );

    while( !queue.empty() )
    {
        vertex_id_t     v = queue.front();
        DIST_TYPE       dist_v = distance[ v ];
        vertex_id_t     shortest_count_v( shortest_count[v] );

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

void betweenness_centrality( graph_t* G, const uint32_t* row_indites, vertex_id_t s,
                             const DIST_TYPE* distance,
                             const vertex_id_t* shortest_count,
                             const wavefront_t& wavefront,
                             double* delta,
                             double* result )
{
    memset( delta, 0, sizeof(double)*G->n );

    const vertex_id_t*  wf_rend( wavefront.rend() );
    for( const vertex_id_t* ii = wavefront.rbegin(); ii != wf_rend; --ii )
    {
        vertex_id_t  w = *ii;
        vertex_id_t* ibegin = G->endV + row_indites[ w ];
        vertex_id_t* iend = G->endV + row_indites[ w + 1 ];
        DIST_TYPE    dist_w_minus_one( distance[w] - 1 );
        double       delta_w( delta[w] );
        double       sc_w( ((double)shortest_count[w]) );

        for( vertex_id_t* e = ibegin; e != iend; ++e )
        {
            vertex_id_t v( *e );
            if( dist_w_minus_one == distance[v] )
            {
                const double sc_v( ((double)shortest_count[v]) );
                delta[v] += (sc_v + sc_v*delta_w)/sc_w;
            }
        }
        if( w != s )
        {
            result[w] += delta_w;
        }
    }
}

struct compute_buffer_t
{
    std::vector<DIST_TYPE>      distance;
    std::vector<vertex_id_t>    shortest_count;
    wavefront_t                 wavefront;
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

        b.distance.resize( G->n );
        b.shortest_count.resize( G->n );
        b.wavefront.resize( G->n );

        b.partial_result.resize( G->n );
        b.delta.resize( G->n );

        std::fill( b.partial_result.begin(), b.partial_result.end(), 0 );
    }

//    std::cout << "omp_get_num_procs   " << omp_get_num_procs() << std::endl;
//    std::cout << "omp_get_num_threads " << omp_get_num_threads() << std::endl;
//    std::cout << "omp_get_max_threads " << omp_get_max_threads() << std::endl;
//    std::cout << "omp_get_num_teams   " << omp_get_num_teams() << std::endl;
//    std::cout << "omp_get_team_num    " << omp_get_team_num() << std::endl;

    #pragma omp parallel for
    for( vertex_id_t s = 0; s < n; ++s )
    {
        compute_buffer_t&   b( buffers[ omp_get_thread_num() ] );
        b.wavefront.reset();
        simplified_dijkstra( G, rows_indices32.data(), s, b.distance.data(), b.shortest_count.data(), b.wavefront );
        betweenness_centrality( G, rows_indices32.data(), s, b.distance.data(), b.shortest_count.data(), b.wavefront, b.delta.data(), b.partial_result.data() );
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
}
