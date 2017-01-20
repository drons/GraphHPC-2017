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
    wavefront_t( vertex_id_t n )
    {
        m_n = n;
        m_p = (vertex_id_t*)malloc( n*sizeof( vertex_id_t ) );
        m_front = m_p;
        m_back = m_p;
    }
    ~wavefront_t()
    {
        free( m_p );
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

void simplified_dijkstra( graph_t* G, vertex_id_t start, std::vector< DIST_TYPE >& distance, std::vector<vertex_id_t>& shortest_count, wavefront_t& queue )
{
    static const DIST_TYPE      INVALID_DISTANCE = std::numeric_limits<DIST_TYPE>::max();

    std::fill( distance.begin(), distance.end(), INVALID_DISTANCE );
    std::fill( shortest_count.begin(), shortest_count.end(), 0 );

    distance[ start ] = 0;
    shortest_count[ start ] = 1;

    queue.push_back( start );

    while( !queue.empty() )
    {
        vertex_id_t     v = queue.front();
        DIST_TYPE       dist_v = distance[ v ];
        vertex_id_t     shortest_count_v( shortest_count[v] );

        queue.pop_front();

        edge_id_t       ibegin = G->rowsIndices[ v ];
        edge_id_t       iend = G->rowsIndices[ v + 1 ];

        for( edge_id_t e = ibegin; e != iend; ++e )
        {
            vertex_id_t w( G->endV[e] );
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

void betweenness_centrality( graph_t* G, vertex_id_t s,
                             const std::vector< DIST_TYPE >& distance,
                             const std::vector<vertex_id_t>& shortest_count,
                             const wavefront_t& wavefront,
                             std::vector<double>& delta,
                             double* result )
{
    std::fill( delta.begin(), delta.end(), 0 );

    const vertex_id_t*  wf_rend( wavefront.rend() );
    for( const vertex_id_t* ii = wavefront.rbegin(); ii != wf_rend; --ii )
    {
        vertex_id_t  w = *ii;
        vertex_id_t* ibegin = G->endV + G->rowsIndices[ w ];
        vertex_id_t* iend = G->endV + G->rowsIndices[ w + 1 ];
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

void run( graph_t* G, double* result )
{
    vertex_id_t n = G->n;
    vertex_id_t stride = 16;

//    std::cout << "omp_get_num_procs   " << omp_get_num_procs() << std::endl;
//    std::cout << "omp_get_max_threads " << omp_get_max_threads() << std::endl;

    #pragma omp parallel for
    for( vertex_id_t s = 0; s < n; s += stride )
    {
        std::vector< DIST_TYPE >    distance;
        std::vector<vertex_id_t>    shortest_count;
        wavefront_t                 wavefront( G->n );
        std::vector<double>         partial_result;
        std::vector<double>         delta;

        distance.resize( G->n );
        shortest_count.resize( G->n );

        partial_result.resize( G->n );
        delta.resize( G->n );

        std::fill( partial_result.begin(), partial_result.end(), 0 );

        vertex_id_t last = s + stride;
        last = std::min( last, n );

        for( vertex_id_t t = s; t < last; ++t )
        {
            wavefront.reset();
            simplified_dijkstra( G, t, distance, shortest_count, wavefront );
            betweenness_centrality( G, t, distance, shortest_count, wavefront, delta, partial_result.data() );
        }

        #pragma omp critical
        for( size_t n = 0; n != G->n; ++n )
        {
            result[n] += partial_result[n];
        }
    }

    for( vertex_id_t s = 0; s < n; ++s )
        result[s] *= 0.5;
}
