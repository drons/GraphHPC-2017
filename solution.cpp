#include "defs.h"

#include <limits>
#include <deque>
#include <map>
#include <algorithm>
#include <omp.h>

typedef vertex_id_t DIST_TYPE;
typedef vertex_id_t PARENT_TYPE;

void simplified_dijkstra( graph_t* G, vertex_id_t start, std::vector< DIST_TYPE >& distance, std::vector<vertex_id_t>& shortest_count, std::vector<vertex_id_t>& wavefront )
{
    static const DIST_TYPE      INVALID_DISTANCE = std::numeric_limits<DIST_TYPE>::max();

    std::fill( distance.begin(), distance.end(), INVALID_DISTANCE );
    std::fill( shortest_count.begin(), shortest_count.end(), 0 );

    distance[ start ] = 0;
    shortest_count[ start ] = 1;

    std::deque< vertex_id_t > queue;

    queue.push_back( start );

    while( !queue.empty() )
    {
        vertex_id_t     v = queue.front();
        DIST_TYPE       dist_v = distance[ v ];
        vertex_id_t     shortest_count_v( shortest_count[v] );

        wavefront.push_back( v );
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
                             const std::vector<vertex_id_t>& wavefront,
                             std::vector<double>& delta,
                             double* result )
{
    std::fill( delta.begin(), delta.end(), 0 );

    for( auto ii = wavefront.rbegin(); ii != wavefront.rend(); ++ii )
    {
        vertex_id_t w = *ii;
        edge_id_t   ibegin = G->rowsIndices[ w ];
        edge_id_t   iend = G->rowsIndices[ w + 1 ];
        DIST_TYPE   dist_w_minus_one( distance[w] - 1 );
        double      delta_w( delta[w] );
        double      sc_w( ((double)shortest_count[w]) );

        for( edge_id_t e = ibegin; e != iend; ++e )
        {
            vertex_id_t v( G->endV[e] );
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
        std::vector<vertex_id_t>    wavefront;
        std::vector<double>         partial_result;
        std::vector<double>         delta;

        distance.resize( G->n );
        shortest_count.resize( G->n );
        wavefront.reserve( G->n );

        partial_result.resize( G->n );
        delta.resize( G->n );

        std::fill( partial_result.begin(), partial_result.end(), 0 );

        vertex_id_t last = s + stride;
        last = std::min( last, n );

        for( vertex_id_t t = s; t < last; ++t )
        {
            wavefront.clear();
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
