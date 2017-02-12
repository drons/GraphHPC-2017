#ifndef SOLUTION_H
#define SOLUTION_H

#include "defs.h"

#include <limits>
#include <algorithm>
#include <omp.h>

typedef uint8_t DIST_TYPE;
typedef uint16_t SCOUNT_TYPE;
typedef double DELTA_TYPE;
typedef double PARTIAL_TYPE;

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
    void swap( wavefront_t& other )
    {
        std::swap( m_n, other.m_n );
        std::swap( m_p, other.m_p );
        std::swap( m_front, other.m_front );
        std::swap( m_back, other.m_back );
    }
    size_t size() const
    {
        return m_back - m_front;
    }
};

struct compute_buffer_t
{
    std::vector<DIST_TYPE>      distance;
    std::vector<SCOUNT_TYPE>    shortest_count;
    std::vector<vertex_id_t>    vertex_on_level_count;
    std::vector<double>         global_vertex_on_level_count;
    std::vector<double>         global_unmarked_vertex_count;
    wavefront_t                 q;
    wavefront_t                 qnext;
    std::vector<PARTIAL_TYPE>   partial_result;
    std::vector<DELTA_TYPE>     delta;

    void dump_bfs_result( vertex_id_t s, vertex_id_t max_distance )
    {
        std::cout << " s = " << s << " ";

        std::cout << "d  = { ";
        for( size_t i = 0; i != distance.size(); ++i )
        {
            std::cout << distance[i] << " ";
        }
        std::cout << "} ";
        std::cout << "sc = { ";
        for( size_t i = 0; i != shortest_count.size(); ++i )
        {
            std::cout << shortest_count[i] << " ";
        }
        std::cout << "} ";
        std::cout << "vc = { ";
        for( size_t i = 0; i != max_distance + 1; ++i )
        {
            std::cout << vertex_on_level_count[i] << " ";
        }
        std::cout << "}" << std::endl;
    }
    void resize( const graph_t* G )
    {
        size_t              max_distance( std::min( (vertex_id_t)std::numeric_limits<DIST_TYPE>::max(), G->n ) );
        distance.resize( G->n );
        shortest_count.resize( G->n );
        vertex_on_level_count.resize( max_distance );
        global_vertex_on_level_count.resize( max_distance );
        global_unmarked_vertex_count.resize( max_distance );
        q.resize( G->n );
        qnext.resize( G->n );
        partial_result.resize( G->n );
        delta.resize( G->n );

        std::fill( global_vertex_on_level_count.begin(), global_vertex_on_level_count.end(), 0 );
        std::fill( global_unmarked_vertex_count.begin(), global_unmarked_vertex_count.end(), G->n );
        std::fill( partial_result.begin(), partial_result.end(), 0 );
    }
};

void simplified_dijkstra( const graph_t* G, const uint32_t* row_indites, vertex_id_t start, DIST_TYPE* distance, SCOUNT_TYPE* shortest_count, wavefront_t& queue );
void bfs( const graph_t* G, const uint32_t* row_indites, vertex_id_t start,
          DIST_TYPE* distance, SCOUNT_TYPE* shortest_count,
          wavefront_t& q, wavefront_t& qnext,
          vertex_id_t* vertex_on_level_count,
          const double* global_vertex_on_level_count,
          const double* global_unmarked_vertex_count, DIST_TYPE& max_distance );
void betweenness_centrality( graph_t* G, const uint32_t* row_indites, vertex_id_t s,
                             const DIST_TYPE* distance,
                             const SCOUNT_TYPE* shortest_count,
                             vertex_id_t* vertex_on_level_count, DIST_TYPE max_distance,
                             DELTA_TYPE* delta,
                             PARTIAL_TYPE* result );

#endif // SOLUTION_H
