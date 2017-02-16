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
    size_t          size;
    size_t          mem_align;
    size_t          max_distance;
    DIST_TYPE*      distance;
    SCOUNT_TYPE*    shortest_count;
    vertex_id_t*    vertex_on_level_count;
    double*         global_vertex_on_level_count;
    double*         global_unmarked_vertex_count;
    PARTIAL_TYPE*   partial_result;
    DELTA_TYPE*     delta;
    wavefront_t     q;
    wavefront_t     qnext;

    compute_buffer_t() :
        size(0),
        mem_align( 16 ),
        distance( NULL ),
        shortest_count( NULL ),
        vertex_on_level_count( NULL ),
        global_vertex_on_level_count( NULL ),
        global_unmarked_vertex_count( NULL ),
        partial_result( NULL ),
        delta( NULL )
    {
    }

    void dump_bfs_result( vertex_id_t s, vertex_id_t max_distance )
    {
        std::cout << " s = " << s << " ";

        std::cout << "d  = { ";
        for( size_t i = 0; i != size; ++i )
        {
            if( distance[i] > max_distance + 1 )
            {
                std::cout << " " << " ";
            }
            else
            {
                std::cout << distance[i] << " ";
            }
        }
        std::cout << "} ";
        std::cout << "sc = { ";
        for( size_t i = 0; i != size; ++i )
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

    bool is_equal( const compute_buffer_t& other )
    {
        for( size_t i = 0; i != size; ++i )
        {
            if( distance[i] != other.distance[i] )
            {
                std::cout << " distance[i] != other.distance[i] " << "i = " << i << " " << distance[i] << " " << other.distance[i] << std::endl;
                return false;
            }
        }
        for( size_t i = 0; i != size; ++i )
        {
            if( shortest_count[i] != other.shortest_count[i] )
            {
                std::cout << " shortest_count[i] != other.shortest_count[i] " << std::endl;
                return false;
            }
        }
        return true;
    }
#ifndef __USE_ISOC11
    void* aligned_alloc (size_t __alignment, size_t __size)
    {
        void* mem = NULL;
        posix_memalign( &mem, __alignment, __size );
        return mem;
    }
#endif
    void resize( const graph_t* G )
    {
        max_distance = std::min( (vertex_id_t)std::numeric_limits<DIST_TYPE>::max(), G->n );
        size = G->n;
        distance = (DIST_TYPE*)aligned_alloc( mem_align, sizeof( DIST_TYPE )*G->n );
        shortest_count = (SCOUNT_TYPE*)aligned_alloc( mem_align, sizeof( SCOUNT_TYPE )*G->n );
        vertex_on_level_count = (vertex_id_t*)aligned_alloc( mem_align, sizeof( vertex_id_t )*max_distance );
        global_vertex_on_level_count = (double*)aligned_alloc( mem_align, sizeof( double )*max_distance );
        global_unmarked_vertex_count = (double*)aligned_alloc( mem_align, sizeof( double )*max_distance );
        partial_result = (PARTIAL_TYPE*)aligned_alloc( mem_align, sizeof( PARTIAL_TYPE )*G->n );
        delta = (DELTA_TYPE*)aligned_alloc( mem_align, sizeof( DELTA_TYPE )*G->n );

        q.resize( G->n );
        qnext.resize( G->n );

        std::fill( global_vertex_on_level_count, global_vertex_on_level_count + max_distance, 0 );
        std::fill( global_unmarked_vertex_count, global_unmarked_vertex_count + max_distance, G->n );
        std::fill( partial_result, partial_result + G->n, 0 );
    }

    void release()
    {
        free( distance );
        free( shortest_count );
        free( vertex_on_level_count );
        free( global_vertex_on_level_count );
        free( global_unmarked_vertex_count );
        free( partial_result );
        free( delta );
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
