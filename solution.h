#ifndef SOLUTION_H
#define SOLUTION_H

#include "defs.h"

#include <limits>
#include <algorithm>
#include <parallel/algorithm>
#include <omp.h>

#define INVALID_DISTANCE 255

typedef uint8_t DIST_TYPE;
typedef uint16_t SCOUNT_TYPE;
typedef double DELTA_TYPE;
typedef double PARTIAL_TYPE;

#define MEMALIGN    32
#ifndef __USE_ISOC11
    inline void* aligned_alloc (size_t __alignment, size_t __size)
    {
        void* mem = NULL;
        posix_memalign( &mem, __alignment, __size );
        return mem;
    }
#endif

class wavefront_t
{
    vertex_id_t     m_n;
    vertex_id_t*    m_p;
    vertex_id_t*    m_front;
    vertex_id_t*    m_back;
private:
	wavefront_t( const wavefront_t& );
	wavefront_t& operator = ( const wavefront_t& );
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
    }
    void resize( vertex_id_t n )
    {
        m_n = n;
        if( m_p != NULL )
        {
            free( m_p );
        }
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
    void release()
    {
        if( m_p != NULL )
        {
            free( m_p );
        }
        m_n = 0;
        m_p = NULL;
        m_front = NULL;
        m_back = NULL;
    }
};

struct compute_buffer_t
{
private:
	compute_buffer_t( const compute_buffer_t& );
	compute_buffer_t& operator = ( const compute_buffer_t& );
public:
    size_t          size;
    size_t          mem_align;
    size_t          max_distance;
    DIST_TYPE*      distance;
    SCOUNT_TYPE*    shortest_count;
    vertex_id_t*    vertex_on_level_count;
    PARTIAL_TYPE*   partial_result;
    DELTA_TYPE*     delta;
    wavefront_t     q;
    wavefront_t     qnext;

    compute_buffer_t() :
        size(0),
        mem_align( MEMALIGN ),
        distance( NULL ),
        shortest_count( NULL ),
        vertex_on_level_count( NULL ),
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
                std::cout << (int)distance[i] << " ";
            }
        }
        std::cout << "} ";
        std::cout << "sc = { ";
        for( size_t i = 0; i != size; ++i )
        {
            std::cout << (int)shortest_count[i] << " ";
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
                std::cout << " distance[i] != other.distance[i] " << "i = " << i << " " << (int)distance[i] << " " << (int)other.distance[i] << std::endl;
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

    void resize( const graph_t* G )
    {
        max_distance = std::min( (vertex_id_t)std::numeric_limits<DIST_TYPE>::max(), G->n );
        size = G->n;
        distance = (DIST_TYPE*)aligned_alloc( mem_align, sizeof( DIST_TYPE )*G->n );
        shortest_count = (SCOUNT_TYPE*)aligned_alloc( mem_align, sizeof( SCOUNT_TYPE )*G->n );
        vertex_on_level_count = (vertex_id_t*)aligned_alloc( mem_align, sizeof( vertex_id_t )*max_distance );
        partial_result = (PARTIAL_TYPE*)aligned_alloc( mem_align, sizeof( PARTIAL_TYPE )*G->n );
        delta = (DELTA_TYPE*)aligned_alloc( mem_align, sizeof( DELTA_TYPE )*G->n );

        q.resize( G->n );
        qnext.resize( G->n );

        std::fill( partial_result, partial_result + G->n, 0 );
    }

    void release()
    {
        free( distance );
        free( shortest_count );
        free( vertex_on_level_count );
        free( partial_result );
        free( delta );

		q.release();
		qnext.release();
    }
};

struct graph_coo_t
{
    struct edge_t
    {
        vertex_id_t s,e;
        edge_t( vertex_id_t _s, vertex_id_t _e ) : s( _s ), e( _e ){}
    };
    edge_t*         edges;
    vertex_id_t     size;
    graph_coo_t() : edges( NULL ), size( 0 )
    {
    }
    ~graph_coo_t()
    {
        release();
    }
    void resize( vertex_id_t size )
    {
        release();
        edges = (edge_t*)aligned_alloc( MEMALIGN, sizeof( edge_t )*size );
    }
    void release()
    {
        free( edges );
        edges = NULL;
    }
    void convert( const graph_t* G )
    {
        resize( G->m );
        vertex_id_t n = 0;
        for( vertex_id_t v = 0; v != G->n; ++v )
        {
            vertex_id_t* ibegin = G->endV + G->rowsIndices[ v ];
            vertex_id_t* iend = G->endV + G->rowsIndices[ v + 1 ];

            for( vertex_id_t* e = ibegin; e != iend; ++e )
            {
                if( v <= *e )//Only graph half
                {
                    edges[n] = edge_t( v, *e );
                    ++n;
                }
            }
        }
        size = n;
    }

    static size_t distance( vertex_id_t v1, vertex_id_t v2 )
    {
        if( v1 > v2 )
            return v1 - v2;
        else
            return v2 - v1;
    }

    static size_t distance( const edge_t& e1, const edge_t& e2 )
    {
        return std::min( distance( e1.s, e2.s ), distance( e1.s, e2.e ) ) +
               std::min( distance( e1.e, e2.s ), distance( e1.e, e2.e ) );
    }

    struct sorter
    {
        edge_t root;
        sorter( edge_t r ) : root( r ) {}
        bool operator () ( const edge_t& e1, const edge_t& e2 )
        {
            return distance( e1, root ) < distance( e2, root );
        }
    };

    double total()
    {
        size_t  total_dist = 0;
        for( size_t e = 0; e != size - 1; ++e )
        {
            total_dist += distance( edges[e], edges[e+1] );
        }
        return (double)total_dist/(double)size;
    }

    void reorder()
    {
        std::cout << std::endl;
        std::cout << "Initial mean distance = " << total() << std::endl;

        size_t  block_size = 8*1024;
        edge_t  root( edges[0] );
        for( size_t n = 0; n < size; n += block_size )
        {
            size_t  ibegin = n;
            size_t  iend = std::min( n + block_size, (size_t)size );
            __gnu_parallel::partial_sort( edges + ibegin, edges + iend, edges + size, sorter( root ) );
            root = edges[iend - 1];
        }
        std::cout << "Reordered mean distance = " << total() << std::endl;
    }
};

void simplified_dijkstra( const graph_t* G, const uint32_t* row_indites, vertex_id_t start, DIST_TYPE* distance, SCOUNT_TYPE* shortest_count, wavefront_t& queue );
void bfs( const graph_t* G, const uint32_t* row_indites, vertex_id_t start,
          DIST_TYPE* distance, SCOUNT_TYPE* shortest_count,
          wavefront_t& q, wavefront_t& qnext,
          vertex_id_t* vertex_on_level_count,
          DIST_TYPE& max_distance );
void betweenness_centrality( graph_t* G, const uint32_t* row_indites, vertex_id_t s,
                             const DIST_TYPE* distance,
                             const SCOUNT_TYPE* shortest_count,
                             vertex_id_t* vertex_on_level_count, DIST_TYPE max_distance,
                             DELTA_TYPE* delta,
                             PARTIAL_TYPE* result );

void bfs( const graph_t* G, const graph_coo_t* C, vertex_id_t start,
          DIST_TYPE* distance, SCOUNT_TYPE* shortest_count,
          wavefront_t& q, wavefront_t& qnext,
          vertex_id_t* vertex_on_level_count,
          const double* global_vertex_on_level_count,
          const double* global_unmarked_vertex_count, DIST_TYPE& max_distance );
void betweenness_centrality( const graph_t* G, const graph_coo_t* C, vertex_id_t s,
                             const DIST_TYPE* distance,
                             const SCOUNT_TYPE* shortest_count,
                             vertex_id_t* vertex_on_level_count, DIST_TYPE max_distance,
                             DELTA_TYPE* delta,
                             PARTIAL_TYPE* result );

std::vector< vertex_id_t> sort_graph( const graph_t * const G,  graph_t* Gwork, int order );

#endif // SOLUTION_H
