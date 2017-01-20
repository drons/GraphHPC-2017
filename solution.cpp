#include "defs.h"

#include <limits>
#include <deque>
#include <map>
#include <algorithm>
#include <omp.h>

template< typename T >
class buckets_t
{
public:
    typedef	std::vector<T>                  bucket_t;
    typedef	std::vector< bucket_t >         page_t;
    typedef size_t                          size_type;
    typedef typename bucket_t::size_type    bsize_type;
    enum { PAGE_SIZE = 1024 };
private:
    std::deque<page_t>              m_pages;
    size_type                       m_offset;
    typename page_t::size_type      m_L;
public:
    buckets_t()
    {
        m_offset = 0;
        m_L = 0;
    }
    ~buckets_t()
    {
    }

    bucket_t* front()
    {
        while( m_pages.size() > 0 )
        {
            page_t& front_page = m_pages.front();
            for( ; m_L != PAGE_SIZE; ++m_L )
            {
                if( !front_page[m_L].empty() )
                    return &front_page[m_L];
            }
            m_pages.pop_front();
            m_L = 0;
            m_offset += PAGE_SIZE;
        }
        return NULL;
    }

    bucket_t& get_bucket( size_type N ) //buckets_on_level[ N ];
    {
        N -= m_offset;
        int page_n = (int)( N/PAGE_SIZE );
        int bucket_n = (int)( N%PAGE_SIZE );
           std::cout << page_n << " ";
        return m_pages[page_n][bucket_n];
    }

    void push( size_type num_of_bucket, const T& point_location )
    {
        num_of_bucket -= m_offset;
        size_type   page_n = num_of_bucket/PAGE_SIZE;
        size_type   bucket_n = num_of_bucket%PAGE_SIZE;
        while( page_n >= m_pages.size() )
        {
            m_pages.emplace_back( page_t( PAGE_SIZE ) );
        }
        m_pages[page_n][bucket_n].push_back( point_location );
    }

    bool empty() const
    {
        return m_pages.empty();
    }
};

typedef vertex_id_t DIST_TYPE;
typedef vertex_id_t PARENT_TYPE;

void dijkstra( graph_t* G, vertex_id_t start, std::vector< DIST_TYPE >& distance, std::vector<vertex_id_t>& shortest_count, std::vector<vertex_id_t>& wavefront )
{
    static const DIST_TYPE      INVALID_DISTANCE = std::numeric_limits<DIST_TYPE>::max();

    distance.resize( G->n );
    shortest_count.resize( G->n );
    wavefront.reserve( G->n );

    std::fill( distance.begin(), distance.end(), INVALID_DISTANCE );
    std::fill( shortest_count.begin(), shortest_count.end(), 0 );

    distance[ start ] = 0;
    shortest_count[ start ] = 1;

    typedef typename buckets_t< vertex_id_t >::bsize_type   bucket_size_t;
    typedef typename buckets_t< vertex_id_t >::size_type    buckets_size_t;
    typedef typename buckets_t< vertex_id_t >::bucket_t     bucket_t;
    buckets_t< vertex_id_t >                                queue;

    queue.push( 0, start );

    while( !queue.empty() )
    {
        bucket_t*       bucketL = queue.front();
        if( bucketL == NULL )
            break;
        vertex_id_t     current = bucketL->back();
        DIST_TYPE       current_dist = distance[ current ];

        wavefront.push_back( current );
        bucketL->pop_back();

        edge_id_t       ibegin = G->rowsIndices[ current ];
        edge_id_t       iend = G->rowsIndices[ current + 1 ];

        for( edge_id_t e = ibegin; e != iend; ++e )
        {
            vertex_id_t next( G->endV[e] );
            DIST_TYPE   next_distance = current_dist + 1;
            DIST_TYPE&  d( distance[ next ] );
            if( next_distance < d )
            {
                // Checked if node(i2,j2) is already in buckets array, it must be removed
                // If node(i2,j2) is already in buckets array, it in array[N]
                if( d != INVALID_DISTANCE )
                {
                    bucket_t&       bucketN = queue.get_bucket( d );
                    bucket_size_t   size = bucketN.size();

                    if( size != 1 )
                    {
                        bucket_size_t n = 0;
                        for( bucket_size_t i = 0; i != size; i++ )
                        {
                            if( bucketN[i] == next )
                                n = i;
                        }
                        bucketN.erase( bucketN.begin() + n );
                    }
                    else
                    {
                        bucketN.clear();
                    }
                }
                d = next_distance;
                shortest_count[ next ] = shortest_count[ current ];
                queue.push( (buckets_size_t)next_distance, next );
            }
            else
            if( next_distance == distance[ next ] )
            {
                shortest_count[ next ] += shortest_count[ current ];
            }
        }
    }
    return;
}

void betweenness_centrality( graph_t* G, vertex_id_t s,
                             const std::vector< DIST_TYPE >& distance,
                             const std::vector<vertex_id_t>& shortest_count,
                             const std::vector<vertex_id_t>& wavefront,
                             double* result )
{
    std::vector<double> delta;

    delta.resize( G->n );
    std::fill( delta.begin(), delta.end(), 0 );

    for( auto ii = wavefront.rbegin(); ii != wavefront.rend(); ++ii )
    {
        vertex_id_t w = *ii;
        edge_id_t   ibegin = G->rowsIndices[ w ];
        edge_id_t   iend = G->rowsIndices[ w + 1 ];

        for( edge_id_t e = ibegin; e != iend; ++e )
        {
            vertex_id_t v( G->endV[e] );
            if( distance[w] == distance[v] + 1 )
            {
                delta[v] += (shortest_count[v] + shortest_count[v]*delta[w])/((double)shortest_count[w]);
            }
        }
        if( w != s )
        {
            result[w] += delta[w];
        }
    }
}

void run( graph_t* G, double* result )
{
    vertex_id_t n = G->n;
    vertex_id_t stride = 16;

    #pragma omp parallel for
    for( vertex_id_t s = 0; s < n; s += stride )
    {
        std::vector< DIST_TYPE >    distance;
        std::vector<vertex_id_t>    shortest_count;
        std::vector<vertex_id_t>    wavefront;
        std::vector<double>         partial_result;

        partial_result.resize( G->n );
        std::fill( partial_result.begin(), partial_result.end(), 0 );

        vertex_id_t last = s + stride;
        last = std::min( last, n );

        for( vertex_id_t t = s; t < last; ++t )
        {
            wavefront.clear();
            dijkstra( G, t, distance, shortest_count, wavefront );
            betweenness_centrality( G, t, distance, shortest_count, wavefront, partial_result.data() );
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
