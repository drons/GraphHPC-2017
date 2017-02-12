#include "defs.h"
#include <mpi.h>

void run_mpi(graph_t *G, double *result);

using namespace std;

char inFilename[FILENAME_LEN];
char outFilename[FILENAME_LEN];

int nIters = 16;
#if defined(CLOCK_MONOTONIC)
#define CLOCK CLOCK_MONOTONIC
#elif defined(CLOCK_REALTIME)
#define CLOCK CLOCK_REALTIME
#else
#error "Failed to find a timing clock."
#endif

/* helper */
void usage(int argc, char **argv)
{
    printf("Usage:\n");
    printf("%s -in <input> [options]\n", argv[0]);
    printf("Options:\n");
    printf("    -in <input> -- input graph filename\n");
    printf("    -nIters <nIters> -- number of iterations. By default 64\n");
    exit(1);
}

/* initialization */
void init(int argc, char **argv, graph_t *G)
{
    int l;
    inFilename[0] = '\0';
    outFilename[0] = '\0';
    bool no_inFilename = true;
    if (argc == 1) {
        usage(argc, argv);
    }
    
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-in")) {
            l = strlen(argv[++i]);
            strncpy(inFilename, argv[i], (l > FILENAME_LEN - 1 ? FILENAME_LEN - 1 : l));
            no_inFilename = false;
        }
        
        if (!strcmp(argv[i], "-nIters")) {
            nIters = (int)atoi(argv[++i]);
        }
    }
    
    if (no_inFilename) {
        usage(argc, argv);
    }
    
    if (strlen(outFilename) == 0) {
        sprintf(outFilename, "%s.res", inFilename);
    }
}

/* write result to the file */
void write_output_information(double *result, vertex_id_t n, char *filename)
{
    FILE *f = fopen(filename, "wb");
    assert(f != NULL);
    assert(fwrite(result, sizeof(double), n, f) == n);
    fclose(f);
}

int main(int argc, char **argv)
{
    graph_t g;
    struct timespec start_ts, finish_ts;
    double *alg_time;
    /* initializing and reading the graph */
    init(argc, argv, &g);

    int     rc = -1;

    rc = MPI_Init( &argc, &argv );
    if( 0 != rc )
    {
        std::cout << "MPI_Init failed" << endl;
        MPI_Abort( MPI_COMM_WORLD, rc );
    }

    readGraph(&g, inFilename);

    int mpi_size = -1;
    int mpi_rank = -1;
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

    //printGraph(&g);
    
    alg_time = new double[nIters];
    assert(alg_time != NULL);
    double *result = new double[g.n];
    assert(result != NULL);
    
    if( mpi_rank == 0 )
    {
        cout << "start algorithm..." << endl;
    }
    for (int i = 0; i < nIters; i++) {
        if( mpi_rank == 0 )
        {
            cout << "iter number = " << i;
        }
        memset(result, 0, g.n * sizeof(double));
        clock_gettime(CLOCK, &start_ts);
        run_mpi(&g, result);
        clock_gettime(CLOCK, &finish_ts);
        double time = (finish_ts.tv_nsec - (double)start_ts.tv_nsec) * 1.0e-9 + (finish_ts.tv_sec - (double)start_ts.tv_sec);
        alg_time[i] = time;
        if( mpi_rank == 0 )
        {
            cout.precision(5);
            cout << " finished. Time is " << fixed << time << " sec." << endl;
        }
    }
    
    if( mpi_rank == 0 )
    {
        cout << "algorithm iterations finished." << endl;

        write_output_information(result, g.n, outFilename);

        /* final print */
        double min_time, max_time, avg_time;
        max_time = avg_time = 0;
        min_time = DBL_MAX;
        for (int i = 0; i < nIters; ++i) {
            avg_time += alg_time[i];
            if (alg_time[i] < min_time) {
                min_time = alg_time[i];
            }

            if (alg_time[i] > max_time) {
                max_time = alg_time[i];
            }
        }

        avg_time /= nIters;

        cout << inFilename << ": vertices = " << g.n << " edges = " << g.m << " nIters = " << nIters
            << " time min = " << min_time << " avg = " << avg_time << " max = " << max_time << endl;
        cout << "Time = " << avg_time << " sec." << endl;
    }

    delete[] alg_time;
    delete[] result;
    freeGraph(&g);

    MPI_Finalize();
    
    return 0;
}
