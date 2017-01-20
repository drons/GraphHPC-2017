#include "defs.h"

char inFilename[FILENAME_LEN];
char resFilename[FILENAME_LEN];

using namespace std;

/* helper */
void usage(int argc, char **argv)
{
    printf("%s\n", argv[0]);
    printf("Usage:\n");
    printf("%s -in <input> -res <result>\n", argv[0]);
    printf("Options:\n");
    printf("    -in <input> -- input graph filename\n");
    printf("    -res <result> -- result filename\n");
    exit(1);
}

/* initialization */
void init(int argc, char **argv)
{
    int l;
    bool no_inFilename = true;
    bool no_resFilename = true;
    inFilename[0] = resFilename[0] = '\0';
    
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-in")) {
            l = strlen(argv[++i]);
            strncpy(inFilename, argv[i], (l > FILENAME_LEN - 1 ? FILENAME_LEN - 1 : l));
            no_inFilename = false;
        }
        
        if (!strcmp(argv[i], "-res")) {
            l = strlen(argv[++i]);
            strncpy(resFilename, argv[i], (l > FILENAME_LEN - 1 ? FILENAME_LEN - 1 : l));
            no_resFilename = false;
        }
    }
    
    if (no_inFilename || no_resFilename) {
        usage(argc, argv);
    }
}

int main(int argc, char **argv)
{
    graph_t g;
    init(argc, argv);
    
    /* read graph from the file */
    readGraph(&g, inFilename);
    
    double *answer = new double[g.n];
    assert(answer != NULL);
    /* get the right answer */
    run(&g, answer);
    
    FILE *f = fopen(resFilename, "r");
    assert(f != NULL);
    double *result = new double[g.n];
    assert(result != NULL);
    /* read the result */
    assert(fread(result, sizeof(double), g.n, f) == g.n);
    fclose(f);
    
    /* comparison the result with right answer */
    for (vertex_id_t i = 0; i < g.n; i++) {
        if (fabs(answer[i] - result[i]) > eps) {
            /* red color */
            cout << "\033[1;31mWrong answer\033[0m\n";
            
            delete[] answer;
            delete[] result;
            freeGraph(&g);
            
            return 0;
        }
    }
    
    /* green color */
    cout << "\033[1;32mAccepted\033[0m\n";
    
    delete[] answer;
    delete[] result;
    freeGraph(&g);
    
    return 0;
}
