#!/bin/sh
GSIZE=12
#export OMP_NUM_THREADS=56
#export KMP_AFFINITY=scatter,granularity=thread

#./solution_mpi.bin -nIters 2 -in /home/sudorgin/g/rmat-3 &&
#./solution_mpi.bin -nIters 2 -in /home/sudorgin/g/rmat-4 &&
#./solution_mpi.bin -nIters 2 -in /home/sudorgin/g/rmat-5 &&
#./solution_mpi.bin -nIters 2 -in /home/sudorgin/g/rmat-6 &&
#./solution_mpi.bin -nIters 2 -in /home/sudorgin/g/rmat-7 &&
#./gen_RMAT -s $GSIZE &&
./solution_mpi.bin -nIters 2 -in /home/sudorgin/g/rmat-$GSIZE &&
#./gen_valid_info -in rmat-$GSIZE &&
#./validation -ans /home/sudorgin/g/rmat-$GSIZE.ans -res /home/sudorgin/g/rmat-$GSIZE.res &&

#./gen_random -s $GSIZE &&
#./ssolution_mpi.bin -nIters 2 -in /home/sudorgin/g/random-$GSIZE &&
#./gen_valid_info -in random-$GSIZE &&
#./validation -ans random-$GSIZE.ans -res random-$GSIZE.res &&
exit 0
