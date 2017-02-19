#!/bin/sh
GSIZE=17
#export OMP_NUM_THREADS=56
#export KMP_AFFINITY=scatter,granularity=thread

for GS in 4 5 6 7 8 9 10 11 12; do
	./gen_RMAT -s $GS &&
	./solution -nIters 2 -in rmat-$GS &&
	./gen_valid_info -in rmat-$GS &&
	./validation -ans rmat-$GS.ans -res rmat-$GS.res &&
	./gen_random -s $GS &&
	./solution -nIters 2 -in random-$GS &&
	./gen_valid_info -in random-$GS &&
	./validation -ans random-$GS.ans -res random-$GS.res ||
	exit 1
done

./gen_RMAT -s $GSIZE &&
./solution -nIters 2 -in rmat-$GSIZE &&
#./gen_valid_info -in rmat-$GSIZE &&
#./validation -ans rmat-$GSIZE.ans -res rmat-$GSIZE.res &&

./gen_random -s $GSIZE &&
./solution -nIters 2 -in random-$GSIZE &&
#./gen_valid_info -in random-$GSIZE &&
#./validation -ans random-$GSIZE.ans -res random-$GSIZE.res &&
exit 0
