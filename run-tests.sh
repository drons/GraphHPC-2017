#!/bin/sh

for GS in 4 5 6 7 8 9 10 11 12; do
	./gen_RMAT -s $GS &&
	./solution -nIters 1 -in rmat-$GS &&
	./gen_valid_info -in rmat-$GS &&
	./validation -ans rmat-$GS.ans -res rmat-$GS.res &&
	./gen_random -s $GS &&
	./solution -nIters 1 -in random-$GS &&
	./gen_valid_info -in random-$GS &&
	./validation -ans random-$GS.ans -res random-$GS.res ||
	exit 1
done

exit 0
