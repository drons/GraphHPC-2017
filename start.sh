#!/bin/sh
make -f Makefile.gcc &&
./gen_RMAT -s 3 &&
./solution -nIters 1 -in rmat-3 &&
./validation -in rmat-3 -res rmat-3.res &&

#./gen_RMAT -s 4 &&
#./solution -nIters 1 -in rmat-4 &&
#./validation -in rmat-4 -res rmat-4.res &&

#./gen_RMAT -s 5 &&
#./solution -nIters 1 -in rmat-5 &&
#./validation -in rmat-5 -res rmat-5.res &&

#./gen_RMAT -s 6 &&
#./solution -nIters 1 -in rmat-6 &&
#./validation -in rmat-6 -res rmat-6.res &&

./gen_RMAT -s 7 &&
./solution -nIters 1 -in rmat-7 &&
./validation -in rmat-7 -res rmat-7.res &&

#./gen_RMAT -s 8 &&
#./solution -nIters 1 -in rmat-8 &&
#./validation -in rmat-8 -res rmat-8.res &&

#./gen_RMAT -s 9 &&
#./solution -nIters 1 -in rmat-9 &&
#./validation -in rmat-9 -res rmat-9.res &&

./gen_RMAT -s 10 &&
./solution -nIters 1 -in rmat-10 &&
./validation -in rmat-10 -res rmat-10.res &&

#./gen_RMAT -s 11 &&
#./solution -nIters 1 -in rmat-11 &&
#./validation -in rmat-11 -res rmat-11.res &&

#./gen_RMAT -s 12 &&
#./solution -nIters 1 -in rmat-12 &&
#./validation -in rmat-12 -res rmat-12.res &&

#./gen_RMAT -s 16 &&
#./solution -nIters 1 -in rmat-16 &&
#./validation -in rmat-16 -res rmat-16.res &&
echo "DONE!"

exit 1

RET="OK"
KEY_FILES_PATH="${HOME}/ssh_keys"
KEY_FILE="${KEY_FILES_PATH}/my-17-key"
REMOTESERVER=sudorgin@access-node.nicevt.ru
REMOTEPATH=/home/sudorgin/bc

#install -d "${KEY_FILES_PATH}"
#ssh-keygen -t rsa -P "" -f ${KEY_FILE} && \
#ssh-copy-id -i ${KEY_FILE}.pub ${REMOTESERVER} || RET="FAILED"

ssh -i "${KEY_FILE}" $REMOTESERVER "rm -r $REMOTEPATH"
ssh -i "${KEY_FILE}" $REMOTESERVER "install -d $REMOTEPATH" &&
scp -i "${KEY_FILE}" *.cpp *.h Makefile Makefile.prof-gen contest-run.sh contest-analize.sh $REMOTESERVER:$REMOTEPATH &&
#scp -i "${KEY_FILE}" graph_tools.cpp main.cpp solution.cpp *.h Makefile contest-run.sh $REMOTESERVER:$REMOTEPATH &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile.prof-gen" &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && srun ./contest-analize.sh" &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile.prof-gen clean" &&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile" &&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && srun ./contest-run.sh" &&
echo "DONE!" || echo "FAIL!"
#ssh -i "${KEY_FILE}" $REMOTESERVER "rm -r $REMOTEPATH"
