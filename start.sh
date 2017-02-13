#!/bin/sh
GSIZE=4
make -f Makefile.gcc clean &&
make -f Makefile.gcc &&
./gen_RMAT -s $GSIZE &&
./solution -nIters 1 -in rmat-$GSIZE &&
./gen_valid_info -in rmat-$GSIZE &&
./validation -ans rmat-$GSIZE.ans -res rmat-$GSIZE.res &&

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
scp -i "${KEY_FILE}" *.cpp *.h Makefile Makefile.prof-gen Makefile.mpi contest-run.sh contest-analize.sh $REMOTESERVER:$REMOTEPATH &&
#scp -i "${KEY_FILE}" graph_tools.cpp main.cpp solution.cpp *.h Makefile contest-run.sh $REMOTESERVER:$REMOTEPATH &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile.prof-gen" &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && srun ./contest-analize.sh" &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile.prof-gen clean" &&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile" &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile clean" &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile.mpi" &&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && srun ./contest-run.sh" &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && srun --nodes=24 --ntasks-per-node=1 ./contest-run.sh" &&
echo "DONE!" || echo "FAIL!"
#ssh -i "${KEY_FILE}" $REMOTESERVER "rm -r $REMOTEPATH"
