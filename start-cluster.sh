#!/bin/sh
RET="OK"
KEY_FILES_PATH="${HOME}/ssh_keys"
KEY_FILE="${KEY_FILES_PATH}/my-17-key"
REMOTESERVER=sudorgin@vertical.nicevt.ru
REMOTEPATH=/home/sudorgin/bc

#install -d "${KEY_FILES_PATH}"
#ssh-keygen -t rsa -P "" -f ${KEY_FILE} && \
#ssh-copy-id -i ${KEY_FILE}.pub ${REMOTESERVER} || RET="FAILED"

make clean
ssh -i "${KEY_FILE}" $REMOTESERVER "rm -r $REMOTEPATH"
ssh -i "${KEY_FILE}" $REMOTESERVER "install -d $REMOTEPATH" &&
scp -i "${KEY_FILE}" solution_mpi *.cpp *.h Makefile* *.sh $REMOTESERVER:$REMOTEPATH &&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile -j6" > /dev/null&&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile clean" &&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && rm Makefile && cp Makefile.mpi Makefile" &&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile -j6" &&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && srun --nodes=8 --ntasks-per-node=1 ./run-cluster.sh" &&
echo "DONE!" || echo "FAIL!"
