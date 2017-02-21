#!/bin/sh
RET="OK"
KEY_FILES_PATH="${HOME}/ssh_keys"
KEY_FILE="${KEY_FILES_PATH}/my-17-key"
REMOTESERVER=sudorgin@access-node.nicevt.ru
REMOTEPATH=/home/sudorgin/bc

#install -d "${KEY_FILES_PATH}"
#ssh-keygen -t rsa -P "" -f ${KEY_FILE} && \
#ssh-copy-id -i ${KEY_FILE}.pub ${REMOTESERVER} || RET="FAILED"

make clean
ssh -i "${KEY_FILE}" $REMOTESERVER "rm -r $REMOTEPATH"
ssh -i "${KEY_FILE}" $REMOTESERVER "install -d $REMOTEPATH" &&
scp -i "${KEY_FILE}" *.cpp *.h Makefile* *.sh $REMOTESERVER:$REMOTEPATH &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile.prof-gen -j6" &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && ./gen_RMAT -s 17 " &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && ./gen_random -s 17 " &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && srun ./contest-analize-rmat.sh" &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && srun ./contest-analize-random.sh" &&
#ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile.prof-gen clean" &&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && make -f Makefile -j6"&&
ssh -i "${KEY_FILE}" $REMOTESERVER "cd $REMOTEPATH && srun ./run-single.sh" &&
echo "DONE!" || echo "FAIL!"
#ssh -i "${KEY_FILE}" $REMOTESERVER "rm -r $REMOTEPATH"
