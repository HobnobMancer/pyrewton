#!/bin/bash

# :args $1: absolute path to pyrewton prediction/tools dir

# Installing the CAZyme prediciton tools dbCAN, eCAMI and CUPP
# Requirements for these tools are installed via requirements.txt

cd $1

# install dbCAN
cd dbcan

test -d db || mkdir db
cd db \
    && wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312019.fa.nr && diamond makedb --in CAZyDB.07312019.fa.nr -d CAZy \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt && mv dbCAN-HMMdb-V8.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm \
    && cd ../ && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff
# To check the installtion of dbCAN has worked, navigate to the dbCAN directory and run:
# run_dbcan.py EscheriaColiK12MG1655.fna prok --out_dir output_EscheriaColiK12MG1655

# download eCAMI
cd $1
git clone https://github.com/zhanglabNKU/eCAMI.git
mv eCAMI ecami

# download CUPP
curl -o CUPP_v1.0.14.tar.gz "https://files.dtu.dk/fss/public/link/public/stream/read/CUPP_v1.0.14.tar.gz?linkToken=hLin6ni4p-SWuKfp&itemName=CUPP_program"
tar -xzf CUPP_v1.0.14.tar.gz
mv CUPP_v1.0.14 cupp
# delete old file
rm CUPP_v1.0.14.tar.gz
