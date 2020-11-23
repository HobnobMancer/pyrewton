# !/bin/bash

# Installing pyrewton and all requirements
# :args $1: name of virutal environment

# build conda virtual environment and install requirements via conda forge
conda create -n $1 python=3.8 diamond hmmer prodical -c conda-forge -c bioconda
conda activate $1

# install all other requirements within requirements.txt
pip3 install -r requirements.txt

# install pyrewton as execurable
pip3 install -e .

# install dbCAN
pip install run-dbcan==2.0.11
cd pyrewton/cazymes/prediction/tools/dbcan
test -d db || mkdir
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

# install eCAMI
cd ..
git clone https://github.com/zhanglabNKU/eCAMI.git
mv eCAMI ecami
# requirements for eCAMI were installed via requirements.txt

# install CUPP
cd cupp
# go to this address: https://files.dtu.dk/userportal/#/shared/public/hLin6ni4p-SWuKfp/CUPP_program/CUPP_v1.0.14
# select 'Download Entire Folder'
# # it's just the how part now