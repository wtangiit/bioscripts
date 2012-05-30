#!/bin/sh

#to install hhsuit at Megallan VM

mkdir hhsuit
cd hhsuit

#download source file

wget ftp://toolkit.genzentrum.lmu.de/HH-suite/hhsuite-latest.tar.gz

tar -xzvf hhsuite-latest.tar.gz

cd hhsuite*

#make and install 
make
sudo make install INSTALL_DIR=/usr/local

#set ENV and PATH
export HHLIB=/usr/local/lib/hh
export PATH=$PATH:/usr/local/bin:$HHLIB/scripts

#configure .bashrc to automate configuration for env and path
echo "export HHLIB=/usr/local/lib/hh
PATH=\$PATH:/usr/local/bin:\$HHLIB/scripts
" >> ~/.bashrc

#downloading database

#SCOP database
cd /mnt
sudo wget ftp://toolkit.genzentrum.lmu.de/HH-suite/databases/hhsuite_dbs/scop70_1.75.tar.gz
tar xzvf scop70_1.75.tar.gz
cd -
