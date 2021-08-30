#!/bin/bash
tar -xzf R361.tar.gz
tar -xzf packages.tar.gz
tar -xzf SLIBS.tar.gz

export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
export R_LIBS=$PWD/packages
export LD_LIBRARY_PATH=$PWD/SS:$LD_LIBRARY_PATH

Rscript getRslt_sim.R $1 $2 $3 $4 $5
