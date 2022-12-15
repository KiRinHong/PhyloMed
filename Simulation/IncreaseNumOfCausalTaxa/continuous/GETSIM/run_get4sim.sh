#!/bin/bash
tar -xzf R413.tar.gz
tar -xzf packages_med.tar.gz
tar -xzf SLIBS-el7.tar.gz

export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
export R_LIBS=$PWD/packages
export LD_LIBRARY_PATH=$PWD/SS:$LD_LIBRARY_PATH

Rscript getRslt_sim.R $1 $2 $3 $4 $5
