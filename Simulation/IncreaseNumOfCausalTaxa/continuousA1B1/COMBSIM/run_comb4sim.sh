#!/bin/bash
tar -xzf R361.tar.gz

export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R

Rscript combRslt_sim.R $1 $2 $3 $4

