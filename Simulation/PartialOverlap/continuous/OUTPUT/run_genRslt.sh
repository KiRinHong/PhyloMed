#!/bin/bash
tar -xzf R361.tar.gz

export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R

Rscript generateRslt.R
