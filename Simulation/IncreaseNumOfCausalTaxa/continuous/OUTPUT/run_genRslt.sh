#!/bin/bash
tar -xzf R413.tar.gz
tar -xzf packages_med.tar.gz

export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
export R_LIBS=$PWD/packages

Rscript generateRslt.R
