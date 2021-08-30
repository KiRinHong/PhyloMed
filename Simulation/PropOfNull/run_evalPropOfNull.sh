#!/bin/bash

# run your script
Rscript evalPropNull.R 1 0 0
Rscript evalPropNull.R 1 1 0
Rscript evalPropNull.R 1 0 1
Rscript evalPropNull.R 2 1 0
Rscript evalPropNull.R 2 0 1
Rscript evalPropNull.R 1 0.001 0
Rscript evalPropNull.R 2 0.001 0
Rscript evalPropNull.R 1 0 0.5
Rscript evalPropNull.R 2 0 0.5

