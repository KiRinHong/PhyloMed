# option 1
# you can pull the docker image from Docker hub or create by yourself
# > docker pull kirinhong/curate:r40
# > cat Dockerfile
##### begin of the Dockerfile #####
# FROM r-base:4.0.1
# RUN apt-get -y update \
# && apt-get install -y libreadline-dev \
# && apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libbz2-dev liblzma-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libgit2-dev
# RUN install2.r --error --skipinstalled --ncpus -1 \
# devtools \
# BiocManager
# RUN R -e 'BiocManager::install(version = '3.11', ask = FALSE, update = FALSE)'
# RUN R -e 'BiocManager::install("curatedMetagenomicData", version = "3.11")'
# RUN R -e 'BiocManager::install("phyloseq", version = "3.11")'
# RUN rm -rf /tmp/downloaded_packages
# RUN apt-get clean && rm -rf /var/cache/apt/lists
##### end of the Dockerfile #####
# > docker build -t kirinhong/curate:r40 .  
# > docker run -it --rm=true kirinhong/curate:r40 /bin/bash
# > R 
# run the following code in the Docker Image
library(curatedMetagenomicData)
library(phyloseq)
zeeviD <- curatedMetagenomicData("ZeeviD_2015.metaphlan_bugs_list.stool", counts = TRUE, dryrun = FALSE)
zeeviD.eset <- zeeviD[[1]]
zeeviD.pseq <- ExpressionSet2phyloseq( zeeviD.eset, phylogenetictree = TRUE)
save(zeeviD.pseq, file="zeeviD_pseq.Rdata")
# transfer the zeeviD_pseq.Rdata from container to local 
# > docker cp <CONTAINER ID>:/zeeviD_pseq.Rdata ./ 
# run the following code in this R script
rm(list = ls())
setwd("~/Documents/Project/PhyloMed/Simulation/")
library(phyloseq)
load("zeeviD_pseq.Rdata")
zeeviD.pseq
# only keep Bacteria
zeeviD.pseq.bact = subset_taxa( zeeviD.pseq,  Kingdom=="Bacteria" )
zeeviD.pseq.bact
save(zeeviD.pseq.bact, file="RareTaxa/continuous/GETSIM/shared/zeeviD_pseq_bacteria.Rdata")

# only keep top 100 species
count = t(otu_table(zeeviD.pseq.bact))
P = count/rowSums(count)
P.avg = colMeans(P)
top100.idx = order(P.avg, decreasing = TRUE)[1:100]
keepotu = rep(FALSE, ncol(count))
keepotu[top100.idx] = TRUE
names(keepotu) = colnames(P)
zeeviD.pseq.bact.com = subset_taxa(zeeviD.pseq.bact, keepotu)

count = t(otu_table(zeeviD.pseq.bact.com))
tax = tax_table(zeeviD.pseq.bact.com)
phy.tree = phy_tree(zeeviD.pseq.bact.com)
zeeviD = list(count=count, tax=tax, phy.tree=phy.tree)
save(zeeviD, file="zeeviD.Rdata")

count.full = t(otu_table(zeeviD.pseq.bact))
sum(colMeans(count.full/rowSums(count.full))[colnames(count)]) # 97.82%
