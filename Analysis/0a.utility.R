get_top_taxa <- function(physeq_obj, n, relative = TRUE, discard_other = FALSE, other_label = "Other"){
  
  #Define a temporary physeq object
  ps_tmp <- physeq_obj
  
  #Check for 0 entries
  smpl_sms <- phyloseq::sample_sums(ps_tmp)
  if (0 %in% smpl_sms){
    stop("Error: some samples contain 0 reads. These have to be removed to avoid
         downstream problems.")
  }
  
  #Extract the otu_table as a data.frame
  otu_tbl <- phyloseq::otu_table(ps_tmp)
  if (!phyloseq::taxa_are_rows(ps_tmp)){
    otu_tbl <- t(otu_tbl)
  }
  
  #Use relative abundances if requested
  if (relative){
    otu_tbl <- apply(otu_tbl, 2, function(x){
      x / sum (x)
    })
  }
  
  #Update the phyloseq object
  phyloseq::otu_table(ps_tmp) <- phyloseq::otu_table(otu_tbl, taxa_are_rows = T)
  
  #Get the top taxa names and discard or merge other taxa
  abun_taxa <- names(sort(phyloseq::taxa_sums(ps_tmp), decreasing = TRUE)[1:n])
  if (discard_other){
    physeq_obj <- phyloseq::prune_taxa(abun_taxa, physeq_obj)
  } else {
    to_merge <- phyloseq::taxa_names(physeq_obj)
    to_merge <- to_merge[!(to_merge %in% abun_taxa)]
    physeq_obj <- merge_taxa(physeq_obj, to_merge)
    tax_tbl <- phyloseq::tax_table(physeq_obj)
    indx <- which(row.names(tax_tbl) %in% to_merge)
    tax_tbl[indx,] <- other_label
    phyloseq::tax_table(physeq_obj) <- tax_tbl
  }
  return(physeq_obj)
}
prepareTree <- function(tree, verbose = FALSE){
  if(class(tree) != "phylo") stop("Input tree is not a phylo class!")
  tree$edge = tree$edge[order(tree$edge[,2]),] # order the tips
  if(.is_binary(tree)){
    if(verbose) cat("The phylogeny tree is already binary!\n")
    if(.is_rooted(tree, .ntaxa(tree))){
      tree.new = tree
    }else{
      outgroup = .pick_new_outgroup(tree)
      if(verbose) cat("Root the tree!\n")
      tree.new = root(tree, outgroup = outgroup, resolve.root = TRUE)
    }
  }else{
    if(.is_rooted(tree, .ntaxa(tree))){
      if(verbose) cat("Resolve the multichotomies in the order they appear on the tree!\n")
      tree.new = multi2di(tree, random = FALSE)
    }else{
      outgroup = .pick_new_outgroup(tree)
      if(verbose) cat("Root the tree!\n")
      tree.new = root(tree, outgroup = outgroup, resolve.root = TRUE)
      if(verbose) cat("Resolve the multichotomies in the order they appear on the tree!\n")
      tree.new = multi2di(tree.new, random = FALSE)
    }
  }
  tree.new = Renumber(tree.new) # conform to certain principle
  return(tree.new)
}

# calculate the number of taxa
.ntaxa <- function(tree){
  length(tree$tip.label)
}
# check whether the tree is not
.is_rooted <- function(tree, K){
  if(!is.null(tree$root.edge)) return(TRUE)
  if(tabulate(tree$edge[,1])[K+1]>2) FALSE else TRUE
}
# check whether the tree is binary
.is_binary <- function(tree){
  .ntaxa(tree)-tree$Nnode+.is_rooted(tree,.ntaxa(tree)) == 2
}
.pick_new_outgroup <- function(tree.unrooted){
  treeDT = cbind(cbind(tree.unrooted$edge, tree.unrooted$edge.length)[1:.ntaxa(tree.unrooted),], tree.unrooted$tip.label)
  colnames(treeDT) = c("from", "to", "length", "id")
  # Take the longest terminal branch as outgroup
  new.outgroup = treeDT[which.max(treeDT[, "length"]),"id"]
  return(new.outgroup)
}