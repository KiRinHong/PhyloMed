

############### Subfunction ############################
## Compute the estimate of pi0 for a fixed value of B ##
########################################################
FixedB <- function(p, B)
{
  ## Input:
  #p: a vector of p-values
  #B: an integer, the interval [0,1] is divided into B equal-length intervals
  
  ## Output:
  #pi0: an estimate of the proportion of true null hypotheses
  
  m <- length(p) 
  t <- seq(0,1,length=B+1)   # equally spaced points in the interval [0,1]
  
  NB <- rep(0,B)		    # number of p-values greater than t_i	
  NBaverage <- rep(0,B)      # average number of p-values in each of the (B-i+1) small intervals on [t_i,1]
  NS <- rep(0,B)             # number of p-values in the interval [t_i, t_(i+1)]
  pi <- rep(0,B)		    # estimates of pi0 	
  for(i in 1:B)
  {
    NB[i] <- length(p[p>=t[i]])
    NBaverage[i] <- NB[i]/(B-(i-1))
    NS[i] <- length(p[p>=t[i]]) - length(p[p>=t[i+1]])
    pi[i] <- NB[i]/(1-t[i])/m 
  }
  
  i <- min(which(NS <= NBaverage))  # Find change point i
  pi0 <- min(1, mean(pi[(i-1):B]))          # average estiamte of pi0
  return(pi0)
}
AverageEstimate <- function(p=NULL, Bvector=c(5, 10, 20, 50, 100))
{
  ## Input:
  #p: a vector of p-values
  #Bvector: a vector of integer values where the interval [0,1] is divided into B equal-length intervals
  #         When Bvector is an integer, number of intervals is consider as fixed. For example Bvector = 10;
  #         When Bvector is a vector, bootstrap method is used to choose an optimal value of B
  #Numboot: number of bootstrap samples
  
  ## Output:
  #pi0: an estimate of the proportion of true null hypotheses
  
  # check if the p-values are valid
  if (min(p)<0 || max(p)>1) print("Error: p-values are not in the interval of [0,1]")
  m <- length(p) 		# Total number p-values
  
  Bvector <- as.integer(Bvector)  # Make sure Bvector is a vector of integers
  
  #Bvector has to be bigger than 1
  if(min(Bvector) <=1) print ("Error: B has to be bigger than 1")
  
  ######## Estimate pi0 ########
  if(length(Bvector) == 1)        # fixed number of numbers, i.e., B is fixed
  {
    pi0 <- FixedB(p, Bvector)
  }
  else
  {
    Numboot <- 100
    OrigPi0Est <- rep(0, length(Bvector))
    for (Bloop in 1:length(Bvector))
    {
      OrigPi0Est[Bloop] <- FixedB(p, Bvector[Bloop])
    }
    
    BootResult <- matrix(0, nrow=Numboot, ncol=length(Bvector)) # Contains the bootstrap results
    
    for(k in 1:Numboot)
    {
      p.boot <- sample(p, m, replace=TRUE)    # bootstrap sample 
      for (Bloop in 1:length(Bvector))
      {
        BootResult[k,Bloop] <- FixedB(p.boot, Bvector[Bloop])
      }	
    }
    
    MeanPi0Est <- mean(OrigPi0Est)             # avearge of pi0 estimates over the range of Bvector
    MSEestimate <- rep(0, length(Bvector))     # compute mean-squared error
    for (i in 1:length(Bvector))
    {
      MSEestimate[i] <- (OrigPi0Est[i]- MeanPi0Est)^2
      for (k in 1:Numboot)
      {
        MSEestimate[i] <- MSEestimate[i]+1/Numboot*(BootResult[k,i] - OrigPi0Est[i])^2	
      }
    }
    pi0 <- OrigPi0Est[MSEestimate==min(MSEestimate)]
  }  # end of else           
  return(mean(pi0))
}

FixedB_mod <- function(p, B)
{
  ## Input:
  #p: a vector of p-values
  #B: an integer, the interval [0,1] is divided into B equal-length intervals
  
  ## Output:
  #pi0: an estimate of the proportion of true null hypotheses
  
  m <- length(p) 
  t <- seq(0,1,length=B+1)   # equally spaced points in the interval [0,1]
  
  NB <- rep(0,B)		    # number of p-values greater than t_i	
  NBaverage <- rep(0,B)      # average number of p-values in each of the (B-i+1) small intervals on [t_i,1]
  NS <- rep(0,B)             # number of p-values in the interval [t_i, t_(i+1)]
  pi <- rep(0,B)		    # estimates of pi0 	
  for(i in 1:B)
  {
    NB[i] <- length(p[p>=t[i]])
    NBaverage[i] <- NB[i]/(B-(i-1))
    NS[i] <- length(p[p>=t[i]]) - length(p[p>=t[i+1]])
    pi[i] <- NB[i]/(1-t[i])^2/m 
  }
  
  i <- min(which(NS <= NBaverage))  # Find change point i
  pi0 <- min(1, mean(pi[(i-1):B]))          # average estiamte of pi0
  return(pi0)
}
AverageEstimate_mod <- function(p=NULL, Bvector=c(5, 10, 20, 50, 100))
{
  ## Input:
  #p: a vector of p-values
  #Bvector: a vector of integer values where the interval [0,1] is divided into B equal-length intervals
  #         When Bvector is an integer, number of intervals is consider as fixed. For example Bvector = 10;
  #         When Bvector is a vector, bootstrap method is used to choose an optimal value of B
  #Numboot: number of bootstrap samples
  
  ## Output:
  #pi0: an estimate of the proportion of true null hypotheses
  
  # check if the p-values are valid
  if (min(p)<0 || max(p)>1) print("Error: p-values are not in the interval of [0,1]")
  m <- length(p) 		# Total number p-values
  
  Bvector <- as.integer(Bvector)  # Make sure Bvector is a vector of integers
  
  #Bvector has to be bigger than 1
  if(min(Bvector) <=1) print ("Error: B has to be bigger than 1")
  
  ######## Estimate pi0 ########
  if(length(Bvector) == 1)        # fixed number of numbers, i.e., B is fixed
  {
    pi0 <- FixedB_mod(p, Bvector)
  }
  else
  {
    Numboot <- 100
    OrigPi0Est <- rep(0, length(Bvector))
    for (Bloop in 1:length(Bvector))
    {
      OrigPi0Est[Bloop] <- FixedB_mod(p, Bvector[Bloop])
    }
    
    BootResult <- matrix(0, nrow=Numboot, ncol=length(Bvector)) # Contains the bootstrap results
    
    for(k in 1:Numboot)
    {
      p.boot <- sample(p, m, replace=TRUE)    # bootstrap sample 
      for (Bloop in 1:length(Bvector))
      {
        BootResult[k,Bloop] <- FixedB_mod(p.boot, Bvector[Bloop])
      }	
    }
    
    MeanPi0Est <- mean(OrigPi0Est)             # avearge of pi0 estimates over the range of Bvector
    MSEestimate <- rep(0, length(Bvector))     # compute mean-squared error
    for (i in 1:length(Bvector))
    {
      MSEestimate[i] <- (OrigPi0Est[i]- MeanPi0Est)^2
      for (k in 1:Numboot)
      {
        MSEestimate[i] <- MSEestimate[i]+1/Numboot*(BootResult[k,i] - OrigPi0Est[i])^2	
      }
    }
    pi0 <- OrigPi0Est[MSEestimate==min(MSEestimate)]
  }  # end of else           
  return(mean(pi0))
}
