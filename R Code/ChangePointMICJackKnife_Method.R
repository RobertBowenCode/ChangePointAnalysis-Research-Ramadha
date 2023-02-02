# testing and locating changes in variance of a process with independent data using the MIC method with JackKnife Likelyhood  
#
#
#
#Author Robert Bowen 10/6/22
#
#MIC Change Point Method for parametric Exponential Distribution
# 
# 





#MIC Change Point Method for parametric Exponential Distribution with JackKnife likelyhood applied



findChangesExponentialMICJackKnife <- function(seq, nom_alpha){ # This should be O(n^2) pretty slow, simulations will take a lot of time unfortunately. 
  #seq is the sequence that is being determined to have a change point or not
  #nom_alpha is the nominal alpha that we would like to use in our hypothesis test
  
  #variables
  MICJa = NULL
  MICJ_probs_alt = NULL
  MICJ_probs_null = NULL
  
  for(j in 1:length(seq)) #apply the JackKnife Method to the alternative likeylhood
  {
    
    jackknife_seq = seq[-c(j)] #remove one value
    
    n = length(jackknife_seq)
    mu = mean(jackknife_seq)
  
    for(i in 2:(n-2)){ #calculate the MIC model at each changepoint i for a removed value of j
      
      
      #parameter estimations for theta_1 and theta_2
      lambda1 =i/sum(seq[1:i]) 
      r=i+1
      lambda2= (n-i) / sum(seq[r:n]) #calc the estimate for the 2nd dist
      
      
      MICJa[i]= -2*(i*log(lambda1) -lambda1*sum(jackknife_seq[1:i]) +(n-i)*log(lambda2)  -lambda2*sum(jackknife_seq[r:n]) ) + 2*log(n) + ((2*i/n-1)^2)*log(n) 
    }           #log likelyhood -2*(L(theta1,theta2,k))                             +         #complexity(ˆθ1k,ˆθ2k, k) = 2dim(ˆθ1k) + (2k/ n − 1)^2 (log(n)), where dim(ˆθ1k) is one
      
    MICJ_probs_alt[j] = min(MICJa, na.rm = TRUE) #store the min for this removed value of k
  }
  
  
  for(j in 1:length(seq)) #apply the JackKnife Method to the null log likelyhood
  {
    jackknife_seq = seq[-c(j)] #remove one value
    
    n = length(jackknife_seq)
    mu = mean(jackknife_seq)
    lambda_n = 1/mu
    
    MICJ_probs_null[j] = -2*(n*log(lambda_n) - (lambda_n)*sum(jackknife_seq[1:n])) + log(n)*1
  }
  
  
  MIC_alt <-(mean(MICJ_probs_alt)) #take the average of the minimized models for Jackknife
  MIC_null <- mean(MICJ_probs_null)  #take the jacknife model for the null hypothesis
  
  
  return(ifelse(MIC_null - MIC_alt + log(n)>qchisq(1 - nom_alpha,3), 1, 0))
  #Sn = MIC(n) − min( 1≤k<n MIC(k)) + dim(θ) log n
  #Sn follows a chi square distribution with d degrees of freedom where d is the dimension of theta
}




findChangesNormalMICJackKnife <- function(seq, nom_alpha){
 
  
  
  
  #seq is the sequence that is being determined to have a change point or not
  #nom_alpha is the nominal alpha that we would like to use in our hypothesis test
  
  #variables
  MICJa = NULL
  MICJ_probs_alt = NULL
  MICJ_probs_null = NULL
  
  for(j in 1:length(seq)) #apply the JackKnife Method to the alternative likeylhood
  {
    
    jackknife_seq = seq[-c(j)] #remove one value
    
    n = length(jackknife_seq)

    
    for(i in 2:(n-2)){ #calculate the MIC model at each changepoint i for a removed value of j
      
      r=i+1
      MICJa[i]= (n)*log(2*pi) + i*log(var(jackknife_seq[1:i]))  +  (n-i)*log(var(jackknife_seq[r:n]))  + (n) +(3+((2*i)/n-1)^2)*log(n)
    }          
    
    MICJ_probs_alt[j] = min(MICJa, na.rm = TRUE) #store the min for this removed value of k
  }
  
  
  for(j in 1:length(seq)) #apply the JackKnife Method to the null log likelyhood
  {
    jackknife_seq = seq[-c(j)] #remove one value
    
    n = length(jackknife_seq)
    
    MICJ_probs_null[j] = n*log(2*pi) + n *log(var(jackknife_seq)) + n + 2*log(n)
  
   
  }
  
  
  MIC_alt <-(mean(MICJ_probs_alt)) #take the average of the minimized models for Jackknife
  MIC_null <- mean(MICJ_probs_null)  #take the jacknife model for the null hypothesis
  
  
  return(ifelse(MIC_null - MIC_alt >qchisq(1 - nom_alpha,2), 1, 0))

}
  







#testing out simulations for type one error to see how well it performs 
n=700
repetition =200
count = NULL
for(i in 1:repetition){
  x=rexp(n,1)
  count[i] = findChangesExponentialMICJackKnife(x, 0.05)
  
}

mean(count)




#testing for JK Normal MIC method

n=700
repetition =200
count = NULL
for(i in 1:repetition){
  x=rexp(n,1)
  count[i] = findChangesNormalMICJackKnife(x, 0.05)
  
}

mean(count)

