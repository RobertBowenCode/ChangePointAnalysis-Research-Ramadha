# testing and locating changes in variance of a process with independent data using the MIC method 
#
#
#
#Author Robert Bowen 10/6/22
#
#MIC Change Point Method for parametric Exponential Distribution
# 
# 


findChangesExponentialMIC <- function(seq, nom_alpha){
  #seq is the sequence that is being determined to have a change point or not
  #nom_alpha is the nominal alpha that we would like to use in our hypothesis test
  
  #variables
  MICa = NULL
  n = length(seq)
  mu = mean(seq)
  lambda_n = 1/mu
  
  for(i in 2:(n-2)){ #calculate the MIC model at each changepoint i
    
    
    #parameter estimations for theta_1 and theta_2
    lambda1 =i/sum(seq[1:i]) 
    r=i+1
    lambda2= (n-i) / sum(seq[r:n]) #calc the estimate for the 2nd dist

    
    MICa[i]= -2*(i*log(lambda1) -lambda1*sum(seq[1:i]) +(n-i)*log(lambda2)  -lambda2*sum(seq[r:n]) ) + 2*log(n) + ((2*i/n-1)^2)*log(n) 
  }           #log likelyhood -2*(L(theta1,theta2,k))                             +         #complexity(ˆθ1k,ˆθ2k, k) = 2dim(ˆθ1k) + (2k/ n − 1)^2 (log(n)), where dim(ˆθ1k) is one
  
  MIC_alt <-(min(MICa, na.rm = TRUE)) #take the best model out of each changepoint
  MIC_null <- -2*(n*log(lambda_n) - (lambda_n)*sum(seq[1:n])) + log(n)*1  #take the model for the null hypothesis
  
  
  return(ifelse(MIC_null - MIC_alt + log(n)>qchisq(1 - nom_alpha,3), 1, 0))
  #Sn = MIC(n) − min( 1≤k<n MIC(k)) + dim(θ) log n
  #Sn follows a chi square distribution with d degrees of freedom where d is the dimension of theta(why is this three here again? shouldn't it be 1, but then it doesn't work strange)
}
