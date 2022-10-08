

##################Method of Hypothesis test ####################################



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
    
    MICJ_probs_null[j] = -2*(n*log(lambda_n) - (lambda_n)*sum(seq[1:n])) + log(n)*1
  }
  
  
  MIC_alt <-(mean(MICJ_probs_alt)) #take the average of the minimized models for Jackknife
  MIC_null <- mean(MICJ_probs_null)  #take the jacknife model for the null hypothesis
  
  
  return(ifelse(MIC_null - MIC_alt + log(n)>qchisq(1 - nom_alpha,3), 1, 0))
  #Sn = MIC(n) − min( 1≤k<n MIC(k)) + dim(θ) log n
  #Sn follows a chi square distribution with d degrees of freedom where d is the dimension of theta
}













##############################################################################


#ChangePoint Simulation with Exp(1) on MIC JackKnife Method
#simulating to calc empirical type 1 error at various sample sizes and nominal values of alpha

#repetitions =500
# testing at nominal alphas = 0.01, 0.05, 0.1

###Variables###
exp_mean = 1

#Store the results
simulations_n_50  = c()
simulations_n_100 = c()
simulations_n_150 = c()
simulations_n_200 = c()



sample_sizes = c(50,100,150,200)
nominal_type_1_error = c(0.01, 0.05, 0.1)
repetitions = 500


for (sizes in sample_sizes) #for each sample size
{
  
  array_index = 1
  
  
  for(type_1 in nominal_type_1_error) #calculate the empirical type 1 error at nominal values of 0.01, 0.05, and 0.1
  {

    
    empirical_type_1 = 0
    count = 0
    
    for(i in 1:repetitions)
    {
      seq <-rexp(sizes, exp_mean) #simulate an exponential dataset with no changes
      count = count + findChangesExponentialMICJackKnife(seq, type_1) #run test on dataset to check for change
    }
    
    empirical_type_1 = count/repetitions
    
    
    #store the results
    if(sizes == 50)
    {
      simulations_n_50[array_index] = empirical_type_1
    }
    else if(sizes == 100)
    {
      simulations_n_100[array_index] = empirical_type_1
    }
    else if (sizes == 150)
    {
      simulations_n_150[array_index] = empirical_type_1
    }
    else if (sizes == 200)
    {
      simulations_n_200[array_index] = empirical_type_1
      
    }
    
    array_index = array_index+1
    
  }
  
  
}













