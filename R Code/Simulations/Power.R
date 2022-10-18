# testing and locating changes in variance of a process with independent data using the MIC method with JackKnife Likelyhood  
#
#
#
#Author Robert Bowen 10/6/22
# Power Simulation at different samples size and locations of changes in the data. 
# MIC Change Point Method for parametric Exponential Distribution
# 
# 

#MIC Change Point Method for parametric Exponential Distribution with JackKnife likelyhood applied




##the Exponential MIC JackKnife method

findChangesExponentialMICJackKnife <- function(seq, nom_alpha){
  # This should be O(n^2) pretty slow, simulations will take a lot of time unfortunately. 
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














#################################################################
#
#
#Power Simulation for JackKnife Method. 
#
#
#
#


alpha=0.05                    #nominal type 1 error
location_set =seq(0,1,by=0.1) # the location of the split in dataset
size_set=c(25,50,75,100)      # the total sample sizes. 
power_index = 1

power_list  = c()

for(sample_size in size_set ){ ##for each sample size
  
  
  power=NULL
  
  for(change_location in location_set){ #for each different location of change point
    
    
    #get samples lengths, at chosen change point split
    firstSamplen=floor(location_set[change_location]*sample_size)
    secondSamplen=sample_size-firstSamplen
    
    
    repetition=50000 #number of repetitions
    
    
    for(i in 1:repetition){ #create the simulations
      
      #do the hypothesis test
      first_sample=rexp(firstSamplen,1)
      second_sample=rexp(secondSamplen,3)
      seq =c(first_sample,second_sample)
      count[i] = findChangesExponentialMICJackKnife(seq, alpha)
      
    }
    print(q[j])
    power=c(power,mean(count))
  }
  
  
  power_list[power_index] = power
  power_index = power_index + 1; 
  
  
  if(k!=1){par(new=TRUE)}
  plot(q,power,type='l',ylim=c(0,1),xlab="Changepoint Location",main=c("Sample 1: Exp(1)",
                                                                       "Sample 2: Exp(3)",
                                                                       "alpha=0.05"),col=k)
}























#Power simulations for MIC method. 

MICa = NULL
count = c()
# Function to estimate single change location with significance level alpha

##this is the method
findChangeLR <- function(seq){
  n = length(seq)
  for(i in 2:(n-2)){
    lambda1=i/sum(x[1:i])
    r=i+1
    lambda2=(n-i)/sum(x[r:n])
    MICa[i]= -2*(i*log(lambda1)   -lambda1*sum(x[1:i])   +(n-i)*log(lambda2)   -lambda2*sum(x[r:n]))   +((2*i/n-1)^2)*log(n)
  }
  return(min(MICa, na.rm = TRUE))
}

alpha=0.05
q=seq(0,1,by=0.1)
n=c(25,50,75,100)

#power is the probability of rejecting the null hypothesis given that its false. 

for(k in 1:4){ ##for each sample size
  power=NULL
  
  for(j in 1:length(q)){ #for each different location of change point
    
    
    #get samples lengths, this is 
    firstSamplen=floor(q[j]*n[k])
    secondSamplen=n[k]-firstSamplen
    
    
    repetition=50000 #number of repitions
    
    
    for(i in 1:repetition){ #create the simulations
      #do the hypothesis test
      z=rexp(firstSamplen,1)
      y=rexp(secondSamplen,3)
      x=c(z,y)
      mu = mean(x)
      lambda=(n[k])/sum(x)
      MICo = -2*(n[k]*log(lambda)-lambda*sum(x))
      MICk = findChangeLR(x)
      count[i] = ifelse(MICo - MICk>qchisq(1-alpha,3), 1, 0)
      
    }
    print(q[j])
    power=c(power,mean(count))
  }
  if(k!=1){par(new=TRUE)}
  plot(q,power,type='l',ylim=c(0,1),xlab="Changepoint Location",main=c("Sample 1: Exp(1)",
                                                                       "Sample 2: Exp(3)",
                                                                       "alpha=0.05"),col=k)
}
abline(h=0)
abline(h=alpha,lty=2)