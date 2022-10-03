# testing and locating changes in variance of a process with independent data using the MIC method 






#MIC Change Point Method for parametric Exponential Distribution

findChangesExponentialMIC <- function(seq){
  #seq is the sequence we're checking for changes in
  #our type 1 error is 0.05
  MICa = NULL
  n = length(seq)
  mu = mean(seq)
  
  for(i in 2:(n-2)){
    #MIC method for exponential distributions
    lambda_1 =i/sum(x[1:i])
    r=i+1
    lambda2= (n-i) / sum(x[r:n])
    MICa[i]= -2*(i*log(lambda1) -lambda1*sum(x[1:i]) +(n-i)*log(lambda2)  -lambda2*sum(x[r:n])) +((2*i/n-1)^2)*log(n)
  }
  
  MIC_alt <-(min(MICa, na.rm = TRUE))
  MIC_null <- n*log(2*pi)+n*log(sum((x-mu)^2)/n)+n+1*log(n)
  
  
  return(ifelse(MICo - MICk +2*log(n)>qchisq(0.95,2), 1, 0))
}


findChangesNormalMIC <- function(seq){
  #write code for the MIC method whens used for exponential parametric case
}


