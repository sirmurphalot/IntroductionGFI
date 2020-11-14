###############################################
######### Fiducial MVN Simulation Study ######
############## Sigma Unknown #################
##############################################
########## Author: Alexander Murph, UNC ######
############## November 10th, 2020 ###########
##############################################

library(mvtnorm)
library(rstan)
library(geigen)

### Uncomment these if you're working with SLURM cloud computing.
# i = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#print(paste("Running sigma value:", i))
# set.seed(i)
### Comment this if you're working with SLURM cloud computing.
i = 1
set.seed(i)
prev_i = i

## Helper Functions (by Murph)
make_cov_matrix = function(n_rows, case){
  # Construct the "truth".  Cases: 1,2,3.
  T_m = diag(n_rows)
  if(case == 1){
    cov_m = T_m
  }else if( case == 2){
    for(i in 1:n_rows){
      for(j in 1:i){
        if(i == j){
          T_m[i,j]=log(i/10+2)
        }else if(j == (i-1)){
          T_m[i,j]=2*(i/n_rows)**2-0.5
        }
      }
    }
    cov_m=T_m%*%t(T_m)
  }else{
    for(i in 1:n_rows){
      for(j in 1:i){
        if(i == j){
          T_m[i,j]=log(i/10+2)
        }else{
          T_m[i,j]=n_rows^(-2)*min(c(i+j,i^(1.5)))*exp(-j/4)
        }
      }
    }
    cov_m=T_m%*%t(T_m)
  }
  return(cov_m)
}
loss_function1 = function(truth, estimate){
  temp_m = solve(truth)%*%estimate
  return( sum(diag(temp_m))-log(det(temp_m))-nrow(truth) )
}
loss_function2 = function(truth, estimate){
  temp_m = solve(truth)%*%estimate
  return( sum(diag( (temp_m-diag(nrow(truth))) ))^2 )
}
FM_distance = function(truth,estimate){
  eigen_values = geigen(truth, estimate, TRUE, only.values=TRUE)$values
  return( sqrt(sum( (log(eigen_values))**2 )) )
}
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


#####################################################
######## Code for drawing from posterior
fiducial_stan = stan_model(file='MVNormal.stan') 
n=100L
truemu=c(1,2,3,1)
truesigma=matrix(c(4,1,0,0, 1,1,0,1, 0,0,9,1, 0,1,1,4),nrow=4)
d=length(truemu)

y=t(rmvnorm(n,truemu,truesigma))
Ybar=rowMeans(y)
Sn=var(t(y))
SSE=Sn*(n-1)
Ypca=eigen(SSE)

#Clunky way to find a good starting points that satisfies sampling constraints 
n_chains=20L
my_count=0
originalsign=sign(det(Ypca$vectors))
for(i in 0:(2^d-1)){
  temp_signs=2*as.numeric(intToBits(i))[1:d]-1
  if(1==originalsign*prod(temp_signs)){
    tempUhat=Ypca$vectors%*%diag(temp_signs)
    tempYhatAm=(diag(d)-tempUhat)%*%solve(diag(d)+tempUhat)
    tempYhatav=tempYhatAm[lower.tri(tempYhatAm)]
    if(1==prod(tempYhatav>-1 & tempYhatav<1)){
      Yhatav=tempYhatav
      Uhat=tempUhat
      my_i=i
      my_count=my_count+1
    }
  }
}
init_stan=list(list(aa=Yhatav))

for(id in 2:n_chains){
  reshuffle=sample(d)
  
  my_count=0
  originalsign=sign(det(Ypca$vectors[,reshuffle]))
  for(i in 0:(2^d-1)){
    temp_signs=2*as.numeric(intToBits(i))[1:d]-1
    if(1==originalsign*prod(temp_signs)){
      tempUhat=Ypca$vectors[,reshuffle]%*%diag(temp_signs)
      tempYhatAm=(diag(d)-tempUhat)%*%solve(diag(d)+tempUhat)
      tempYhatav=tempYhatAm[lower.tri(tempYhatAm)]
      if(1==prod(tempYhatav>-1 & tempYhatav<1)){
        Yhatav=tempYhatav
        Uhat=tempUhat
        my_i=i
        my_count=my_count+1
      }
    }
  }
  init_stan[[id]]=list(aa=Yhatav)
}

data_stan=list(d=d,n=n,ybar=Ybar,SSE=SSE)
FIDfit=sampling(fiducial_stan, data = data_stan, init = init_stan, chains = n_chains, iter = 5000L, cores = 4L)
posterior_samples = extract(FIDfit)



#####################################################
########## Analysis Code
## Aligning our different parameter naming:
d_value=d
true_sigma = truesigma
# Start by getting the mean covariance matrix:
mean_cov_mat = matrix(0,ncol=d_value,nrow=d_value)
num_samples = nrow(posterior_samples$Sigma[,,1])
for(j in 1:d_value){
  mean_cov_mat[j,] = colSums(posterior_samples$Sigma[,,j])/num_samples
}


# Gather posterior distance metrics
posterior_data = data.frame(value = NA, description = NA)
#################
# Get all metrics for the true value of sigma
diff_from_truth_loss1 = loss_function1(true_sigma,mean_cov_mat)
diff_from_truth_loss2 = loss_function2(true_sigma,mean_cov_mat)
det_of_truth = determinant(truesigma)$modulus
diff_from_truth_FMdist = FM_distance(true_sigma,mean_cov_mat)
diff_from_truth_Frobenius_dist = norm(true_sigma-mean_cov_mat, type = "F")
diff_from_truth_SingularMax_dist = norm(true_sigma-mean_cov_mat,type="2")
SingularMax_of_truth = norm(true_sigma,type="2")
Frobenius_of_truth = norm(true_sigma,type="F")
# Do the same for mu
meanmu = colSums(posterior_samples$mu)/nrow(posterior_samples$mu)
diff_from_truth_mu = euc.dist(truemu, meanmu)
#
### Add in all true values
posterior_data = rbind(posterior_data,
                       data.frame(value=diff_from_truth_loss1,description="Loss_1_True_Dist"),
                       data.frame(value=diff_from_truth_loss2,description="Loss_2_True_Dist"),
                       data.frame(value=det_of_truth,description="LogDet_Truth"),
                       data.frame(value=diff_from_truth_FMdist,description="FM_True_Dist"),
                       data.frame(value=SingularMax_of_truth,description="SingularMax_Truth"),
                       data.frame(value=Frobenius_of_truth,description="Frobenius_Truth"),
                       data.frame(value=diff_from_truth_mu,description="Mu_True_Dist"))
### Add in all norms of the mean
SingularMax_mean = max(sqrt(eigen(mean_cov_mat)$values))
Frobenius_mean = norm(mean_cov_mat,type="F")
LogDet_mean = determinant(mean_cov_mat)$modulus
posterior_data = rbind(posterior_data,
                       data.frame(value=SingularMax_mean,description="SingularMax_Mean"),
                       data.frame(value=Frobenius_mean,description="Frobenius_Mean"),
                       data.frame(value=LogDet_mean,description="LogDet_Mean"))
write.csv(posterior_data,file = paste("Logs/posterior_data", prev_i, ".csv", sep="_"))

#### CALCULATE VALUES FROM THE POSTERIOR DISTRIBUTION
for(i in 1:num_samples){
  print(paste("Gathering data for sample",i))
  # First, construct the matrix from the posterior output:
  temp_matrix=matrix(0,ncol=d_value,nrow=d_value)
  for(j in 1:d_value){
    temp_matrix[j,] = posterior_samples$Sigma[,,j][i,]
  }
  temp_mu = posterior_samples$mu[i,]
  #### Loss functions (NOT norms)
  #################
  # Get Loss1 Data: 95% CI around mean & qqplot
  diff_from_mean = loss_function1(mean_cov_mat,temp_matrix)
  temp_row1 = data.frame(value=diff_from_mean, description="Loss_1_Difference")
  temp_row2 = data.frame(value=(diff_from_mean>diff_from_truth_loss1), description="Loss_1_CI_Coverage")
  posterior_data = rbind(posterior_data,temp_row1,temp_row2)
  #################
  # Get Loss2 Data: 95% CI around mean & qqplot
  diff_from_mean = loss_function2(mean_cov_mat,temp_matrix)
  temp_row1 = data.frame(value=diff_from_mean, description="Loss_2_Difference")
  temp_row2 = data.frame(value=(diff_from_mean>diff_from_truth_loss2), description="Loss_2_CI_Coverage")
  posterior_data = rbind(posterior_data,temp_row1,temp_row2)
  #################
  # Get FM Distance Data: 95% CI around mean & qqplot
  diff_from_mean = FM_distance(mean_cov_mat,temp_matrix)
  temp_row1 = data.frame(value=diff_from_mean, description="FM_Difference")
  temp_row2 = data.frame(value=(diff_from_mean>diff_from_truth_FMdist), description="FM_CI_Coverage")
  posterior_data = rbind(posterior_data,temp_row1,temp_row2)
  #################
  # Get Frobenius Distance Data: 95% CI around mean & qqplot
  diff_from_mean = norm(mean_cov_mat - temp_matrix, type = "F")
  temp_row1 = data.frame(value=diff_from_mean, description="Frobenius_Difference")
  temp_row2 = data.frame(value=(diff_from_mean>diff_from_truth_Frobenius_dist), description="Frobenius_Diff_CI_Coverage")
  posterior_data = rbind(posterior_data,temp_row1,temp_row2)
  #################
  # Get Spectral Norm Distance Data: 95% CI around mean & qqplot
  diff_from_mean = norm(mean_cov_mat - temp_matrix, type = "2")
  temp_row1 = data.frame(value=diff_from_mean, description="SpectralNorm_Difference")
  temp_row2 = data.frame(value=(diff_from_mean>diff_from_truth_SingularMax_dist), description="Spectral_Diff_CI_Coverage")
  posterior_data = rbind(posterior_data,temp_row1,temp_row2)
  #################
  # Get Distance from Mu Data: 95% CI around mean & qqplot
  diff_from_mean = euc.dist(meanmu,temp_mu)
  temp_row1 = data.frame(value=diff_from_mean, description="Mu_Difference")
  temp_row2 = data.frame(value=(diff_from_mean>diff_from_truth_mu), description="Mu_CI_Coverage")
  posterior_data = rbind(posterior_data,temp_row1,temp_row2)
  
  #### Norms
  #################
  # Get LogDet Data: 95% CI around mean & qqplot
  logDet_of_temp = determinant(temp_matrix)$modulus
  temp_row1 = data.frame(value=logDet_of_temp, description="LogDet_value")
  temp_row2 = data.frame(value=(logDet_of_temp>det_of_truth), description="LogDet_UpperCI_Coverage")
  posterior_data = rbind(posterior_data,temp_row1,temp_row2)
  #################
  # Get SingularMax Data: 95% CI around mean & qqplot
  SingularMax_of_temp = norm(temp_matrix, type = "2")
  temp_row1 = data.frame(value=SingularMax_of_temp, description="SingularMax_value")
  temp_row2 = data.frame(value=(SingularMax_of_temp>SingularMax_of_truth), description="SingularMax_UpperCI_Coverage")
  posterior_data = rbind(posterior_data,temp_row1,temp_row2)
  #################
  # Get Frobenius Data: 95% CI around mean & qqplot
  Frobenius_of_temp = norm(temp_matrix, type = "F")
  temp_row1 = data.frame(value=Frobenius_of_temp, description="Frobenius_value")
  temp_row2 = data.frame(value=(Frobenius_of_temp>Frobenius_of_truth), description="Frobenius_UpperCI_Coverage")
  posterior_data = rbind(posterior_data,temp_row1,temp_row2)
}
posterior_data=posterior_data[-1,]

# Use Posterior Data to calculate coverage probabilities:
CI_loss1 = mean(posterior_data$value[posterior_data$description=="Loss_1_CI_Coverage"])
Containment_loss1=quantile(posterior_data$value[posterior_data$description=="Loss_1_Difference"],
                           probs=0.95)>diff_from_truth_loss1 
CI_loss2 = mean(posterior_data$value[posterior_data$description=="Loss_2_CI_Coverage"])
Containment_loss2=quantile(posterior_data$value[posterior_data$description=="Loss_2_Difference"],
                           probs=0.95)>diff_from_truth_loss2 
CI_FM = mean(posterior_data$value[posterior_data$description=="FM_CI_Coverage"])
Containment_FM=quantile(posterior_data$value[posterior_data$description=="FM_Difference"],
                           probs=0.95)>diff_from_truth_FMdist 
CI_Frobenius_Diff = mean(posterior_data$value[posterior_data$description=="Frobenius_Diff_CI_Coverage"])
Containment_Frobenius_Diff=quantile(posterior_data$value[posterior_data$description=="Frobenius_Difference"],
                        probs=0.95)>diff_from_truth_Frobenius_dist
CI_Spectral_Diff = mean(posterior_data$value[posterior_data$description=="Spectral_Diff_CI_Coverage"])
Containment_Spectral_Diff=quantile(posterior_data$value[posterior_data$description=="SpectralNorm_Difference"],
                        probs=0.95)>diff_from_truth_SingularMax_dist


CI_mu = mean(posterior_data$value[posterior_data$description=="Mu_CI_Coverage"])
Containment_mu=quantile(posterior_data$value[posterior_data$description=="Mu_Difference"],
                        probs=0.95)>diff_from_truth_mu 
CI_LogDet = mean(posterior_data$value[posterior_data$description=="LogDet_UpperCI_Coverage"])
Containment_LogDet = (quantile(posterior_data$value[posterior_data$description=="LogDet_value"],
                        probs=0.025)<det_of_truth)&(quantile(posterior_data$value[posterior_data$description=="LogDet_value"],
                                                             probs=0.975)>det_of_truth)
CI_SingularMax = mean(posterior_data$value[posterior_data$description=="SingularMax_UpperCI_Coverage"])
Containment_SingularMax= (quantile(posterior_data$value[posterior_data$description=="SingularMax_value"],
                              probs=0.025)<SingularMax_of_truth)&(quantile(posterior_data$value[posterior_data$description=="SingularMax_value"],
                                                                   probs=0.975)>SingularMax_of_truth)
CI_Frobenius = mean(posterior_data$value[posterior_data$description=="Frobenius_UpperCI_Coverage"])
Containment_Frobenius= (quantile(posterior_data$value[posterior_data$description=="Frobenius_value"],
                              probs=0.025)<Frobenius_of_truth)&(quantile(posterior_data$value[posterior_data$description=="Frobenius_value"],
                                                                   probs=0.975)>Frobenius_of_truth)

#### COMPILE DATA AND WRITE TO CSV
final_data = data.frame(value = c(CI_loss1, Containment_loss1,
                                  CI_loss2, Containment_loss2,
                                  CI_FM, Containment_FM,
                                  CI_LogDet, Containment_LogDet,
                                  CI_SingularMax, Containment_SingularMax,
                                  CI_Frobenius, Containment_Frobenius,
                                  CI_mu, Containment_mu,
                                  CI_Frobenius_Diff, Containment_Frobenius_Diff,
                                  CI_Spectral_Diff, Containment_Spectral_Diff),
                        description = c("QQ_loss1", "Loss1_95CI_ContainsTruth",
                                        "QQ_loss2", "Loss2_95CI_ContainsTruth",
                                        "QQ_FM", "FM_95CI_ContainsTruth",
                                        "QQ_LogDet", "LogDet_95CI_ContainsTruth",
                                        "QQ_SingularMax", "SingularMax_95CI_ContainsTruth",
                                        "QQ_Frobenius", "Frobenius_95CI_ContainsTruth",
                                        "QQ_Mu", "Mu_95CI_ContainsTruth",
                                        "QQ_Frobenius_Diff", "Frobenius_Diff_95CI_ContainsTruth",
                                        "QQ_Spectral_Diff", "Spectral_Diff_95CI_ContainsTruth"),
                        iteration = rep(prev_i, times = 18))

write.csv(final_data, file = paste("Data/MVNSTAN", prev_i, ".csv", sep="_"))
