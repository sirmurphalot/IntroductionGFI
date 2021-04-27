###########################################################
######### Fiducial Random Effects Simulation Study ########
############## Sigma, Beta unknown ########################
###########################################################
########## Author: Alexander Murph, UNC ###################
############## November 10th, 2020 ########################
###########################################################
library(ggplot2)
library(grid)
library(latex2exp)
library(xtable)

#### SCRIPT FOR CREATING THE VISUALIZATIONS FOR THE RANDOM EFFECTS SECTIONS
#### Note: Assumes one has simulated the data already using RESimulation.R

# Identity covariance matrix:
full_data = data.frame(value = NA, data=NA,metric=NA)
all_acceptance_ratios = c()
for(name in list.files(path = "Data/")[-1]){
  temp_data = read.csv(paste("Data/",name, sep =""))[,-1]
  # temp_data = temp_data[complete.cases(temp_data),]
  full_data = rbind(full_data, temp_data)
}
full_data = full_data[-1,]

#### GATHER ALL THE RESULTS
results_siga = data.frame(data = NA, coverage95 = NA, aveLength = NA)
results_sige = data.frame(data = NA, coverage95 = NA, aveLength = NA)
for(data_type in unique(temp_data$data)){
  temp_temp_data = full_data[which( (full_data$data==data_type)&(full_data$metric%in%c("Length of 95% CI on Sigma_a2"))),]
  mean_length = median(temp_temp_data$value)
  temp_temp_data = full_data[which( (full_data$data==data_type)&(full_data$metric%in%c("95% Sigma_a2 Containment"))),]
  coverage = mean(temp_temp_data$value)
  new_row1 = data.frame(data=data_type,coverage95=coverage,aveLength = mean_length)
  
  temp_temp_data = full_data[which( (full_data$data==data_type)&(full_data$metric%in%c("Length of 95% CI on Sigma_e2"))),]
  mean_length = median(temp_temp_data$value)
  temp_temp_data = full_data[which( (full_data$data==data_type)&(full_data$metric%in%c("95% Sigma_e2 Containment"))),]
  coverage = mean(temp_temp_data$value)
  new_row2 = data.frame(data=data_type,coverage95=coverage,aveLength = mean_length)
  
  results_siga = rbind(results_siga,new_row1)
  results_sige = rbind(results_sige,new_row2)
}

results_siga = results_siga[-1,]
results_sige = results_sige[-1,]
results_siga = rbind(results_siga,results_sige)
results_siga$eta = "eta>=1"
eta_low = c("(1, 0.1, 10)",
            "(2, 0.1, 10)",
            "(3, 0.1, 10)",
            "(4, 0.1, 10)",
            "(5, 0.1, 10)",
            "(6, 0.1, 10)",
            "(7, 0.1, 10)",
            "(1, 0.5, 10)",
            "(2, 0.5, 10)",
            "(3, 0.5, 10)",
            "(4, 0.5, 10)",
            "(5, 0.5, 10)",
            "(6, 0.5, 10)", 
            "(7, 0.5, 10)",
            "(1, 0.5, 2)",
            "(2, 0.5, 2)",
            "(3, 0.5, 2)", 
            "(4, 0.5, 2)",
            "(5, 0.5, 2)",
            "(6, 0.5, 2)",
            "(7, 0.5, 2)")
results_siga$eta[which(results_siga$data%in%eta_low)]="eta<1"
results_siga$eta = as.character(results_siga$eta)

#### CREATE PLOTS
ggplot(results_siga, aes(x = eta, y = coverage95)) + geom_boxplot(fill = "skyblue")+ 
  scale_x_discrete(labels = c('eta<1' = expression(eta<1),'eta>=1'= expression(eta>=1))) +
  xlab("") + ylab("Empirical Coverage of 95% Confidence Interval")+ 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=17), legend.title = element_text(size=20), legend.text = element_text(size=17))

ggplot(results_siga, aes(x = log(aveLength+1), y = data)) + geom_boxplot(fill = "cyan") + 
  ylab(TeX("data source:(Group, $\\sigma^2_a$, $\\sigma^2_e$)")) +
  xlab(TeX("$\\log_{10}(x+1)$ of CI Lengths"))


CI_length_data = full_data[which(full_data$metric == "Length of 95% CI on Sigma_a2"),]
desc = function(x){print(paste("IQR",IQR(x),", Median", median(x), "\n"))}
tapply(CI_length_data$value, CI_length_data$data, desc)

