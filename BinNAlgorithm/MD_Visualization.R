###################################################
########## Fiducial Binomial MD Visuals ###########
###################################################
########### Author: Alexander Murph, UNC ##########
###################################################
library(ggplot2)
list.files(path = ".")
full_data = data.frame(P = NA, MD = NA, parameter = NA, data_size = NA,
                       upper_p = NA, lower_p = NA)

for(name in list.files(path = "MD_Data/")[-1]){
  temp_data = read.csv(paste("MD_Data/",name, sep =""))[,-1]
  full_data = rbind(full_data, temp_data)
}
full_data = full_data[-which(is.na(full_data$P)),]
full_data$P = as.factor(full_data$P)
full_data$parameter[which(full_data$parameter=='bayesian-n')] = "bayesian"
full_data$parameter[which(!(full_data$parameter=='bayesian'))] = "fiducial"

full_data$approach = full_data$parameter
full_data$data_size = as.character(full_data$data_size)
full_data$data_size[which(full_data$data_size=="10")] = "size = 10"
full_data$data_size[which(full_data$data_size=="50")] = "size = 50"
full_data$data_size[which(full_data$data_size=="100")] = "size = 100"
full_data$data_size = factor(full_data$data_size, levels=c("size = 10","size = 50","size = 100"))

# All sizes MD:

ggplot(full_data, aes(x = P, y = MD, color = approach)) + ylab("MAD") + 
  geom_boxplot() + facet_grid( data_size ~ .) + 
  theme(text = element_text(size=30), axis.text.x = element_text(angle = 90), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20))


full_data = data.frame(P = NA, MD = NA, parameter = NA, data_size = NA,
                       upper_p = NA, lower_p = NA)

for(name in list.files(path = "MD_Data/")[-1]){
  temp_data = read.csv(paste("MD_Data/",name, sep =""))[,-1]
  full_data = rbind(full_data, temp_data)
}
full_data = full_data[-which(is.na(full_data$P)),]
full_data$P = as.factor(full_data$P)
full_data$parameter[which(full_data$parameter=='bayesian-n')] = "bayesian"
full_data$parameter[which(!(full_data$parameter=='bayesian'))] = "fiducial"

size = 100
my_data = full_data[which(full_data$data_size ==size),]
ggplot(my_data, aes(x = P, y = MD, color = parameter)) + geom_boxplot() + 
  ggtitle(paste("Mean Difference; Data Size = ", size, sep="")) + 
  theme(text = element_text(size=30), axis.text.x = element_text(angle = 90), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20))

ggplot(my_data, aes(x = upper_p, color = P, linetype = parameter)) + stat_ecdf() + 
  xlab("Nominal Coverage") + ylab("Empirical Coverage")+ 
  guides(linetype=guide_legend(title="Approach")) +
  scale_color_manual(values=unlist(lapply(c(0.4, 0.5, 0.55, 0.6, 0.75, 0.8, 0.85, 0.9, 0, 0.05, 0.12, .15), colors))) + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20))


my_data = read.csv("binomial_simulation_data.csv")
P_values = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 
             0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
TrueN = 10
num_of_draws = 100
size_of_data = c(10, 50, 100)
full_data = my_data[-1,]
full_data$P = as.factor(full_data$P)
error_bar_data = data.frame(P = NA, n_min = NA, n_max = NA, parameter = NA, data_size = NA)
for(P in P_values){
  for(size in size_of_data){
    ymin = quantile(full_data[which(full_data$P == P & full_data$data_size == size & full_data$parameter == "bayesian-n"),]$n, 
                    prob = 0.05)
    ymax = quantile(full_data[which(full_data$P == P & full_data$data_size == size & full_data$parameter == "bayesian-n"),]$n, 
                    prob = 0.95)
    new_row1 = data.frame(P = as.factor(P), n_min = ymin, n_max = ymax, parameter = "bayesian-n", data_size = size)
    ymin = quantile(full_data[which(full_data$P == P & full_data$data_size == size & full_data$parameter == "fiducial-n"),]$n, 
                    prob = 0.05)
    ymax = quantile(full_data[which(full_data$P == P & full_data$data_size == size & full_data$parameter == "fiducial-n"),]$n, 
                    prob = 0.95)
    new_row2 = data.frame(P = as.factor(P), n_min = ymin, n_max = ymax, parameter = "fiducial-n", data_size = size)
    error_bar_data = rbind(error_bar_data, new_row1, new_row2)
  }
}
error_bar_data = error_bar_data[-1,]
error_bar_data$parameter[which(error_bar_data$parameter=='bayesian-n')] = "bayesian"
error_bar_data$parameter[which(!(error_bar_data$parameter=='bayesian'))] = "fiducial"
error_bar_data$approach = error_bar_data$parameter
error_bar_data$data_size = as.character(error_bar_data$data_size)
error_bar_data$data_size[which(error_bar_data$data_size=="10")] = "size = 10"
error_bar_data$data_size[which(error_bar_data$data_size=="50")] = "size = 50"
error_bar_data$data_size[which(error_bar_data$data_size=="100")] = "size = 100"
error_bar_data$data_size = factor(error_bar_data$data_size, levels=c("size = 10","size = 50","size = 100"))
ggplot() + geom_errorbar(data = error_bar_data, aes(x = P, ymin = n_min,
                                                    ymax = n_max, color = parameter), position = position_dodge(0.75)) +
  geom_hline(yintercept = 10,  linetype="dotted", color= "black") +  facet_grid( data_size ~ .)+ theme(text = element_text(size=20)) + ylab("n")+ 
  theme(text = element_text(size=30), axis.text.x = element_text(angle = 90), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20))


