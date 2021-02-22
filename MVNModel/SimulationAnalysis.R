# Analysing Covariance Simulations
library(ggplot2)
library(grid)

# Identity covariance matrix:
full_data = data.frame(value = NA, description=NA,iteration=NA)

all_acceptance_ratios = c()
for(name in list.files(path = "Data/")[-1]){
  temp_data = read.csv(paste("Data/",name, sep =""))[,-1]
  # temp_data = temp_data[complete.cases(temp_data),]
  full_data = rbind(full_data, temp_data)
}
full_data = full_data[-1,]

# Checking the QQ plots:
temp_data = full_data[which(full_data$description%in%c("QQ_loss1",
                                                       "QQ_FM","QQ_LogDet",
                                                       "QQ_SingularMax","QQ_Frobenius",
                                                       "QQ_Spectral_Diff","QQ_Frobenius_Diff") ),]
ggplot(temp_data, aes(x=value,color=description)) + stat_ecdf() + 
  scale_color_discrete(name = "description")+
  xlab("Nominal Coverage")+ylab("Empirical Coverage")+ 
  scale_color_discrete(name = "description", labels = c("FM Distance","Frobenius Parameter", "Frobenius Distance",
                                                        "LogDet Parameter", "Stein's Loss", "Spectral Norm Parameter",
                                                        "Spectral Norm Distance")) +theme(text = element_text(size=20))
temp_data = full_data[which(full_data$description%in%c("QQ_Mu") ),]
ggplot(temp_data, aes(x=value,color=description)) + stat_ecdf() + 
  scale_color_discrete(name = "description", labels = c("Upper CI on Mu"))+
  xlab("Nominal Coverage")+ylab("Empirical Coverage")+theme(text = element_text(size=20))


temp_data = full_data[which(full_data$description%in%c("Loss1_95CI_ContainsTruth",
                                                       "FM_95CI_ContainsTruth","LogDet_95CI_ContainsTruth",
                                                       "SingularMax_95CI_ContainsTruth", "Frobenius_95CI_ContainsTruth",
                                                       "Spectral_Diff_95CI_ContainsTruth", "Frobenius_Diff_95CI_ContainsTruth") ),]
tapply(temp_data$value, temp_data$description, summary)

