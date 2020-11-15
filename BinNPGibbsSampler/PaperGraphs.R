#####################################################
######### Fiducial Binomial Simulation Study ########
############## n unknown, p unknown #################
#####################################################
########## Author: Alexander Murph, UNC #############
############## November 15th, 2020 ##################
#####################################################

### This script has the code for the Bin n,p unknown graphs in the paper.

library(ggplot2)
library(ggExtra)
library(gridExtra)

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
x = read.csv("FinalPosteriors/all_data0151.csv")
x = x[which(x$iteration_number %in% floor(number_of_iterations/2):number_of_iterations),]
error_bar_data = x
error_bar_data$lower_mu = error_bar_data$lower_p*error_bar_data$n_value
error_bar_data$upper_mu = error_bar_data$upper_p*error_bar_data$n_value
# error_bar_data$color_scale = as.numeric(error_bar_data$iteration_number - 
                                          # floor(number_of_iterations/2))/floor(number_of_iterations/2)
n_hat = read.csv("FinalPosteriors/first_n0151.csv")$x
first_mu = read.csv("FinalPosteriors/first_mu0151.csv")$x
n = 15
p = 0.1
number_of_iterations = 10000
Mid = floor(3*number_of_iterations/4)
pmainTest = ggplot(data = error_bar_data, aes(x = n_value, ymin = lower_mu,
                                              ymax = upper_mu, y = (lower_mu + upper_mu)/2,
                                              color = iteration_number)) +
  coord_flip() +
  geom_errorbar(position = position_dodge(1)) + ylab(TeX("$\\mu$")) +
  scale_color_gradient2(midpoint = Mid, low="white", mid = "skyblue", high="navyblue") +
  theme(legend.position = "none", text = element_text(size=20)) + 
  geom_point(size = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"))+
  geom_vline(xintercept= c(n_hat, n), color = c("darkorchid2", "red")) +
  geom_hline(yintercept=c(first_mu, p*n), color = c("darkorchid2", "red"))
pmainTest

range_mus = c()
range_ns = c()
error_bar_data$iteration_number = as.factor(error_bar_data$iteration_number)
for(value in unique(error_bar_data$iteration_number)){
  temp_data = error_bar_data[which(error_bar_data$iteration_number == value),]
  max_n = max(temp_data$n_value)
  min_n = min(temp_data$n_value)
  max_mu = max(temp_data$upper_mu)
  min_mu = min(temp_data$lower_mu)
  if(max_n > 1500){
    max_n = Inf
  }
  
  if(max_n == min_n){
    range_ns = c(range_ns, max_n)
  } else {
    range_ns = c(range_ns, min_n, max_n)
  }
  if(max_mu == min_mu){
    range_mus = c(range_mus, max_mu)
  } else {
    range_mus = c(range_mus, min_mu, max_mu)
  }
}


# Chart creation
# We want this density chart to be on top margin
min_mu = layer_scales(pmainTest)$y$range$range[1]
max_mu = layer_scales(pmainTest)$y$range$range[2]
min_n = layer_scales(pmainTest)$x$range$range[1]
max_n = layer_scales(pmainTest)$x$range$range[2]

confidence_curve_data = data.frame(mu = seq(from = min_mu, to = max_mu, length.out = 100 ))

post_dist_mu = ecdf(range_mus)

conf.value = function(mu) {
  diff = abs( 1/2 - post_dist_mu(mu)) * 2
  return(diff)
}
confidence_curve_data$Confidence = sapply(confidence_curve_data$mu, conf.value)

top_density <- ggplot(confidence_curve_data, aes(x = mu, y = Confidence)) +
  geom_line( color = "blue") +
  scale_x_continuous(expand = c(0, 0)) +
  expand_limits(x = c(min_mu, max_mu)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"), text = element_text(size=20)) +
  ylab('') +
  xlab('')


# Empty chart. If you place a chart here, it will appear in the top right
# margin corner, on top of the box plot.
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank())


confidence_curve_data = data.frame(n_val = seq(from = min_n, to = max_n, length.out = 100 ))
post_dist_mu = ecdf(range_ns)
conf.value = function(mu) {
  diff = abs( 1/2 - post_dist_mu(mu)) * 2
  return(diff)
}
confidence_curve_data$Confidence = sapply(confidence_curve_data$n, conf.value)

# We want this boxplot to be on the right margin
right_boxplot <- ggplot(confidence_curve_data, aes(y = Confidence, x = n_val)) +
  geom_line( color = "blue") +
  coord_flip() +
  theme(text = element_text(size=20)) +
  scale_x_continuous(expand = c(0, 0)) +
  expand_limits(x = c(min_n, max_n)) +
  ylab('') +
  xlab('')


grid.arrange(top_density, empty, pmainTest, right_boxplot, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))



###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
x = read.csv("FinalPosteriors/all_data0751.csv")
x = x[which(x$iteration_number %in% floor(number_of_iterations/2):number_of_iterations),]
error_bar_data = x
error_bar_data$lower_mu = error_bar_data$lower_p*error_bar_data$n_value
error_bar_data$upper_mu = error_bar_data$upper_p*error_bar_data$n_value
# error_bar_data$color_scale = as.numeric(error_bar_data$iteration_number - 
# floor(number_of_iterations/2))/floor(number_of_iterations/2)
n_hat = read.csv("FinalPosteriors/first_n0751.csv")$x
first_mu = read.csv("FinalPosteriors/first_mu0751.csv")$x
n = 75
p = 0.1
number_of_iterations = 10000
Mid = floor(3*number_of_iterations/4)
pmainTest = ggplot(data = error_bar_data, aes(x = n_value, ymin = lower_mu,
                                              ymax = upper_mu, y = (lower_mu + upper_mu)/2,
                                              color = iteration_number)) +
  coord_flip() +
  geom_errorbar(position = position_dodge(1)) + ylab(TeX("$\\mu$")) +
  scale_color_gradient2(midpoint = Mid, low="white", mid = "skyblue", high="navyblue") +
  theme(legend.position = "none") + 
  geom_point(size = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"))+
  geom_vline(xintercept= c(n_hat, n), color = c("darkorchid2", "red")) +
  geom_hline(yintercept=c(first_mu, p*n), color = c("darkorchid2", "red"))
pmainTest

range_mus = c()
range_ns = c()
error_bar_data$iteration_number = as.factor(error_bar_data$iteration_number)
for(value in unique(error_bar_data$iteration_number)){
  temp_data = error_bar_data[which(error_bar_data$iteration_number == value),]
  max_n = max(temp_data$n_value)
  min_n = min(temp_data$n_value)
  max_mu = max(temp_data$upper_mu)
  min_mu = min(temp_data$lower_mu)
  if(max_n > 1500){
    max_n = Inf
  }
  
  if(max_n == min_n){
    range_ns = c(range_ns, max_n)
  } else {
    range_ns = c(range_ns, min_n, max_n)
  }
  if(max_mu == min_mu){
    range_mus = c(range_mus, max_mu)
  } else {
    range_mus = c(range_mus, min_mu, max_mu)
  }
}


# Chart creation
# We want this density chart to be on top margin
min_mu = layer_scales(pmainTest)$y$range$range[1]
max_mu = layer_scales(pmainTest)$y$range$range[2]
min_n = layer_scales(pmainTest)$x$range$range[1]
max_n = layer_scales(pmainTest)$x$range$range[2]

confidence_curve_data = data.frame(mu = seq(from = min_mu, to = max_mu, length.out = 100 ))

post_dist_mu = ecdf(range_mus)

conf.value = function(mu) {
  diff = abs( 1/2 - post_dist_mu(mu)) * 2
  return(diff)
}
confidence_curve_data$Confidence = sapply(confidence_curve_data$mu, conf.value)

top_density <- ggplot(confidence_curve_data, aes(x = mu, y = Confidence)) +
  geom_line( color = "blue") +
  scale_x_continuous(expand = c(0, 0)) +
  expand_limits(x = c(min_mu, max_mu)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines")) +
  ylab('') +
  xlab('')


# Empty chart. If you place a chart here, it will appear in the top right
# margin corner, on top of the box plot.
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank())


confidence_curve_data = data.frame(n_val = seq(from = min_n, to = max_n, length.out = 100 ))
post_dist_mu = ecdf(range_ns)
conf.value = function(mu) {
  diff = abs( 1/2 - post_dist_mu(mu)) * 2
  return(diff)
}
confidence_curve_data$Confidence = sapply(confidence_curve_data$n, conf.value)

# We want this boxplot to be on the right margin
right_boxplot <- ggplot(confidence_curve_data, aes(y = Confidence, x = n_val)) +
  geom_line( color = "blue") +
  coord_flip() +
  scale_x_continuous(expand = c(0, 0)) +
  expand_limits(x = c(min_n, max_n)) +
  ylab('') +
  xlab('')


grid.arrange(top_density, empty, pmainTest, right_boxplot, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))




###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
x = read.csv("FinalPosteriors/all_data2159.csv")
x = x[which(x$iteration_number %in% floor(number_of_iterations/2):number_of_iterations),]
error_bar_data = x
error_bar_data$lower_mu = error_bar_data$lower_p*error_bar_data$n_value
error_bar_data$upper_mu = error_bar_data$upper_p*error_bar_data$n_value
# error_bar_data$color_scale = as.numeric(error_bar_data$iteration_number - 
# floor(number_of_iterations/2))/floor(number_of_iterations/2)
n_hat = read.csv("FinalPosteriors/first_n2159.csv")$x
first_mu = read.csv("FinalPosteriors/first_mu2159.csv")$x
n = 15
p = 0.9
number_of_iterations = 10000
Mid = floor(3*number_of_iterations/4)

error_bar_data$iteration_number = as.factor(error_bar_data$iteration_number)
pale = colorRampPalette(c("white", "navyblue"))( nrow(error_bar_data))
colScale <- scale_colour_manual(name = "iteration_number",values = pale)

pmainTest = ggplot(data = error_bar_data, aes(x = n_value, ymin = lower_mu,
                                              ymax = upper_mu, y = (lower_mu + upper_mu)/2,
                                              color = iteration_number)) +
  coord_flip() +
  geom_errorbar(position=position_dodge(1)) + ylab(TeX("$\\mu$")) +
  theme(legend.position = "none") + 
  geom_point(size = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  colScale +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"))+
  geom_vline(xintercept= c(n_hat, n), color = c("darkorchid2", "red")) +
  geom_hline(yintercept=c(first_mu, p*n), color = c("darkorchid2", "red"))
pmainTest

range_mus = c()
range_ns = c()
error_bar_data$iteration_number = as.factor(error_bar_data$iteration_number)
for(value in unique(error_bar_data$iteration_number)){
  temp_data = error_bar_data[which(error_bar_data$iteration_number == value),]
  max_n = max(temp_data$n_value)
  min_n = min(temp_data$n_value)
  max_mu = max(temp_data$upper_mu)
  min_mu = min(temp_data$lower_mu)
  if(max_n > 1500){
    max_n = Inf
  }
  
  if(max_n == min_n){
    range_ns = c(range_ns, max_n)
  } else {
    range_ns = c(range_ns, min_n, max_n)
  }
  if(max_mu == min_mu){
    range_mus = c(range_mus, max_mu)
  } else {
    range_mus = c(range_mus, min_mu, max_mu)
  }
}


# Chart creation
# We want this density chart to be on top margin
min_mu = layer_scales(pmainTest)$y$range$range[1]
max_mu = layer_scales(pmainTest)$y$range$range[2]
min_n = layer_scales(pmainTest)$x$range$range[1]
max_n = layer_scales(pmainTest)$x$range$range[2]

confidence_curve_data = data.frame(mu = seq(from = min_mu, to = max_mu, length.out = 100 ))

post_dist_mu = ecdf(range_mus)

conf.value = function(mu) {
  diff = abs( 1/2 - post_dist_mu(mu)) * 2
  return(diff)
}
confidence_curve_data$Confidence = sapply(confidence_curve_data$mu, conf.value)

top_density <- ggplot(confidence_curve_data, aes(x = mu, y = Confidence)) +
  geom_line( color = "blue") +
  scale_x_continuous(expand = c(0, 0)) +
  expand_limits(x = c(min_mu, max_mu)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines")) +
  ylab('') +
  xlab('')


# Empty chart. If you place a chart here, it will appear in the top right
# margin corner, on top of the box plot.
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank())


confidence_curve_data = data.frame(n_val = seq(from = min_n, to = max_n, length.out = 100 ))
post_dist_mu = ecdf(range_ns)
conf.value = function(mu) {
  diff = abs( 1/2 - post_dist_mu(mu)) * 2
  return(diff)
}
confidence_curve_data$Confidence = sapply(confidence_curve_data$n, conf.value)

# We want this boxplot to be on the right margin
right_boxplot <- ggplot(confidence_curve_data, aes(y = Confidence, x = n_val)) +
  geom_line( color = "blue") +
  coord_flip() +
  scale_x_continuous(expand = c(0, 0)) +
  expand_limits(x = c(min_n, max_n)) +
  ylab('') +
  xlab('')


grid.arrange(top_density, empty, pmainTest, right_boxplot, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))



###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
x = read.csv("FinalPosteriors/all_data0759.csv")
x = x[which(x$iteration_number %in% floor(number_of_iterations/2):number_of_iterations),]
error_bar_data = x
error_bar_data$lower_mu = error_bar_data$lower_p*error_bar_data$n_value
error_bar_data$upper_mu = error_bar_data$upper_p*error_bar_data$n_value
# error_bar_data$color_scale = as.numeric(error_bar_data$iteration_number - 
# floor(number_of_iterations/2))/floor(number_of_iterations/2)
n_hat = read.csv("FinalPosteriors/first_n0759.csv")$x
first_mu = read.csv("FinalPosteriors/first_mu0759.csv")$x
n = 75
p = 0.9
number_of_iterations = 10000
Mid = floor(3*number_of_iterations/4)

error_bar_data$iteration_number = as.factor(error_bar_data$iteration_number)
pale = colorRampPalette(c("white", "navyblue"))( nrow(error_bar_data))
colScale <- scale_colour_manual(name = "iteration_number",values = pale)

pmainTest = ggplot(data = error_bar_data, aes(x = n_value, ymin = lower_mu,
                                              ymax = upper_mu, y = (lower_mu + upper_mu)/2,
                                              color = iteration_number)) +
  coord_flip() +
  geom_errorbar(position=position_dodge(1)) + ylab(TeX("$\\mu$")) +
  theme(legend.position = "none") + 
  geom_point(size = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  colScale +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"))+
  geom_vline(xintercept= c(n_hat, n), color = c("darkorchid2", "red")) +
  geom_hline(yintercept=c(first_mu, p*n), color = c("darkorchid2", "red"))
pmainTest

range_mus = c()
range_ns = c()
error_bar_data$iteration_number = as.factor(error_bar_data$iteration_number)
for(value in unique(error_bar_data$iteration_number)){
  temp_data = error_bar_data[which(error_bar_data$iteration_number == value),]
  max_n = max(temp_data$n_value)
  min_n = min(temp_data$n_value)
  max_mu = max(temp_data$upper_mu)
  min_mu = min(temp_data$lower_mu)
  if(max_n > 1500){
    max_n = Inf
  }
  
  if(max_n == min_n){
    range_ns = c(range_ns, max_n)
  } else {
    range_ns = c(range_ns, min_n, max_n)
  }
  if(max_mu == min_mu){
    range_mus = c(range_mus, max_mu)
  } else {
    range_mus = c(range_mus, min_mu, max_mu)
  }
}


# Chart creation
# We want this density chart to be on top margin
min_mu = layer_scales(pmainTest)$y$range$range[1]
max_mu = layer_scales(pmainTest)$y$range$range[2]
min_n = layer_scales(pmainTest)$x$range$range[1]
max_n = layer_scales(pmainTest)$x$range$range[2]

confidence_curve_data = data.frame(mu = seq(from = min_mu, to = max_mu, length.out = 100 ))

post_dist_mu = ecdf(range_mus)

conf.value = function(mu) {
  diff = abs( 1/2 - post_dist_mu(mu)) * 2
  return(diff)
}
confidence_curve_data$Confidence = sapply(confidence_curve_data$mu, conf.value)

top_density <- ggplot(confidence_curve_data, aes(x = mu, y = Confidence)) +
  geom_line( color = "blue") +
  scale_x_continuous(expand = c(0, 0)) +
  expand_limits(x = c(min_mu, max_mu)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines")) +
  ylab('') +
  xlab('')


# Empty chart. If you place a chart here, it will appear in the top right
# margin corner, on top of the box plot.
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank())


confidence_curve_data = data.frame(n_val = seq(from = min_n, to = max_n, length.out = 100 ))
post_dist_mu = ecdf(range_ns)
conf.value = function(mu) {
  diff = abs( 1/2 - post_dist_mu(mu)) * 2
  return(diff)
}
confidence_curve_data$Confidence = sapply(confidence_curve_data$n, conf.value)

# We want this boxplot to be on the right margin
right_boxplot <- ggplot(confidence_curve_data, aes(y = Confidence, x = n_val)) +
  geom_line( color = "blue") +
  coord_flip() +
  scale_x_continuous(expand = c(0, 0)) +
  expand_limits(x = c(min_n, max_n)) +
  ylab('') +
  xlab('')


grid.arrange(top_density, empty, pmainTest, right_boxplot, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

