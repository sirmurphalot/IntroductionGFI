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
library(latex2exp)

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
x = read.csv("FinalPosteriors/all_data0151.csv")
n = 15
p = 0.1
number_of_iterations = 10000
Mid = floor(3*number_of_iterations/4)
x = x[which(x$iteration_number %in% floor(number_of_iterations/2):number_of_iterations),]
error_bar_data = x
error_bar_data$lower_mu = error_bar_data$lower_p*error_bar_data$n_value
error_bar_data$upper_mu = error_bar_data$upper_p*error_bar_data$n_value
# error_bar_data$color_scale = as.numeric(error_bar_data$iteration_number - 
                                          # floor(number_of_iterations/2))/floor(number_of_iterations/2)
n_hat = read.csv("FinalPosteriors/first_n0151.csv")$x
first_mu = read.csv("FinalPosteriors/first_mu0151.csv")$x
pmainTest = ggplot(data = error_bar_data, aes(x = n_value, ymin = lower_mu,
                                              ymax = upper_mu, y = (lower_mu + upper_mu)/2,
                                              color = iteration_number)) +
  coord_flip() +
  geom_errorbar(position = position_dodge(1)) + xlab(TeX("$n$")) + ylab(TeX("$\\mu$")) +
  scale_color_gradient2(midpoint = Mid, low="white", mid = "skyblue", high="navyblue") +
  theme(legend.position = "none") + 
  geom_point(size = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"))+
  geom_vline(xintercept= c(n_hat, n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1)) +
  geom_hline(yintercept=c(first_mu, p*n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))+ 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))
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
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"), text = element_text(size=20), axis.text.x = element_text(size=20)) +
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

# Updating this 05/04/2021 wrt Jan's note from 04/30/2021:
# conf.value = function(mu) {
#   diff = abs( 1/2 - post_dist_mu(mu)) * 2
#   return(diff)
# }
conf.value = function(n) {
  less_or_equal = mean(range_ns<=n)
  greater_or_equal = mean(range_ns>=n)
  diff = 2*(max(c(less_or_equal,greater_or_equal))-0.5)
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
  xlab('')+ 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), axis.text.x = element_text(angle=270,size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))


grid.arrange(top_density, empty, pmainTest, right_boxplot, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

## Marginal Plot
summary(error_bar_data)
marginal_ns = NULL
marginal_mus = NULL
for(i in unique(error_bar_data$iteration_number)){
  draw_indices = which(error_bar_data$iteration_number == i)
  mus = c( min(error_bar_data$lower_mu[draw_indices]), max(error_bar_data$upper_mu[draw_indices]))
  ns = c( min(error_bar_data$n_value[draw_indices]), max(error_bar_data$n_value[draw_indices]))
  
  bern_1 = 2 - rbinom(1,size=1,prob=0.5)
  bern_2 = 2 - rbinom(1,size=1,prob=0.5)
  
  marginal_mus = c(marginal_mus,mus[bern_1])
  marginal_ns = c(marginal_ns,ns[bern_2])
}

graph_data = data.frame(x=5000:10000,mus=marginal_mus)
p1 = ggplot(graph_data,aes(x=x,y=mus)) + geom_line(color="turquoise3") +
  theme(legend.position = "none", text = element_text(size=25)) + xlab("")+ ylab(TeX("$\\mu$")) +
  geom_hline(yintercept=c(first_mu, p*n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))
graph_data2 = data.frame(x=5000:10000,ns=marginal_ns)
p2 = ggplot(graph_data2,aes(x=x,y=ns)) + geom_line(color="darkblue") +
  theme(legend.position = "none", text = element_text(size=25)) + xlab("iteration") + ylab("n") +
  geom_hline(yintercept= c(n_hat, n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))
grid.arrange(p1,p2, ncol=1, nrow=2)


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
n = 75
p = 0.1
number_of_iterations = 10000
x = read.csv("FinalPosteriors/all_data0751.csv")
x = x[which(x$iteration_number %in% floor(number_of_iterations/2):number_of_iterations),]
error_bar_data = x
error_bar_data$lower_mu = error_bar_data$lower_p*error_bar_data$n_value
error_bar_data$upper_mu = error_bar_data$upper_p*error_bar_data$n_value
# error_bar_data$color_scale = as.numeric(error_bar_data$iteration_number - 
# floor(number_of_iterations/2))/floor(number_of_iterations/2)
n_hat = read.csv("FinalPosteriors/first_n0751.csv")$x
first_mu = read.csv("FinalPosteriors/first_mu0751.csv")$x

Mid = floor(3*number_of_iterations/4)
pmainTest = ggplot(data = error_bar_data, aes(x = n_value, ymin = lower_mu,
                                              ymax = upper_mu, y = (lower_mu + upper_mu)/2,
                                              color = iteration_number)) +
  coord_flip() +
  geom_errorbar(position = position_dodge(1)) + xlab(TeX("$n$")) + ylab(TeX("$\\mu$")) +
  scale_color_gradient2(midpoint = Mid, low="white", mid = "skyblue", high="navyblue") +
  theme(legend.position = "none") +
  geom_point(size = -1) +
  labs( color ="Iteration Number\n")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"),
        legend.key.height = unit(1, "inch"))+
  geom_vline(xintercept= c(n_hat, n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1)) +
  geom_hline(yintercept=c(first_mu, p*n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))+ 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=30))
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
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"), text = element_text(size=20), axis.text.x = element_text(size=20)) +
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
# Updating this 05/04/2021 wrt Jan Hannig's note from 04/30/2021:
# conf.value = function(mu) {
#   diff = abs( 1/2 - post_dist_mu(mu)) * 2
#   return(diff)
# }
conf.value = function(n) {
  less_or_equal = mean(range_ns<=n)
  greater_or_equal = mean(range_ns>=n)
  diff = 2*(max(c(less_or_equal,greater_or_equal))-0.5)
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
  xlab('')+ 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), axis.text.x = element_text(angle=270,size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))


grid.arrange(top_density, empty, pmainTest, right_boxplot, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

## Marginal Plot
summary(error_bar_data)
marginal_ns = NULL
marginal_mus = NULL
for(i in unique(error_bar_data$iteration_number)){
  draw_indices = which(error_bar_data$iteration_number == i)
  mus = c( min(error_bar_data$lower_mu[draw_indices]), max(error_bar_data$upper_mu[draw_indices]))
  ns = c( min(error_bar_data$n_value[draw_indices]), max(error_bar_data$n_value[draw_indices]))
  
  bern_1 = 2 - rbinom(1,size=1,prob=0.5)
  bern_2 = 2 - rbinom(1,size=1,prob=0.5)
  
  marginal_mus = c(marginal_mus,mus[bern_1])
  marginal_ns = c(marginal_ns,ns[bern_2])
}

graph_data = data.frame(x=5000:10000,mus=marginal_mus)
p1 = ggplot(graph_data,aes(x=x,y=mus)) + geom_line(color="turquoise3") +
  theme(legend.position = "none", text = element_text(size=25)) + xlab("")+ ylab(TeX("$\\mu$")) +
  geom_hline(yintercept=c(first_mu, p*n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))
graph_data2 = data.frame(x=5000:10000,ns=marginal_ns)
p2 = ggplot(graph_data2,aes(x=x,y=ns)) + geom_line(color="darkblue") +
  theme(legend.position = "none", text = element_text(size=25)) + xlab("iteration") + ylab("n") +
  geom_hline(yintercept= c(n_hat, n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))
grid.arrange(p1,p2, ncol=1, nrow=2)


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
n = 15
p = 0.9
number_of_iterations = 10000
x = read.csv("FinalPosteriors/all_data2159.csv")
x = x[which(x$iteration_number %in% floor(number_of_iterations/2):number_of_iterations),]
error_bar_data = x
error_bar_data$lower_mu = error_bar_data$lower_p*error_bar_data$n_value
error_bar_data$upper_mu = error_bar_data$upper_p*error_bar_data$n_value
# error_bar_data$color_scale = as.numeric(error_bar_data$iteration_number - 
# floor(number_of_iterations/2))/floor(number_of_iterations/2)
n_hat = read.csv("FinalPosteriors/first_n2159.csv")$x
first_mu = read.csv("FinalPosteriors/first_mu2159.csv")$x

Mid = floor(3*number_of_iterations/4)

error_bar_data$iteration_number = as.factor(error_bar_data$iteration_number)
pale = colorRampPalette(c("white", "navyblue"))( nrow(error_bar_data))
colScale <- scale_colour_manual(name = "iteration_number",values = pale)

pmainTest = ggplot(data = error_bar_data, aes(x = n_value, ymin = lower_mu,
                                              ymax = upper_mu, y = (lower_mu + upper_mu)/2,
                                              color = iteration_number)) +
  coord_flip() +
  geom_errorbar(position=position_dodge(1)) + ylab(TeX("$\\mu$")) + xlab(TeX("$n$")) +
  theme(legend.position = "none") +
  geom_point(size = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  colScale +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"))+
  geom_vline(xintercept= c(n_hat, n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1)) +
  geom_hline(yintercept=c(first_mu, p*n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))+ 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))
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
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"), text = element_text(size=20), axis.text.x = element_text(size=20)) +
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

# Updating this 05/04/2021 wrt Jan Hannig's note from 04/30/2021:
# conf.value = function(mu) {
#   diff = abs( 1/2 - post_dist_mu(mu)) * 2
#   return(diff)
# }
conf.value = function(n) {
  less_or_equal = mean(range_ns<=n)
  greater_or_equal = mean(range_ns>=n)
  diff = 2*(max(c(less_or_equal,greater_or_equal))-0.5)
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
  xlab('') + 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), axis.text.x = element_text(angle=270,size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))



grid.arrange(top_density, empty, pmainTest, right_boxplot, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

## Marginal Plot
summary(error_bar_data)
marginal_ns = NULL
marginal_mus = NULL
for(i in unique(error_bar_data$iteration_number)){
  draw_indices = which(error_bar_data$iteration_number == i)
  mus = c( min(error_bar_data$lower_mu[draw_indices]), max(error_bar_data$upper_mu[draw_indices]))
  ns = c( min(error_bar_data$n_value[draw_indices]), max(error_bar_data$n_value[draw_indices]))
  
  bern_1 = 2 - rbinom(1,size=1,prob=0.5)
  bern_2 = 2 - rbinom(1,size=1,prob=0.5)
  
  marginal_mus = c(marginal_mus,mus[bern_1])
  marginal_ns = c(marginal_ns,ns[bern_2])
}

graph_data = data.frame(x=5000:10000,mus=marginal_mus)
p1 = ggplot(graph_data,aes(x=x,y=mus)) + geom_line(color="turquoise3") +
  theme(legend.position = "none", text = element_text(size=25)) + xlab("")+ ylab(TeX("$\\mu$")) +
  geom_hline(yintercept=c(first_mu, p*n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))
graph_data2 = data.frame(x=5000:10000,ns=marginal_ns)
p2 = ggplot(graph_data2,aes(x=x,y=ns)) + geom_line(color="darkblue") +
  theme(legend.position = "none", text = element_text(size=25)) + xlab("iteration") + ylab("n") +
  geom_hline(yintercept= c(n_hat, n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))
grid.arrange(p1,p2, ncol=1, nrow=2)


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
x = read.csv("FinalPosteriors/all_data0759.csv")
n = 75
p = 0.9
number_of_iterations = 10000
x = x[which(x$iteration_number %in% floor(number_of_iterations/2):number_of_iterations),]
error_bar_data = x
error_bar_data$lower_mu = error_bar_data$lower_p*error_bar_data$n_value
error_bar_data$upper_mu = error_bar_data$upper_p*error_bar_data$n_value
# error_bar_data$color_scale = as.numeric(error_bar_data$iteration_number - 
# floor(number_of_iterations/2))/floor(number_of_iterations/2)
n_hat = read.csv("FinalPosteriors/first_n0759.csv")$x
first_mu = read.csv("FinalPosteriors/first_mu0759.csv")$x

Mid = floor(3*number_of_iterations/4)

error_bar_data$iteration_number = as.factor(error_bar_data$iteration_number)
pale = colorRampPalette(c("white", "navyblue"))( nrow(error_bar_data))
colScale <- scale_colour_manual(name = "iteration_number",values = pale)

pmainTest = ggplot(data = error_bar_data, aes(x = n_value, ymin = lower_mu,
                                              ymax = upper_mu, y = (lower_mu + upper_mu)/2,
                                              color = iteration_number)) +
  coord_flip() +
  geom_errorbar(position=position_dodge(1)) + ylab(TeX("$\\mu$")) + xlab(TeX("$n$")) +
  theme(legend.position = "none") +
  geom_point(size = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  colScale +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"))+
  geom_vline(xintercept= c(n_hat, n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1)) +
  geom_hline(yintercept=c(first_mu, p*n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))+ 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))
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
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"), text = element_text(size=20), axis.text.x = element_text(size=20)) +
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
# Updating this 05/04/2021 wrt Jan Hannig's note from 04/30/2021:
# conf.value = function(mu) {
#   diff = abs( 1/2 - post_dist_mu(mu)) * 2
#   return(diff)
# }
conf.value = function(n) {
  less_or_equal = mean(range_ns<=n)
  greater_or_equal = mean(range_ns>=n)
  diff = 2*(max(c(less_or_equal,greater_or_equal))-0.5)
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
  xlab('') + 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), axis.text.x = element_text(angle=270,size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))



grid.arrange(top_density, empty, pmainTest, right_boxplot, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))


## Marginal Plot
summary(error_bar_data)
marginal_ns = NULL
marginal_mus = NULL
for(i in unique(error_bar_data$iteration_number)){
  draw_indices = which(error_bar_data$iteration_number == i)
  mus = c( min(error_bar_data$lower_mu[draw_indices]), max(error_bar_data$upper_mu[draw_indices]))
  ns = c( min(error_bar_data$n_value[draw_indices]), max(error_bar_data$n_value[draw_indices]))
  
  bern_1 = 2 - rbinom(1,size=1,prob=0.5)
  bern_2 = 2 - rbinom(1,size=1,prob=0.5)
  
  marginal_mus = c(marginal_mus,mus[bern_1])
  marginal_ns = c(marginal_ns,ns[bern_2])
}

graph_data = data.frame(x=5000:10000,mus=marginal_mus)
p1 = ggplot(graph_data,aes(x=x,y=mus)) + geom_line(color="turquoise3") +
  theme(legend.position = "none", text = element_text(size=25)) + xlab("")+ ylab(TeX("$\\mu$")) +
  geom_hline(yintercept=c(first_mu, p*n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))
graph_data2 = data.frame(x=5000:10000,ns=marginal_ns)
p2 = ggplot(graph_data2,aes(x=x,y=ns)) + geom_line(color="darkblue") +
  theme(legend.position = "none", text = element_text(size=25)) + xlab("iteration") + ylab("n") +
  geom_hline(yintercept= c(n_hat, n), linetype = c("dotted", "longdash"), color = c("red","red"), size = c(1,1))
grid.arrange(p1,p2, ncol=1, nrow=2)



###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
### Example mu vs p parameterization figure
## p parameterization
error_bar_data = data.frame(upper_mu = -50*((1:30)/1000)^2 + 0.3, 
                            lower_mu = 50*((1:30)/1000)^2 + 0.1, 
                            n_value = rep(1:30),
                            iteration_number = rep(1,times=30))
error_bar_data$lower_p = error_bar_data$lower_mu / error_bar_data$n_value
error_bar_data$upper_p = error_bar_data$upper_mu / error_bar_data$n_value

error_bar_data$iteration_number = as.factor(error_bar_data$iteration_number)

pmainTest = ggplot() +
  geom_errorbar(data = error_bar_data,
                mapping = aes(x = n_value, ymin = lower_p,
                              ymax = upper_p, y = (lower_p + upper_p)/2),color = "blue") +
  ylab(TeX("$p$")) +
  xlab(TeX("$n$")) + ylim(0,0.4) + xlim(0,33) +
  coord_flip() +
  geom_point(size = -1) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines")) + 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))
pmainTest

## mu parameterization
error_bar_data = data.frame(upper_mu = -50*((1:30)/1000)^2 + 0.3, 
                            lower_mu = 50*((1:30)/1000)^2 + 0.1, 
                            n_value = rep(1:30),
                            iteration_number = rep(1,times=30))

error_bar_data$iteration_number = as.factor(error_bar_data$iteration_number)

pmainTest = ggplot() +
  geom_errorbar(data = error_bar_data,
                mapping = aes(x = n_value, ymin = lower_mu,
                              ymax = upper_mu, y = (lower_mu + upper_mu)/2),color = "blue") +
  ylab(TeX("$\\mu$")) +
  xlab(TeX("$n$")) + ylim(0,0.4) + xlim(0,33) +
  coord_flip() +
  geom_point(size = -1) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines")) + 
  theme(text = element_text(size=30), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))
pmainTest


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
##################################
## Large Start Stuff

x = read.csv("Datalogs/large_start_all_data0759.csv")
n = 75
p = 0.9
number_of_iterations = 10000
# x = x[which(x$iteration_number %in% floor(number_of_iterations/2):number_of_iterations),]
error_bar_data = x
error_bar_data$lower_mu = error_bar_data$lower_p*error_bar_data$n_value
error_bar_data$upper_mu = error_bar_data$upper_p*error_bar_data$n_value
# error_bar_data$color_scale = as.numeric(error_bar_data$iteration_number - 
# floor(number_of_iterations/2))/floor(number_of_iterations/2)
n_hat = read.csv("Datalogs/large_start_first_n0759.csv")$x
first_mu = read.csv("Datalogs/large_start_first_mu0759.csv")$x

## Marginal Plot
summary(error_bar_data)
marginal_ns = NULL
marginal_mus = NULL
for(i in unique(error_bar_data$iteration_number)){
  draw_indices = which(error_bar_data$iteration_number == i)
  mus = c( min(error_bar_data$lower_mu[draw_indices]), max(error_bar_data$upper_mu[draw_indices]))
  ns = c( min(error_bar_data$n_value[draw_indices]), max(error_bar_data$n_value[draw_indices]))
  
  bern_1 = 2 - rbinom(1,size=1,prob=0.5)
  bern_2 = 2 - rbinom(1,size=1,prob=0.5)
  
  marginal_mus = c(marginal_mus,mus[bern_1])
  marginal_ns = c(marginal_ns,ns[bern_2])
}

graph_data = data.frame(x=1:10000,mus=marginal_mus)
p1 = ggplot(graph_data,aes(x=x,y=mus)) + geom_line(color="turquoise3") +
  theme(legend.position = "none", text = element_text(size=25)) + xlab("Iteration")+ ylab(TeX("$\\mu$")) +
  geom_hline(yintercept=c(first_mu, p*n), color = c("darkorchid2", "red"))+ ggtitle(TeX("Traceplot of $\\mu$"))
p1
graph_data2 = data.frame(x=1:10000,ns=marginal_ns) 
p2 = ggplot(graph_data2,aes(x=x,y=ns)) + geom_line(color="darkblue") +
  theme(legend.position = "none", text = element_text(size=25)) + xlab("Iteration") + ylab("n") +
  geom_hline(yintercept= c(n_hat, n), color = c("darkorchid2", "red")) + ggtitle(TeX("Traceplot of $n$"))
p2


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
# Marginal Coverages Information for Table:

simulation_data = data.frame(TrueMu = NA,
                             TrueN = NA,
                             TrueP = NA,
                             DataSize = NA,
                             SimulationNumber = NA,
                             pois_belief_lower_bound_mu = NA,
                             pois_belief_upper_bound_mu = NA,
                             pois_belief_lower_bound_n = NA,
                             pois_belief_upper_bound_n = NA,
                             pois_belief_contains_truth = NA,
                             pois_plausability_lower_bound_mu = NA,
                             pois_plausability_upper_bound_mu = NA,
                             pois_plausability_lower_bound_n = NA,
                             pois_plausability_upper_bound_n = NA,
                             pois_plausability_contains_truth = NA,
                             nopois_belief_lower_bound_mu = NA,
                             nopois_belief_upper_bound_mu = NA,
                             nopois_belief_lower_bound_n = NA,
                             nopois_belief_upper_bound_n = NA,
                             nopois_belief_contains_truth = NA,
                             nopois_plausability_lower_bound_mu = NA,
                             nopois_plausability_upper_bound_mu = NA,
                             nopois_plausability_lower_bound_n = NA,
                             nopois_plausability_upper_bound_n = NA,
                             nopois_plausability_contains_truth = NA,
                             n_upper_CI = NA,
                             n_lower_CI = NA,
                             mu_upper_CI = NA,
                             mu_lower_CI = NA)

#### SET PARAMETERS
p_values = c(0.1, 0.5, 0.9)
data_sizes = c(100)
true_ns  = c(15, 75)

for(name in list.files(path = "BinData/")){
  new_rows = read.csv(paste("BinData/", name, sep = ""))
  new_rows$X = NULL
  new_rows = new_rows[-1,]
  simulation_data = rbind(simulation_data, new_rows)
}
simulation_data = simulation_data[-1,]
simulation_data$nopois_plausability_contains_truth = NULL
simulation_data$nopois_plausability_lower_bound_mu = NULL
simulation_data$nopois_plausability_upper_bound_mu = NULL
simulation_data$nopois_plausability_lower_bound_n = NULL
simulation_data$nopois_plausability_upper_bound_n = NULL
for(p_val in p_values){
  for(n_val in true_ns){
    simulation_data_temp = simulation_data[which( (simulation_data$TrueP==p_val) &
                                                    (simulation_data$TrueN==n_val) ),]
    print(paste("(",p_val,",",n_val,") mu marginal coverage: ", mean((simulation_data_temp$mu_lower_CI<=0.975)&(simulation_data_temp$mu_upper_CI<=0.975)) ))
    print(paste("(",p_val,",",n_val,") n marginal coverage: ", mean(simulation_data_temp$n_lower_CI<=0.95) ))
    print(paste("(",p_val,",",n_val,") belief coverage: ", mean(simulation_data_temp$pois_belief_contains_truth) ))
    print(paste("(",p_val,",",n_val,") plausability coverage: ", mean(simulation_data_temp$pois_plausability_contains_truth) ))
  }
}
