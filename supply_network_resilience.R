library(igraph)
library(NetSim)
library(ggplot2)
library(tidyr)
#library(readr)
#library(dplyr)

source('network_construction.r')
source('disruption_status.r')

sim_round = 200
tier = 5
#sup_num_tier = c(1,round(runif(tier, 3, 6),0))
sup_num_tier = c(1, 10, 10, 10, 10, 10)
p_edge = 0.3
p_edge_tier = rep(p_edge, tier)  # prob can be different in diff tiers
same_tier=0
p_edge_same_tier = 0.05  # also can be different
direction=T
exam_t = 10   # the examination cycle, can be 1
disruption_times = 3 # how many disruptions happen -- the first always happens at the beginning
disrupt_t = c(1,3,5) # a user-specified disruption time points...
prct_node_down = 0.35   # percentage nodes down
node_down_to_capacity = 0.2   # only a portion of products will go through...
recovery_node_each_time = 0.1
#prct_edge_down = 0.3
#edge_down_to_capacity = 0.7
#recovery_edge_each_time = 0.1

focal_prod_matrix = matrix(NA, sim_round, exam_t)
colnames(focal_prod_matrix) = paste0("t", 1:exam_t)

written_row = 1
sim_now = 1

while (written_row <= sim_round){
  
  cat("this is sim", sim_now, ", writing row", written_row, "\n")
  
  ######### Section 1. Setup and Simulate networks #########
  
  seed = sim_now * 100 + 12345
  set.seed(seed)
  
  sn = network_construction(tier=tier, sup_num_tier=sup_num_tier, p_edge_tier=p_edge_tier, same_tier=same_tier,
                            p_edge_same_tier=p_edge_same_tier, direction=direction)
  
  ######### Section 2. Dynamic disruptions #########
  
  network = sn$sn_reduced
  
  disrupt_r = disruption_status(network=network, exam_t=exam_t, disruption_times=disruption_times,
                                disrupt_t=disrupt_t, prct_node_down=prct_node_down,
                                node_down_to_capacity=node_down_to_capacity,
                                recovery_node_each_time=recovery_node_each_time, seed=seed)
  
  if (!(is.null(disrupt_r$focal_prod))) {
    focal_prod_matrix[written_row,] = disrupt_r$focal_prod
    written_row = written_row + 1
  }
  
  ######### Section 3. Searching / Navigation costs #########
  
  sim_now = sim_now + 1

}

focal_prod_matrix = cbind.data.frame(sim = 1:sim_round, focal_prod_matrix)

######### CUG test on the network

######### Plot the longitudinal trend
focal_prod_matrix_long = focal_prod_matrix %>% gather(time, prod_level, t1:t10)
focal_prod_matrix_long$time = sapply(focal_prod_matrix_long$time, FUN=substr, 2, 3)
focal_prod_matrix_long$time = as.numeric(focal_prod_matrix_long$time)

ggplot(data = focal_prod_matrix_long, aes(x = time, y = prod_level)) +
  geom_line(aes(group = sim), alpha = 0.3) +
  stat_summary(group = 1, fun.y = mean, geom = "line", lwd = 1.5) +
  stat_summary(group = 1, fun.y = mean, geom = "point", size = 3) +
  theme_bw() +
  xlab("Time") +
  ylab("Focal node product level") +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))

######### Additional. Plot the graph
sn_1 = sn$sn_ori
sn_2 = sn$sn_reduced

plot(sn_1, vertex.label = NA, vertex.size=4, edge.arrow.size=0.2, vertex.color=V(sn_1)$color)
par(oma=c(0, 0, 0, 5))
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend=paste("tier", 0:tier), col=categorical_pal(tier+1), pch=16)

deg.dist = degree_distribution(sn_1, cumulative=T)
plot(x=0:max(degree(sn_1)), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")

plot(sn_2, vertex.label = NA, vertex.size=4, edge.arrow.size=0.2, vertex.color=(1+V(sn_2)$tier))
par(oma=c(0, 0, 0, 5))
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend=paste("tier", 0:tier), col=categorical_pal(tier+1), pch=16)
