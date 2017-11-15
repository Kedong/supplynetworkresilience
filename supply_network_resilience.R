library(igraph)
library(NetSim)

######### Section 1. Setup and Simulate networks #########
source('network_construction.r')
# network_construction(tier, sup_num_tier, p_edge_tier, same_tier=0, p_edge_same_tier, direction=T, seed=1)

seed = 1
set.seed(seed)
tier = 3
sup_num_tier = c(1,round(runif(tier, 3, 6),0))
p_edge = 0.4
p_edge_tier = rep(p_edge, tier)  # prob can be different in diff tiers
same_tier=0
p_edge_same_tier = 0.05  # also can be different
direction=T

sn = network_construction(tier, sup_num_tier, p_edge_tier, same_tier=0, p_edge_same_tier, direction=T)

## plot the graph
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

######### Section 2. Dynamic disruptions #########
source('disruption_status.r')
network = sn$sn_reduced
exam_t = 1 # the examination cycle
disruption_times = 1 # how many disruptions happen -- the first always happens at the beginning
#disrupt_t = # a user-specified disruption time points...
prct_node_down = 0.3   # percentage nodes down
node_down_to_capacity = 0.7   # only a portion of products will go through...
recovery_node_each_time = 0.1
#prct_edge_down = 0.3
#edge_down_to_capacity = 0.7
#recovery_edge_each_time = 0.1

disrupt_r = disruption_status(network, tier, exam_t, disruption_times, disrupt_t=NULL, prct_node_down, node_down_to_capacity, recovery_node_each_time, seed)
disrupt_r$focal_prod



