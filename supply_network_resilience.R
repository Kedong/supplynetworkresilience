library(igraph)
library(NetSim)

######### Section 1. Setup and Simulate networks #########
source('network_construction.r')
# network_construction(tier, sup_num_tier, p_edge_tier, same_tier=0, p_edge_same_tier, direction=T, seed=1)

set.seed(1)
tier = 4
sup_num_tier = c(1,round(runif(tier, 3, 6),0))
p_edge = 0.35
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
# All nodes are with initial health score 1
# All last-tier suppliers have products 1
# If the network is healthy, the focal node should receive 1*last-tier-suppliers
# When disruption happens, the number of products can be discounted at nodes or at routes
# Let's see how much the focal can receive...
# "Dynamic" means we allow time to recover...
# We can set the number of disruption in an examination cycle

exam_t = 10 # the examination cycle
disruption_times = 3 # how many disruptions happen -- the first always happens at the beginning
if (disruption_times == 1) {disrupt_t = 1} else { disrupt_t = c(1,sort(sample(2:exam_t, (disruption_times-1))))}
# disrupt_t is the time points that disruptions happen
focal_total_prod = numeric(exam_t)

# order and production state
prod_adj = sn$sn_reduced_adj
prod_adj[,1] = 1 / sum(prod_adj[,1]) * prod_adj[,1]
for (i in 2:gorder(sn$sn_reduced)){
  if (sum(prod_adj[,i])){
    prod_adj[,i] = sum(prod_adj[i,]) / sum(prod_adj[,i]) * prod_adj[,i]
  }
}

prct_node_down = 0.3   # percentage nodes down
node_down_to_capacity = 0.7   # only a portion of products will go through...
recovery_node_each_time = 0.1

#prct_edge_down = 0.3
#edge_down_to_capacity = 0.7
#recovery_edge_each_time = 0.1

# start the time...
node_health = rep(1, gorder(sn$sn_reduced))
for (t in 1:exam_t){
  if (t %in% disrupt_t){  # means there is a disruption
    # randomly choose nodes
    # focal node won't down...
    unhealth = 1 - (1 - c(1,1*(runif(gorder(sn$sn_reduced)-1) > prct_node_down))) * (1 - node_down_to_capacity)
    # here the unhealth is totally random -- in reality, the one got hurt may not be hurt again...
    node_health = pmin(unhealth, node_health)
  }
  # check capacity backwards...
  # last-tier capacity
  prod_temp_tier = colSums(node_health[which(V(sn$sn_reduced)$tier == tier)] * prod_adj[which(V(sn$sn_reduced)$tier == tier),])
  for (i in (tier-1):1){
    prod_temp_tier = colSums(prod_temp_tier[which(V(sn$sn_reduced)$tier == i)] * node_health[which(V(sn$sn_reduced)$tier == i)] * prod_adj[which(V(sn$sn_reduced)$tier == i),] / rowSums(prod_adj[which(V(sn$sn_reduced)$tier == i),]))
  }
  focal_total_prod[t] = prod_temp_tier[1]
  node_health = sapply(node_health + recovery_node_each_time, min, 1)
}
