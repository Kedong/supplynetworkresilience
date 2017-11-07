######### Section 1. Setup and Simulate networks #########
library(igraph)
library(NetSim)

## We simulate a focal company's supply network
## For three types: random, scale-free, small-world

# 1. Random
n=100
p=0.05
# p is the probability for drawing an edge between two arbitrary vertices (G(n,p) graph)
rdm_1 = erdos.renyi.game(n=n, p, directed=T, type = "gnp")
plot(rdm_1, vertex.label = NA, vertex.size=3, edge.arrow.size=0.17, xlab = "Random graph model")

# 2. Scale-free
n=100
p=0.1
sf_1 = barabasi.game(n=n, power=p, directed=F) 
plot(sf_1, vertex.label = NA, vertex.size=3, edge.arrow.size=0.17, xlab = "Scale free model") 

# 3. small-world
# We use NetSim to generate directed small-world network
# a regular ring lattice structure with 19 nodes is defined as the starting point
n=100
p=0.8
nei=2
# p is the rewiring probability
sw_1 = watts.strogatz.game(dim=1, size=n, nei=1, p=p, loops = FALSE, multiple = FALSE, directed=T) #undirected
plot(sw_1, vertex.label = NA, vertex.size=3, edge.arrow.size=0.17, xlab = "Small world model")

# 4. Create supply network
# First, determine how many tiers, without focal company (tier-0)
set.seed(1)
tier = 5
# Second, determine how many suppliers in each tier, including focal company
sup_num_tier = c(1,round(runif(tier, 10, 20),0))  # uniform distribution
#edge_p_tier = runif(tier)
p_edge = 0.35  # used in random graph
p_same_tier_edge = 0.05
edge_p_tier = rep(p_edge, tier)
# Third, construct adjacency matrix
# Scenario 1. Same weight
adj_m = matrix(0, sum(sup_num_tier), sum(sup_num_tier))
same_tier = 0 #0: don't allow same-tier connection
strt_pos = 1
for (i in 2:(tier+1)){
  #adj_temp = matrix(ifelse(runif(sup_num_tier[i-1]*sup_num_tier[i])<=edge_p_tier[i-1],1,0), sup_num_tier[i], sup_num_tier[i-1])
  rdm_bipar_temp = sample_bipartite(n1=sup_num_tier[i-1], n2=sup_num_tier[i], type="gnp", p=edge_p_tier[i-1], directed=T, mode="in")
  adj_temp = as_adjacency_matrix(rdm_bipar_temp)[(1+sup_num_tier[i-1]):(sup_num_tier[i-1]+sup_num_tier[i]),(1:sup_num_tier[i-1])]
  adj_m[(strt_pos+1):(strt_pos+sup_num_tier[i]),((strt_pos-sup_num_tier[i-1]+1):strt_pos)] = as.matrix(adj_temp)
  if (same_tier) {
    # add same-tier links
    for (j in (strt_pos+2):(strt_pos+sup_num_tier[i])){
      for (k in (strt_pos+1):(j-1)){
        adj_m[j,k] = ifelse(runif(1)<=p_same_tier_edge,1,0)
      }
    }
  }
  strt_pos = strt_pos + sup_num_tier[i]
}
# Generate vertex names
sn_1 = graph_from_adjacency_matrix(adj_m)
V(sn_1)$name = paste0("v",0:(sum(sup_num_tier)-1))
V(sn_1)$tier = rep(c(0:tier),sup_num_tier)
#V(sn_1)$color = ifelse(V(sn_1)$name %in% "v0", "red", "blue")
plot(sn_1, vertex.label = NA, vertex.size=4, edge.arrow.size=0.2, vertex.color=(1+V(sn_1)$tier))
par(oma=c(0, 0, 0, 5))
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend=paste("tier", 0:tier), col=categorical_pal(tier+1), pch=16)

deg.dist = degree_distribution(sn_1, cumulative=T)
plot(x=0:max(degree(sn_1)), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")

# need to check random bipartite graph properties
# Q: is the multi-layer random bipartite still random?

# Q: how to remove links to make the graph scale-free or small-world?

######### Section 2. Resilience in the supply network #########
# sn_2 is the focal company's extended ego network
sn_2 = make_ego_graph(sn_1, order=tier, nodes="v0", mode="in")
plot(sn_2[[1]], vertex.label = NA, vertex.size=4, edge.arrow.size=0.2, vertex.color=(1+V(sn_2[[1]])$tier))
par(oma=c(0, 0, 0, 5))
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend=paste("tier", 0:tier), col=categorical_pal(tier+1), pch=16)



######### Section 3. Searching cost in the supply network #########