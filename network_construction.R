## This function constructs a supply network with random links, a combination of multiple bipartite graphs...
## Inputs: tier (focal=0), number of suppliers in each tier, starting with 1 as the focal
## Inputs: probability of link, if allow same tier connection (0 or 1), 
## Inputs: probability of same-tier link, if directed graph, seed
## Output: igraph object
## Developed by Kedong Chen
## Last update: 11/14/2017

network_construction = function (tier, sup_num_tier, p_edge_tier, same_tier=0, p_edge_same_tier, direction=T){

require(igraph)

# Construct adjacency matrix
# Scenario 1. Same weight
adj_m = matrix(0, sum(sup_num_tier), sum(sup_num_tier))
strt_pos = 1
if (direction) {g_mode="in"} else {g_mode="all"}
for (i in 2:(tier+1)){
  #adj_temp = matrix(ifelse(runif(sup_num_tier[i-1]*sup_num_tier[i])<=edge_p_tier[i-1],1,0), sup_num_tier[i], sup_num_tier[i-1])
  rdm_bipar_temp = sample_bipartite(n1=sup_num_tier[i-1], n2=sup_num_tier[i], type="gnp", p=p_edge_tier[i-1], directed=direction, mode=g_mode)
  adj_temp = as_adjacency_matrix(rdm_bipar_temp)[(1+sup_num_tier[i-1]):(sup_num_tier[i-1]+sup_num_tier[i]),(1:sup_num_tier[i-1])]
  adj_m[(strt_pos+1):(strt_pos+sup_num_tier[i]),((strt_pos-sup_num_tier[i-1]+1):strt_pos)] = as.matrix(adj_temp)
  if (same_tier) {
    # add same-tier links
    for (j in (strt_pos+2):(strt_pos+sup_num_tier[i])){
      for (k in (strt_pos+1):(j-1)){
        adj_m[j,k] = ifelse(runif(1)<=p_edge_same_tier,1,0)
      }
    }
  }
  strt_pos = strt_pos + sup_num_tier[i]
}
# Generate original network, even tier-1 are not connected to focal
sn_ori = graph_from_adjacency_matrix(adj_m)
V(sn_ori)$name = paste0("v",0:(sum(sup_num_tier)-1))
V(sn_ori)$tier = rep(c(0:tier),sup_num_tier)
V(sn_ori)$color = 1+V(sn_ori)$tier
# Generate reduced network, all tiers are eventually connected to focal
sn_reduced = make_ego_graph(sn_ori, order=tier, nodes="v0", mode=g_mode)[[1]]
adj_ori = as.matrix(get.adjacency(sn_ori))
adj_red = as.matrix(get.adjacency(sn_reduced))
  
return(list(sn_ori=sn_ori, sn_reduced=sn_reduced, sn_ori_adj=adj_ori, sn_reduced_adj=adj_red))
}