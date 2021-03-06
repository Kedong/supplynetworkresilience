## This function disrupts the supply network at the NODE level and return product status...
## This function is dynamic and reflects multi-period. To show one-time disruption, set parameters accordingly
## This function as of 11/15/2017 does NOT consider same-tier transporation...
## All nodes are with initial health score 1. Focal company orders "1" and allocation proportionally
## If the nodes are healthy, they can ship all incoming products, otherwise # products * health status
## Final outcome: how much the focal can receive...
## We allow nodes to recover from disruption...
## Inputs: iGraph network, tier, examination period (how long), # disruptions (the first always happens at the beginning)
## Inputs: Disruption time points -- random if not specified, prct_node_down, node_down_to_capacity, recovery_node_each_time
## Output: End node health, focal node prod in the examination period...
## Developed by Kedong Chen
## Last update: 11/23/2017

disruption_status = function (network, exam_t, disruption_times, disrupt_t=NULL, prct_node_down, node_down_to_capacity, recovery_node_each_time, seed){

  require(igraph)
  set.seed(seed)
  
  if (gorder(network)>1){
    tier = max(V(network)$tier)
    if (is.null(disrupt_t)){
      if (disruption_times == 1) {disrupt_t = 1} else { disrupt_t = c(1,sort(sample(2:exam_t, (disruption_times-1))))}
    }
    
    focal_total_prod = numeric(exam_t)
    net_health_equal = numeric(exam_t)
    net_health_weighted_d = numeric(exam_t)
    net_health_weighted_bet = numeric(exam_t)
    net_health_weighted_eigen = numeric(exam_t)
    node_health = rep(1, gorder(network))
    
    # initial production state
    prod_adj = as.matrix(get.adjacency(network))
    prod_adj[,1] = 1 / sum(prod_adj[,1]) * prod_adj[,1]
    for (i in 2:gorder(network)){
      if (sum(prod_adj[,i])){
        prod_adj[,i] = sum(prod_adj[i,]) / sum(prod_adj[,i]) * prod_adj[,i]
      }
    }
    
    for (t in 1:exam_t){
      node_health = pmin((node_health + recovery_node_each_time), 1)  # to avoid output as "next-time"...
      if (t %in% disrupt_t){  # means there is a disruption
        # randomly choose nodes
        # focal node won't be down...
        # this is RANDOM ATTACK
        unhealth = 1 - (1 - c(1,1*(runif(gorder(network)-1) > prct_node_down))) * (1 - node_down_to_capacity)
        # here the unhealth is totally random -- in reality, the one got hurt may not be hurt again...
        node_health = pmin(unhealth, node_health)
      }
      # check capacity backwards...
      # last-tier capacity
      prod_temp_tier = colSums(node_health[which(V(network)$tier == tier),drop = FALSE] * prod_adj[which(V(network)$tier == tier),,drop = FALSE])
      if (tier > 1) {
        for (i in (tier-1):1){
          # check if second-last tier becomes the final supplier...
          for (j in 1:length(which(V(network)$tier == i))){
            if (prod_temp_tier[which(V(network)$tier == i)][j] == 0){
              prod_temp_tier[which(V(network)$tier == i)][j] = rowSums(prod_adj[which(V(network)$tier == i),,drop = FALSE])[j]
            } else {
              prod_temp_tier[which(V(network)$tier == i)][j] = min(prod_temp_tier[which(V(network)$tier == i)][j],
                       rowSums(prod_adj[which(V(network)$tier == i),,drop = FALSE])[j])
            }
          }
          prod_temp_tier = colSums(prod_temp_tier[which(V(network)$tier == i),drop = FALSE] * 
                  node_health[which(V(network)$tier == i),drop = FALSE] * 
                  prod_adj[which(V(network)$tier == i),,drop = FALSE] /
                  rowSums(prod_adj[which(V(network)$tier == i),,drop = FALSE]))
        }
      }
      focal_total_prod[t] = prod_temp_tier[1]
      net_health_equal[t] = mean(node_health)
      net_health_weighted_d[t] = mean(node_health*igraph::degree(network, mode="all", normalized=T)) #ignore direction
      net_health_weighted_bet[t] = mean(node_health*betweenness(network, direct=F, normalized=T)) #ignore direction
      net_health_weighted_eigen[t] = mean(node_health*eigen_centrality(network, direct=F)$vector) #ignore direction
    }
    
    return(list(disrupt_t = disrupt_t, focal_prod = focal_total_prod, net_health_equal=net_health_equal,
                net_health_weighted_d=net_health_weighted_d, net_health_weighted_bet=net_health_weighted_bet,
                net_health_weighted_eigen=net_health_weighted_eigen))
  }
}