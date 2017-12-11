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
## Last update: 12/3/2017

disruption_status_rp = function (network, same_tier=0, exam_t, disruption_times, disrupt_t=NULL, prct_node_down, node_down_to_capacity, if_risk_diff=F, risk_diff_rate=1, risk_direction=T, recovery_node_each_time, seed){

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
    all_time_affected_index = NULL
    
    # One product situation -- only make sense when same_tier=0
    if (same_tier == 0) {
      # initial production state
      prod_adj = as.matrix(get.adjacency(network))
      prod_adj[,1] = 1 / sum(prod_adj[,1]) * prod_adj[,1]
      for (i in 2:gorder(network)){
        if (sum(prod_adj[,i])){
          prod_adj[,i] = sum(prod_adj[i,]) / sum(prod_adj[,i]) * prod_adj[,i]
        }
      }
      
      for (t in 1:exam_t){
        node_health = pmin((node_health + recovery_node_each_time), 1)  # to avoid output as in "next-time"...
        if (t %in% disrupt_t){  # means there is a disruption
          # randomly choose nodes
          # focal node won't be down...
          # this is RANDOM ATTACK
          not_affected = c(1,1*(runif(gorder(network)-1) > prct_node_down))
          this_time_affected_index = which((1 - not_affected) %in% 1)   # NOT name, the INDEX
          all_time_affected_index = unique(c(this_time_affected_index, all_time_affected_index))
          health_after_disrupt = 1 - (1 - not_affected) * (1 - node_down_to_capacity)
          # we model that if a firm has double disruption, the health score will be drive to the lowest one...
          # the effect of two disruption does not stack
          # i.e. the lowest is still "node_down_to_capacity"
          # we can do the percentage -- the effects can stack
          node_health = pmin(health_after_disrupt, node_health)
        }
        # this is the risk propagation part
        if (if_risk_diff & t>1) {
          severity = (1 - node_down_to_capacity)*(risk_diff_rate^(t-1))  # node_down_to = 1 - severity
          if (risk_direction) {
            mode = "in"
          } else {
            mode = "all"
          }
          this_time_vulnerable_index = unique(unlist(adjacent_vertices(network, 
                    v=V(network)[all_time_affected_index], mode = mode), use.names = F))  # NOT name, the INDEX
          # delete focal node
          if (sum(this_time_vulnerable_index %in% 1)) {
            this_time_vulnerable_index = this_time_vulnerable_index[-which(this_time_vulnerable_index %in% 1)]
          }
          vul_not_aff_index = setdiff(this_time_vulnerable_index, all_time_affected_index) # vulnerable but not affected
          # focal node will not be affected...
          node_health[vul_not_aff_index] = 1 - severity
          all_time_affected_index = unique(c(this_time_vulnerable_index, all_time_affected_index))
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
    }
    
    return(list(disrupt_t = disrupt_t, focal_prod = focal_total_prod, net_health_equal=net_health_equal,
                net_health_weighted_d=net_health_weighted_d, net_health_weighted_bet=net_health_weighted_bet,
                net_health_weighted_eigen=net_health_weighted_eigen))
  }
}