## This function constructs a directed supply network under preferential attachment
## Inputs: tier (focal=0), number of suppliers in each tier, starting with 1 as the focal
## Inputs: initial connection m0, if allow same tier connection (0 or 1), 
## Inputs: m is # edges at each time, size of tier -- can be a vector and non-integer
## Inputs: power is by default 1 -- the power law alpha, 1=linear
## Inputs: zero.appeal -- The ‘attractiveness’ of the vertices with no adjacent edges, default=0
## Inputs: use_in_degree is by default 0 -- if 1, preferential attachment is based on in-degree only (Price, 1976)
## if same_tier == 0, m[1] is automatically transformed to be 1
## if m is a number, then avg edge is m; if vector, then avg tier edge is m; if matrix, strict integer
## if transpass-tier connection is allowed, the probability of connection to further tiers
## (starting from tier-2-to-0) is "trans_tier"^tier-diff
## if transpass-tier is allowed, the tier-k supplier must have at least one connection to tier-k-1
## if transpass-tier is allowed, start status matters -- open triad is 1 (default, flexible) and 2 (forced)
## while closed triad is 0
## number of suppliers in tier 1 should be at least 2
## Output: igraph object
## Developed by Kedong Chen
## Last update: 12/3/2017

network_construction_ba = function (tier, sup_num_tier, m, sigma_m=0.1, power=1, zero.appeal=0, use_in_degree=0, same_tier=0, trans_tier=1, start_status=1, seed){

  require(igraph)
  set.seed(seed)
  
  if (length(m)!=1 & length(m)!=tier) {
    cat("Wrong length of m \n")
    return()
  }
  
  # adjacency matrix
  adj_m = matrix(0, sum(sup_num_tier), sum(sup_num_tier))
  tier_start_index = cumsum(sup_num_tier) # cumulative sum standing for the end index of tier

  if (same_tier==0) {  # strict tiered supply network
    adj_m[(1+tier_start_index[1]):tier_start_index[2],1] = 1  # implying m[1]=1
    # from tier 2, start the preferential attachment with edges m
    for (i in 2:tier){
      # assumption: # edges m at each tier follows normal distribution with mean m[i] sd sigma
      if (length(m)!=1) {
        size_temp = round(rnorm(n = sup_num_tier[i+1], mean = m[i], sd = sigma_m), digits=0)
      } else {
        size_temp = round(rnorm(n = sup_num_tier[i+1], mean = m, sd = sigma_m), digits=0)
      }
      for (j in 1:sup_num_tier[i+1]) {
        colsum = colSums(adj_m) # for in-degree
        rowcolsum = rowSums(adj_m) + colsum
        if (use_in_degree) {
          prob_temp = colsum[(1+tier_start_index[i-1]):tier_start_index[i]]^power + zero.appeal
        } else {
          prob_temp = rowcolsum[(1+tier_start_index[i-1]):tier_start_index[i]]^power + zero.appeal
        }
        if (sum(prob_temp)!=0) {
          prob_pref_attach = prob_temp / sum(prob_temp)
        } else {
          prob_pref_attach = rep(1/length(prob_temp), length(prob_temp))
        }
        # assumption: preferential attachment is on both incoming & outgoing degrees...
        adj_m[j+tier_start_index[i],(1+tier_start_index[i-1]):tier_start_index[i]][sample(c(1:sup_num_tier[i]),
                  size=size_temp[j], replace=F, prob=prob_pref_attach+1e-10)] = 1
      }
    }
    sn_ba = graph_from_adjacency_matrix(adj_m)
  } else {
    # tier-1 special, starting from open triad (flexible)
    adj_m[2:3,1] = 1
    if (length(m)!=1) {
      size_temp = round(rnorm(sup_num_tier[2], mean = m[1], sd = sigma_m), digits=0)
    } else {
      size_temp = round(rnorm(sup_num_tier[2], mean = m, sd = sigma_m), digits=0)
    }
    if (size_temp[1]>1 | size_temp[2]>1) {
      adj_m[3,2] = 1
    }
    if (start_status == 0) { # forced closed triad
      adj_m[3,2] = 1
    }
    if (start_status == 2) { # forced open triad
      adj_m[3,2] = 0
    }

    if (sup_num_tier[2]>2) {
      # we assume tier-k supplier MUST connect to at least one tier-k-1 supplier
      for (j in 3:sup_num_tier[2]) {
        colsum = colSums(adj_m)
        rowcolsum = rowSums(adj_m) + colsum
        if (use_in_degree) {
          prob_temp = colsum[1:j]^power + zero.appeal
        } else {
          prob_temp = rowcolsum[1:j]^power + zero.appeal
        }
        prob_pref_attach = prob_temp / sum(prob_temp)
        
        # Part I. the "must" portion -- connect to one tier-k-1
        # Part II. the "optional" portion -- may not need to connect if m[1] is 1
        if (size_temp[j]>1) {
          selected_temp = sample(c(1:j), size=size_temp[j], replace=F, prob=prob_pref_attach)
          if (1 %in% selected_temp) {
            adj_m[1+j,selected_temp] = 1
          } else {
            adj_m[1+j,1] = 1
            if (sum(prob_pref_attach[2:j])==0) {
              new_prob = rep(1/(j-1), j-1)
            } else {
              new_prob = prob_pref_attach[2:j]
            }
            adj_m[1+j,sample(c(2:j), size=size_temp[j]-1, replace=F, prob=new_prob)] = 1
          }
        } else {
          adj_m[1+j,1] = 1
        }
      }
    }
    # now begin tier-2 and so on...
    prob_trans_index = c(0,0)
    for (i in 2:tier){
      if (length(m)!=1) {
        size_temp = round(rnorm(n = sup_num_tier[i+1], mean = m[i], sd = sigma_m), digits=0)
      } else {
        size_temp = round(rnorm(n = sup_num_tier[i+1], mean = m, sd = sigma_m), digits=0)
      }
      prob_trans_index = c(i-1,prob_trans_index)
      for (j in 1:sup_num_tier[i+1]) {
        colsum = colSums(adj_m)
        rowcolsum = rowSums(adj_m) + colsum
        if (use_in_degree) {
          prob_temp = colsum[1:(tier_start_index[i]+j-1)]^power + zero.appeal
        } else {
          prob_temp = rowcolsum[1:(tier_start_index[i]+j-1)]^power + zero.appeal
        }
        # consider the reduced prob for transpass-connection
        trans_prob = rep(trans_tier^prob_trans_index, times=sup_num_tier[1:length(trans_tier^prob_trans_index)])
        prob_pref_attach = prob_temp/sum(prob_temp)*trans_prob[1:(tier_start_index[i]+j-1)]
        
        # Part I. the "must" portion -- connect to one tier-k-1
        # Part II. the "optional" portion -- may not need to connect if m[1] is 1
        if (size_temp[j]>1) {
          selected_temp = sample(c(1:(tier_start_index[i]+j-1)), size=size_temp[j], replace=F, prob=prob_pref_attach)
          if (sum(c((1+tier_start_index[i-1]):tier_start_index[i]) %in% selected_temp)>0) { # if any previous tier is selected...
            adj_m[j+tier_start_index[i],selected_temp] = 1
          } else {
            # tier-k-1 pick one
            if (sum(prob_pref_attach[c((1+tier_start_index[i-1]):tier_start_index[i])])==0) {
              new_prob = rep(1/(sup_num_tier[i]),sup_num_tier[i])
            } else {
              new_prob = prob_pref_attach[c((1+tier_start_index[i-1]):tier_start_index[i])]
            }
            tier_k_minus_1_picked = sample(c((1+tier_start_index[i-1]):tier_start_index[i]), size=1, replace=F, prob=new_prob)
            adj_m[j+tier_start_index[i], tier_k_minus_1_picked] = 1
            # all nodes pick one
            # delete that picked node first
            delete_that_node = c(1:(tier_start_index[i]+j-1))[-match(tier_k_minus_1_picked, c(1:(tier_start_index[i]+j-1)))]
            other_picked = sample(delete_that_node, size=size_temp[j]-1, replace=F, prob=prob_pref_attach[-match(tier_k_minus_1_picked, c(1:(tier_start_index[i]+j-1)))])
            adj_m[j+tier_start_index[i],other_picked] = 1
          }
        } else {
          selected_temp = sample(c(1:(tier_start_index[i]+j-1)), size=size_temp[j], replace=F, prob=prob_pref_attach)
          if (sum(c((1+tier_start_index[i-1]):tier_start_index[i]) %in% selected_temp)>0) {
            adj_m[j+tier_start_index[i], selected_temp] = 1
          } else {
            if (sum(prob_pref_attach[c((1+tier_start_index[i-1]):tier_start_index[i])])==0) {
              new_prob = rep(1/(sup_num_tier[i]),sup_num_tier[i])
            } else {
              new_prob = prob_pref_attach[c((1+tier_start_index[i-1]):tier_start_index[i])]
            }
            tier_k_minus_1_picked = sample(c((1+tier_start_index[i-1]):tier_start_index[i]), size=size_temp[j], replace=F, prob=new_prob)
            adj_m[j+tier_start_index[i],tier_k_minus_1_picked] = 1
          }
        }
      }
    }
    sn_ba = graph_from_adjacency_matrix(adj_m)
  }
  V(sn_ba)$name = paste0("v",0:(sum(sup_num_tier)-1))
  V(sn_ba)$tier = rep(c(0:tier),sup_num_tier)
  V(sn_ba)$color = 1+V(sn_ba)$tier
  return(list(sn_ba=sn_ba))
}