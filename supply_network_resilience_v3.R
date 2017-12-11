library(igraph)
library(NetSim)
library(ggplot2)
library(tidyr)
library(ergm) # to test the randomness of the graph
#library(readr)
#library(dplyr)

source('network_construction_psuedo_random.r')
source('disruption_status.r')

sim_round = 200
tier = 5
#sup_num_tier = c(1,round(runif(tier, 3, 6),0))
setting = 4
sup_num_tier = c(1, 15, 10, 5, 8, 12)
p_edge = 0.3
p_edge_tier = rep(p_edge, tier)  # prob can be different in diff tiers
same_tier=0
p_edge_same_tier = 0.05  # also can be different
direction=T
exam_t = 10   # the examination cycle, can be 1
disruption_times = 3 # how many disruptions happen -- the first always happens at the beginning
disrupt_t = c(1,3,5) # a user-specified disruption time points...
prct_node_down = 0.3   # percentage nodes down
node_down_to_capacity = 0.5   # only a portion of products will go through...
recovery_node_each_time = 0.1
#prct_edge_down = 0.3
#edge_down_to_capacity = 0.7
#recovery_edge_each_time = 0.1

focal_prod_matrix = matrix(NA, sim_round, exam_t)
colnames(focal_prod_matrix) = paste0("t", 1:exam_t)
network_health_equal = matrix(NA, sim_round, exam_t)
network_health_d = matrix(NA, sim_round, exam_t)
network_health_bet = matrix(NA, sim_round, exam_t)
network_health_eigen = matrix(NA, sim_round, exam_t)
colnames(network_health_equal) = paste0("t", 1:exam_t)
colnames(network_health_d) = paste0("t", 1:exam_t)
colnames(network_health_bet) = paste0("t", 1:exam_t)
colnames(network_health_eigen) = paste0("t", 1:exam_t)
network_property = matrix(NA, sim_round, 7)
colnames(network_property) = c("order", "size", "power_law_alpha", "centra_d", "centra_bet", "centra_eigen", 
                               "density")
# order is # nodes, size is # edge

written_row = 1
sim_now = 1

while (written_row <= sim_round){
  
  cat("this is sim", sim_now, ", writing row", written_row, "\n")
  
  ######### Section 1. Setup and Simulate networks #########
  
  seed = sim_now * 100 + 12345
  set.seed(seed)
  
  sn = network_construction_er(tier=tier, sup_num_tier=sup_num_tier, p_edge_tier=p_edge_tier, same_tier=same_tier,
                            p_edge_same_tier=p_edge_same_tier, direction=direction)
  
  ######### Section 2. Dynamic disruptions #########
  
  network = sn$sn_reduced
  
  disrupt_r = disruption_status(network=network, exam_t=exam_t, disruption_times=disruption_times,
                                disrupt_t=disrupt_t, prct_node_down=prct_node_down,
                                node_down_to_capacity=node_down_to_capacity,
                                recovery_node_each_time=recovery_node_each_time, seed=seed)
  
  if (!(is.null(disrupt_r$focal_prod))) {
    focal_prod_matrix[written_row,1:exam_t] = disrupt_r$focal_prod
    network_health_equal[written_row,1:exam_t] = disrupt_r$net_health_equal
    network_health_d[written_row,1:exam_t] = disrupt_r$net_health_weighted_d
    network_health_bet[written_row,1:exam_t] = disrupt_r$net_health_weighted_bet
    network_health_eigen[written_row,1:exam_t] = disrupt_r$net_health_weighted_eigen
    network_property[written_row,1:7] = c(gorder(network), gsize(network), 
                                        coef(fit_power_law(degree(network),implementation = "R.mle")),
                                        centr_degree(network, mode="all")$centralization, 
                                        centr_betw(network, directed=F)$centralization,
                                        centr_eigen(network, directed=F)$centralization, edge_density(network))
    written_row = written_row + 1
  }
  
  ######### Section 3. Searching / Navigation costs #########
  
  sim_now = sim_now + 1

}

# transform to long data #

focal_prod_matrix = cbind.data.frame(sim = 1:sim_round, focal_prod_matrix)
focal_prod_matrix_long = focal_prod_matrix %>% gather(time, prod_level, t1:t10)
focal_prod_matrix_long$time = sapply(focal_prod_matrix_long$time, FUN=substr, 2, 3)
focal_prod_matrix_long$time = as.numeric(focal_prod_matrix_long$time)

network_health_equal = cbind.data.frame(sim = 1:sim_round, network_health_equal)
network_health_equal_long = network_health_equal %>% gather(time, net_hlth_eql, t1:t10)
network_health_equal_long$time = sapply(network_health_equal_long$time, FUN=substr, 2, 3)
network_health_equal_long$time = as.numeric(network_health_equal_long$time)

network_health_d = cbind.data.frame(sim = 1:sim_round, network_health_d)
network_health_d_long = network_health_d %>% gather(time, net_hlth_deg, t1:t10)
network_health_d_long$time = sapply(network_health_d_long$time, FUN=substr, 2, 3)
network_health_d_long$time = as.numeric(network_health_d_long$time)

network_health_bet = cbind.data.frame(sim = 1:sim_round, network_health_bet)
network_health_bet_long = network_health_bet %>% gather(time, net_hlth_bet, t1:t10)
network_health_bet_long$time = sapply(network_health_bet_long$time, FUN=substr, 2, 3)
network_health_bet_long$time = as.numeric(network_health_bet_long$time)

network_health_eigen = cbind.data.frame(sim = 1:sim_round, network_health_eigen)
network_health_eigen_long = network_health_eigen %>% gather(time, net_hlth_eigen, t1:t10)
network_health_eigen_long$time = sapply(network_health_eigen_long$time, FUN=substr, 2, 3)
network_health_eigen_long$time = as.numeric(network_health_eigen_long$time)

network_property = cbind.data.frame(sim = 1:sim_round, network_property, setting = setting)

# merge data
sim_result = merge(merge(merge(merge(merge(focal_prod_matrix_long, network_health_equal_long, by=c("sim","time"), sort=F),
                   network_health_d_long, by=c("sim","time"), sort=F), network_health_bet_long, by=c("sim","time"), sort=F),
                   network_health_eigen_long, by=c("sim","time"), sort=F), network_property, by="sim", all.x=T)

write.csv(sim_result, "setting4.csv", row.names = F)

######### CUG test on the network
# CUG test determines if certain graph-level properties
# (e.g., clustering, average path length, centralization, homophily) result from chance.
# The process is as follows:
# 1.Take measurement f from the observed graph
# 2. Generate a random graph that controls for certain properties of the observed graph
# (e.g., size, number of edges, degree distribution, etc)
# 3. Take measurement f from the random graph
# 4. Repeat steps 2 and 3 many times (e.g., 1000) to produce a null distribution
# 5. Compare the observed measurement to the null distribution
# The null hypothesis of the CUG test is that the observed GLI (or function thereof) 
# was drawn from a distribution equivalent to that of said GLI evaluated (uniformly)
# on the space of all graphs conditional on one or more features.

######### Check if the graph is random or scale-free

sn.gof = gof(sn_2 ~ model)


d <- degree(sn_2, mode="in") # can be "all"
fit1 <- fit_power_law(d, implementation="plfit")
fit2 <- fit_power_law(d, implementation="R.mle")

fit1$alpha
stats4::coef(fit2)
fit1$logLik
stats4::logLik(fit2)

######### Plot the longitudinal trend

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

deg.dist = degree_distribution(sn_2, cumulative=T)
plot(x=0:max(degree(sn_2)), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")

plot(sn_2, vertex.label = NA, vertex.size=4, edge.arrow.size=0.2, vertex.color=(1+V(sn_2)$tier))
par(oma=c(0, 0, 0, 5))
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend=paste("tier", 0:tier), col=categorical_pal(tier+1), pch=16)


sn_bara = sample_pa(n = gorder(sn_2), power = 2, start.graph = sn_2)
gorder(sn_bara)
gsize(sn_bara)
plot(sn_bara, vertex.label = NA, vertex.size=4, edge.arrow.size=0.2, vertex.color=(1+V(sn_2)$tier))
par(oma=c(0, 0, 0, 5))
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend=paste("tier", 0:tier), col=categorical_pal(tier+1), pch=16)
