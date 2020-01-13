library(statnet)
library(ggplot2)
library(ggnetwork)
library(openxlsx)
library(reshape2)



################################### PART1: NW CONSTRUCT ########################################

# 30% low ses
# overall edge prob: p0 = 0.5
# assorative mixing: 50% edges nodematched  
  # 50% edges in adj[1:3,4:10]; the other 50% in the rest
  # prob of an edge within group: p0 * 50% / [(3*3+7*7)/(10*10)] = 0.5 * 0.5 / 0.58 = 0.431
  # prob of an edge across group: p0 * 50% / [(3*7*2)/(10*10)] = 0.5 * 0.5 / 0.42 = 0.595
# edges with weights: nodematched edges twice stronger as across-group edges 
# for this scenario: create 10 networks with the same setting 

# container for adjacency matrices
adjacency_list <- list()

# container for networks
nw_list <- list()


for(k in 1:10){
  
  # initialize adjacency matrix 
  adj <- matrix(0, nc=10, nr=10)
  
  # node matched within group: adj[1:3,1:3] & adj[4:10, 4:10]
  for(i in 1:3){
    for(j in 1:3){
      adj[i,j] <- rbinom(1,1,0.431)
    }
  }

  for(i in 4:10){
    for(j in 4:10){
      adj[i,j] <- rbinom(1,1,0.431)
    }
  }
  
  # node unmatched across group: adj[4:10,1:3]
  for(i in 4:10){
    for(j in 1:3){
      adj[i,j] <- rbinom(1,1,0.595)
    }
  }
  
  # always true within node: diagnal = 1
  diag(adj) <- 1
  
  # add edge weights: contact*0.5 across group
  adj[1:3,4:10] <- 0.5*adj[1:3,4:10]
  adj[4:10,1:3] <- 0.5*adj[4:10,1:3]
  
  # make adj symmetric using only lower triangle 
  adj[upper.tri(adj)] <- t(adj)[upper.tri(adj)]
  
  # create df to record 10 communities, 3 of lower ses
  df <- data.frame(id=1:10, ses=rep(c(0,1),c(7,3)))
  
  # attach vertex attributes to network
  attr <- c("id","ses")
  for(i in 1:2){
    set.vertex.attribute(nw, attr[i], df[,i])
  }
  
  # save adj as the first adjacency matrix for this scenario
  adjacency_list[[k]] <- adj
  
  # initialize network and save to the list
  nw <- network(adj, directed=F)
  nw_list[[k]] <- nw
  
}



################################### PART2: TRANSMISSION #######################################
# still 10 repeats for a same scenario

for(k in 1:1){
  
  # number of communities(nodes)
  num_nodes <- 10
  
  # designate the first infected
  init_I <- 1
  init_nodes <- sample(1:num_nodes, init_I)
  init_nodes
  
  # create var: index_I; "column join" to the dataset df
  index_I <- data.frame(id=init_nodes, starting_I=1)
  df <- dplyr::left_join(df, index_I, by="id", all=T)
  
  # create var: node_pop, population per community, to "column join" to the dataset df
  df$node_pop <- 1000
  
  # initialize states for the 10 communities(nodes)
  I <- rep(0, num_nodes)
  S <- rep(1000, num_nodes)
  
  # modify the first infected community/communities
  I[init_nodes] <- I[init_nodes] + init_I
  S[init_nodes] <- S[init_nodes] - init_I
  
  # define parameters
  R0 <- 16
  gamma <- 0.5
  beta <- R0*gamma
  
  # simulate for 52 weeks
  t <- 52
  dt <- 1
  
  # incidence trace for each community over each step
  incidence <- matrix(0, nr=t, nc=num_nodes)
  
  # trace overall incidence 
  overall_incidence <- rep(0,t)
  
  # trace total S and I
  total_I <- rep(0,t)
  total_S <- rep(0,t)
  total_I[1] <- sum(I)
  total_S[1] <- sum(S)
  
  # SIS iterations
  for(i in 1:t){
    
    # nodewise force of infection & recovery
    node_foi <- beta * adjacency[k] %*% I # vector length 10, force of infection for each node
    node_for <- gamma * I # fixed number, force of recovery for each node
    
    # nodewise probability of infection (CDF of the Exponential Dist.)
    node_poi <- 1-exp(-node_foi)
    node_por <- 1-exp(-node_for)
    
    # sample incidence and new recovery using probability aboce
    incidence[i,] <- rbinom(length(S), S, node_poi) # row i, all 10 cols of matrix "incidence"
    new_recovery <- rbinom(length(I), I, node_por) # vector length 10
    
    # update states
    S <- S - incidence[i] + new_recovery
    I <- I + incidence[i] - new_recovery
    
    # overall incidence across communities
    overall_incidence[i] <- sum(incidence[i,])
    
    # update prevalence: total_I
    total_I[i] <- sum(I)
    total_S[i] <- sum(S)
    
    # incidence of each community saved as a 3rd vertex attribute, and copied in the dataframe: df_plot_nw
    df_plot_nw <- set.vertex.attribute(nw, "incidence", incidence[i,])
    
    # mark the first infected community as a 4th vertex attribute, and copied in the dataframe: df_plot_nw
    df_plot_nw <- set.vertex.attribute(df_plot_nw, "first", index_I$starting_I)
    
    # plots:
    
  }
  
}

