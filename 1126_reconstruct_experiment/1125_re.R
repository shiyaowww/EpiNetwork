library(statnet)
library(ggplot2)
library(ggnetwork)
library(gridExtra)
library(grid)
library(openxlsx)
library(reshape2)

setwd("/Users/shiyaow/Documents/【Summer 2019】/notebook/1125_reconstruct")

# 10 communities, 3 of lower ses
df <- data.frame(id=1:10, ses=rep(c(1,0),c(3,7)))

# create adjacency matrix, let diag=1
adj <- matrix(rbinom(10*10,1,0.4), nc=10, nr=10)
adj

# diagnal=1
diag(adj) <- 1      

# initialize adjacency matrix with edge weights, across-group contacts governed by assortative mixing scenarios
adj[1:3,4:10] <- adj[1:3,4:10] * 0.5
adj[4:10,1:3] <- adj[4:10,1:3] * 0.5

# make adj symmetric using only lower triangle 
adj[upper.tri(adj)] <- t(adj) [upper.tri(adj)]

# check adj
adj

# create network using adj matrix
nw <- network(adj, directed=F)

# attach vertex attributes to network
attr <- c("id","ses")
for(i in 1:2){
  set.vertex.attribute(nw, attr[i],df[,i])
}

# attach weights to edges
set.edge.value(nw, "weight", adj)

# nw check
nw
adj
summary(nw)
get.edge.value(nw, "weight")
as.sociomatrix(nw)
as.matrix(nw,matrix.type="edgelist")

# plot the initialized network
plot(nw,vertex.col="ses", label="id")

# summarize contact features using ergm 
model1 <- ergm(nw ~ edges + nodematch("ses", diff=F),verbose =T)
sim1 <- simulate(model1, burnin=1e+5, nsim=100, verbose=T)
nwsim <- sim1[[100]]
nwsim
nw
######## edge weights not included in adjacency matrix
######## valued network might work
######## try ergm ranks package


### SIS ###

# number of communities(nodes)
num_nodes <- 10

# designate first infected 
init_I <- 1
init_nodes <- sample(1:num_nodes,init_I)
init_nodes

# create var: index_I, to "column join" to the dataset df
index_I <- data.frame(id=init_nodes, starting_I=1)
df <- dplyr::left_join(df,index_I, by="id",all=T)

# create var: node_pop, population per community, to "column join" to the dataset df
df$node_pop <- 1000
df

# initialize states for the 10 communities
I <- rep(0,num_nodes)
S <- rep(1000,num_nodes)

# change a bit to account for the first infected
I[init_nodes] <- I[init_nodes] + 1
S[init_nodes] <- S[init_nodes] - 1

# record total population
total_pop = df$node_pop * num_nodes

# define parameters
beta <- 0.5
gamma <- 0.5

# simulate for 8 weeks
t <- 8
dt <- 1

# track incidence for each community over each step
incidence  <- matrix(0, nr=t, nc=num_nodes)

# record overall incidence over each step
overall_incidence <- rep(0,t)

# track total S and I, vector length = 8
total_I <- rep(0,t)
total_S <- rep(0,t)
total_I[1] <- sum(I)
total_S[1] <- sum(S)

# list of plots
nw_seq <- list()

# iterations
for(i in 1:t){
  
  # nodewise force of infection & recovery
  node_foi <- beta * adj %*% I # vector length 10, force of infection for each node(community)
  node_for <- gamma * adj # fixed number, force of recovery for each node
  
  # nodewise probability of infection (CDF of the Exponential Dist.)
  node_poi <- 1-exp(-node_foi)
  node_por <- 1-exp(-node_for)
  
  # sample incidence and new recovery using probability above
  incidence[i,] <- rbinom(length(S), S, node_poi) # row i, all 10 cols of matrix "incidence"
  new_recovery <- rbinom(length(I), I, node_por) # vector length 10
  
  # update states 
  S <- S - incidence[i,] + new_recovery
  I <- I + incidence[i,] - new_recovery
  
  # overall incidence across all communities
  overall_incidence[i] <- sum(incidence[i,]) 
  
  # update prevalence: total_I
  total_I[i] <- sum(I)
  total_S[i] <- sum(S)
  
  ### create one frame in the network png sequence ###
  
  # incidence of each community saved as a 3rd vertex attribute, and copied in the dataframe: plot_nwsim 
  plot_nw <- set.vertex.attribute(nw, "incidence", incidence[i,])
  
  # mark the first infected community as a 4th vertex attribute, and copied in the dataframe: plot_nwsim 
  plot_nw <- set.vertex.attribute(plot_nw, "first", index_I$starting_I)
  
  ######### df or nw here ???
  # plot_nw <- data.frame(node=1:10,incidence=nw %v% "incidence", first=nw %v% "first")
  
  # plot 3 layers
  ######### use edge attribute here
  p_nw <- 
    ggplot(plot_nw, aes(x=x,y=y,xend=xend,yend=yend)) +
    geom_edges(aes(color=as.factor(weight)),size=2, alpha=0.3,show.legend=F) + 
    geom_nodes(aes(color=as.factor(ses)), size=7, alpha=0.3) +
    geom_nodes(size=7, color="#0021A0", show.legend=F, data=function(x){x[x$first==1,]}) +
    geom_nodes(aes(alpha=incidence), size=7, color="#3F00EB", show.legend=F) +
    theme_blank() +
    ggtitle(paste("at",i,"weeks"))
    
  nw_seq[[i]] <- p_nw
}

# plot prevalence(total_I) over each step
total_I
df_prev <- data.frame(t=1:8, prev=total_I)
p_prev <- 
  ggplot(df_prev, aes(x=t,y=prev)) +
  geom_line() + 
  ggtitle("overall prevalence over time")

# plot overall incidence over each step
overall_incidence
df_incidence <- data.frame(t=1:8, incidence=overall_incidence)
p_incidence <- 
  ggplot(df_incidence, aes(x=t,y=incidence)) +
  geom_line() +
  ggtitle("overall incidence over time")

# create png sequence for transmission pattern

p_prev
p_incidence
p_nw
nw_seq
nw_seq[[3]]
