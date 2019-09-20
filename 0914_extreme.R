library(EpiModel)
library(statnet)
library(ggplot2)

###############################################################
# to see whether contact structure alone affect endemic level
###############################################################

##################################################
# extreme example 1: 250 edges all around 1 node
##################################################
nw1 <- network.initialize(500, directed=F)
mark <- rep(0,500)
mark[1] = 1
set.vertex.attribute(nw1, "mark",mark)

# formation predicted by 2: (1)edges (2)nodefactor: mark = 1 or 0
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodefactor: mark1 mean degree * size of mark1 population= 250*1 = 250 (number of nodes of mark1 that in an edge)
formation1 <- ~edges + nodefactor("mark")
target.stats1 <- c(250,250)
coef.diss1 <- dissolution_coefs(dissolution= ~offset(edges), duration=36)
est1 <- netest(nw1, formation1, target.stats1, coef.diss1)

dx1 <- netdx(est1, nsim=5, nsteps=100, 
             nwstats.formula = ~edges + nodefactor("mark", base=0), verbose=F, keep.tnetwork=T)
dx1
# plot(dx1)
nw_at10_1 <- get_network(dx1,sim=1, collapse=T, at=10)
plot(nw_at10_1, vertex.col="mark", main="network snapshop for extreme situation 1, overall meandeg=1")

# SIS over extreme 1
# contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param1 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assigne baseline prevalence of 1%
status.vector1 <- rbinom(500,1,0.01)
status.vector1 <- ifelse(status.vector1==1,"i","s")
init1 <- init.net(status.vector=status.vector1)
control1 <- control.net(type="SIS", nsims=50, nsteps=480, verbose.int=0)
sim1 <- netsim(est1, param1, init1, control1)
sim1
plot(sim1, popfrac=F, main="state sizes extreme situation 1, nsim=50", sim.lines=T, qnts=F, mean.smooth=F)



##################################################
# extreme example 2: 250 edges homogeneous
##################################################
nw2 <- network.initialize(500, directed=F)

# formation predicted by edges only
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
formation2 <- ~edges 
target.stats2 <- 250
coef.diss2 <- dissolution_coefs(dissolution= ~offset(edges), duration=36)
est2 <- netest(nw2, formation2, target.stats2, coef.diss2)

dx2 <- netdx(est2, nsim=5, nsteps=100, 
             nwstats.formula = ~edges, verbose=F, keep.tnetwork=T)
dx2
# plot(dx2)
nw_at10_2 <- get_network(dx2,sim=1, collapse=T, at=10)
plot(nw_at10_2, main="network snapshop for extreme situation 2, overall meandeg=1")

# SIS over extreme 1
# contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param2 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assigne baseline prevalence of 1%
status.vector2 <- rbinom(500,1,0.01)
status.vector2 <- ifelse(status.vector2==1,"i","s")
init2 <- init.net(status.vector=status.vector2)
control2 <- control.net(type="SIS", nsims=50, nsteps=480, verbose.int=0)
sim2 <- netsim(est2, param2, init2, control2)
sim2
plot(sim2, popfrac=F, main="state sizes extreme situation 1,nsim=50", sim.lines=T, qnts=F, mean.smooth=F)



##################################################
# extreme example 3: 250 edges all around 2 nodes
##################################################
nw3 <- network.initialize(500, directed=F)
mark3 <- rep(0,500)
mark3[1] = mark3[2] = 1
set.vertex.attribute(nw3, "mark",mark3)

# formation predicted by 2: (1)edges (2)nodefactor: mark = 1 or 0
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodefactor: mark1 mean degree * size of mark1 population= 125*2 = 250 (number of nodes of mark1 that in an edge)
formation3 <- ~edges + nodefactor("mark")
target.stats3 <- c(250,250)
coef.diss3 <- dissolution_coefs(dissolution= ~offset(edges), duration=36)
est3 <- netest(nw3, formation3, target.stats3, coef.diss3)

dx3 <- netdx(est3, nsim=5, nsteps=100, 
             nwstats.formula = ~edges + nodefactor("mark", base=0), verbose=F, keep.tnetwork=T)
dx3
# plot(dx3)
nw_at10_3 <- get_network(dx3,sim=1, collapse=T, at=10)
plot(nw_at10_3, vertex.col="mark", main="network snapshop for extreme situation 3, overall meandeg=1")

# SIS over extreme 3
# contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param3 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assigne baseline prevalence of 1%
status.vector3 <- rbinom(500,1,0.01)
status.vector3 <- ifelse(status.vector3==1,"i","s")
init3 <- init.net(status.vector=status.vector3)
control3 <- control.net(type="SIS", nsims=50, nsteps=480, verbose.int=0)
sim3 <- netsim(est3, param3, init3, control3)
sim3
plot(sim3, popfrac=F, main="state sizes extreme situation 3, nsim=50", sim.lines=T, qnts=F, mean.smooth=F)


