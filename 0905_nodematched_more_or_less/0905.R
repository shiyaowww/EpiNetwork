library(EpiModel)
library(ggplot2)

#####################################
# 30% low ses; 10% edges nodematched
#####################################
nw1 <- network.initialize(500, directed=F)
set.vertex.attribute(nw1,"ses",rbinom(500,1,0.3))
# formation predicted by 3: (1)edges (2)mixing by ses (3)mean degree by ses
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodematch: let 10% edges nodematched, 250*0.1 = 25
# nodefactor: ses1 mean degree * size of ses1 population = 1*500*0.3 = 150 (number of nodes of ses1 that in an edge)
formation1 <- ~edges + nodematch("ses", diff=F) + nodefactor("ses")
target.stats1 <- c(250,25,150)
coef.diss1 <- dissolution_coefs(dissolution = ~offset(edges), duration=36)

# only very limited types of dissolution models are supported in EpiModel
# homogeneous, conditional on existence of the edge
est1 <- netest(nw1, formation1, target.stats1, coef.diss1)

# Diagnostics
dx1 <- netdx(est1, nsims=5, nsteps=500, 
             nwstats.formula = ~edges + nodematch("ses") + nodefactor("ses", base=0))
dx1
plot(dx1)
par(mfrow=c(1,2))
plot(dx1, type="duration")
plot(dx1, type="dissolution")
par(mfrow=c(1,1))

# SIS over Network
# # contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param1 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assign a baseline prevalence of 1%
status.vector1 <- rbinom(500, 1, 0.01)
status.vector1 <- ifelse(status.vector1 == 1, "i", "s")
init1 <- init.net(status.vector1 = status.vector1)
# simulation control
control1 <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "ses", verbose.int=0)
# simulation
sim1 <- netsim(est1, param1, init1, control1)
sim1
plot(sim1, main = "State Prevalences 30% low ses 10% nodematched")
plot(sim1, popfrac = FALSE, main = "State Sizes 30% low ses 10% nodematched", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)



#####################################
# 30% low ses; 25% edges nodematched
#####################################
nw2 <- network.initialize(500, directed=F)
set.vertex.attribute(nw2,"ses",rbinom(500,1,0.3))
# formation predicted by 3: (1)edges (2)mixing by ses (3)mean degree by ses
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodematch: let 25% edges nodematched, 250*0.25 = 62.5
# nodefactor: ses1 mean degree * size of ses1 population = 1*500*0.3 = 150 (number of nodes of ses1 that in an edge)
formation2 <- ~edges + nodematch("ses", diff=F) + nodefactor("ses")
target.stats2 <- c(250,62.5,150)
coef.diss2 <- dissolution_coefs(dissolution = ~offset(edges), duration=36)

# only very limited types of dissolution models are supported in EpiModel
# homogeneous, conditional on existence of the edge
est2 <- netest(nw2, formation2, target.stats2, coef.diss2)

# Diagnostics
dx2 <- netdx(est2, nsims=5, nsteps=500, 
             nwstats.formula = ~edges + nodematch("ses") + nodefactor("ses", base=0))
dx2
plot(dx2)
par(mfrow=c(1,2))
plot(dx2, type="duration")
plot(dx2, type="dissolution")
par(mfrow=c(1,1))

# SIS over Network
# # contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param2 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assign a baseline prevalence of 1%
status.vector2 <- rbinom(500, 1, 0.01)
status.vector2 <- ifelse(status.vector2 == 1, "i", "s")
init2 <- init.net(status.vector2 = status.vector2)
# simulation control
control2 <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "ses", verbose.int=0)
# simulation
sim2 <- netsim(est2, param2, init2, control2)
sim2
plot(sim2, main = "State Prevalences 30% low ses 25% nodematched")
plot(sim2, popfrac = FALSE, main = "State Sizes 30% low ses 25% nodematched", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)



#####################################
# 30% low ses; 50% edges nodematched
#####################################
nw3 <- network.initialize(500, directed=F)
set.vertex.attribute(nw3,"ses",rbinom(500,1,0.3))
# formation predicted by 3: (1)edges (2)mixing by ses (3)mean degree by ses
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodematch: let 50% edges nodematched, 250*0.5 = 125
# nodefactor: ses1 mean degree * size of ses1 population = 1*500*0.3 = 150 (number of nodes of ses1 that in an edge)
formation3 <- ~edges + nodematch("ses", diff=F) + nodefactor("ses")
target.stats3 <- c(250,125,150)
coef.diss3 <- dissolution_coefs(dissolution = ~offset(edges), duration=36)

# only very limited types of dissolution models are supported in EpiModel
# homogeneous, conditional on existence of the edge
est3 <- netest(nw3, formation3, target.stats3, coef.diss3)

# Diagnostics
dx3 <- netdx(est3, nsims=5, nsteps=500, 
             nwstats.formula = ~edges + nodematch("ses") + nodefactor("ses", base=0))
dx3
plot(dx3)
par(mfrow=c(1,2))
plot(dx3, type="duration")
plot(dx3, type="dissolution")
par(mfrow=c(1,1))

# SIS over Network
# # contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param3 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assign a baseline prevalence of 1%
status.vector3 <- rbinom(500, 1, 0.01)
status.vector3 <- ifelse(status.vector3 == 1, "i", "s")
init3 <- init.net(status.vector3 = status.vector3)
# simulation control
control3 <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "ses", verbose.int=0)
# simulation
sim3 <- netsim(est3, param3, init3, control3)
sim3
plot(sim3, main = "State Prevalences 30% low ses 50% nodematched")
plot(sim3, popfrac = FALSE, main = "State Sizes 30% low ses 50% nodematched", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)



#####################################
# 30% low ses; 75% edges nodematched
#####################################
nw4 <- network.initialize(500, directed=F)
set.vertex.attribute(nw4,"ses",rbinom(500,1,0.3))
# formation predicted by 3: (1)edges (2)mixing by ses (3)mean degree by ses
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodematch: let 75% edges nodematched, 250*0.75 = 187.5
# nodefactor: ses1 mean degree * size of ses1 population = 1*500*0.3 = 150 (number of nodes of ses1 that in an edge)
formation4 <- ~edges + nodematch("ses", diff=F) + nodefactor("ses")
target.stats4 <- c(250,187.5,150)
coef.diss4 <- dissolution_coefs(dissolution = ~offset(edges), duration=36)

# only very limited types of dissolution models are supported in EpiModel
# homogeneous, conditional on existence of the edge
est4 <- netest(nw4, formation4, target.stats4, coef.diss4)

# Diagnostics
dx4 <- netdx(est4, nsims=5, nsteps=500, 
             nwstats.formula = ~edges + nodematch("ses") + nodefactor("ses", base=0))
dx4
plot(dx4)
par(mfrow=c(1,2))
plot(dx4, type="duration")
plot(dx4, type="dissolution")
par(mfrow=c(1,1))

# SIS over Network
# # contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param4 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assign a baseline prevalence of 1%
status.vector4 <- rbinom(500, 1, 0.01)
status.vector4 <- ifelse(status.vector4 == 1, "i", "s")
init4 <- init.net(status.vector4 = status.vector4)
# simulation control
control4 <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "ses", verbose.int=0)
# simulation
sim4 <- netsim(est4, param4, init4, control4)
sim4
plot(sim4, main = "State Prevalences 30% low ses 75% nodematched")
plot(sim4, popfrac = FALSE, main = "State Sizes 30% low ses 75% nodematched", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)



#####################################
# 30% low ses; 90% edges nodematched
#####################################
nw5 <- network.initialize(500, directed=F)
set.vertex.attribute(nw5,"ses",rbinom(500,1,0.3))
# formation predicted by 3: (1)edges (2)mixing by ses (3)mean degree by ses
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodematch: let 90% edges nodematched, 250*0.9 = 225
# nodefactor: ses1 mean degree * size of ses1 population = 1*500*0.3 = 150 (number of nodes of ses1 that in an edge)
formation5 <- ~edges + nodematch("ses", diff=F) + nodefactor("ses")
target.stats5 <- c(250,225,150)
coef.diss5 <- dissolution_coefs(dissolution = ~offset(edges), duration=36)

# only very limited types of dissolution models are supported in EpiModel
# homogeneous, conditional on existence of the edge
est5 <- netest(nw5, formation5, target.stats5, coef.diss5)

# Diagnostics
dx5 <- netdx(est5, nsims=5, nsteps=500, 
             nwstats.formula = ~edges + nodematch("ses") + nodefactor("ses", base=0))
dx5
plot(dx5)
par(mfrow=c(1,2))
plot(dx5, type="duration")
plot(dx5, type="dissolution")
par(mfrow=c(1,1))

# SIS over Network
# # contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param5 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assign a baseline prevalence of 1%
status.vector5 <- rbinom(500, 1, 0.01)
status.vector5 <- ifelse(status.vector5 == 1, "i", "s")
init5 <- init.net(status.vector5 = status.vector5)
# simulation control
control5 <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "ses", verbose.int=0)
# simulation
sim5 <- netsim(est5, param5, init5, control5)
sim5
plot(sim5, main = "State Prevalences 30% low ses 90% nodematched")
plot(sim5, popfrac = FALSE, main = "State Sizes 30% low ses 90% nodematched", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)



#####################################
# 12% low ses; 10% edges nodematched
#####################################
nw6 <- network.initialize(500, directed=F)
set.vertex.attribute(nw6,"ses",rbinom(500,1,0.12))
# formation predicted by 3: (1)edges (2)mixing by ses (3)mean degree by ses
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodematch: let 10% edges nodematched, 250*0.1 = 25
# nodefactor: ses1 mean degree * size of ses1 population = 1*500*0.12 = 60 (number of nodes of ses1 that in an edge)
formation6 <- ~edges + nodematch("ses", diff=F) + nodefactor("ses")
target.stats6 <- c(250,25,60)
coef.diss6 <- dissolution_coefs(dissolution = ~offset(edges), duration=36)

# only very limited types of dissolution models are supported in EpiModel
# homogeneous, conditional on existence of the edge
est6 <- netest(nw6, formation6, target.stats6, coef.diss6)

# Diagnostics
dx6 <- netdx(est6, nsims=5, nsteps=500, 
             nwstats.formula = ~edges + nodematch("ses") + nodefactor("ses", base=0))
dx6
plot(dx6)
par(mfrow=c(1,2))
plot(dx6, type="duration")
plot(dx6, type="dissolution")
par(mfrow=c(1,1))

# SIS over Network
# # contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param6 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assign a baseline prevalence of 1%
status.vector6 <- rbinom(500, 1, 0.01)
status.vector6 <- ifelse(status.vector6 == 1, "i", "s")
init6 <- init.net(status.vector6 = status.vector6)
# simulation control
control6 <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "ses", verbose.int=0)
# simulation
sim6 <- netsim(est6, param6, init6, control6)
sim6
plot(sim6, main = "State Prevalences 12% low ses 10% nodematched")
plot(sim6, popfrac = FALSE, main = "State Sizes 12% low ses 10% nodematched", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)



#####################################
# 12% low ses; 25% edges nodematched
#####################################
nw7 <- network.initialize(500, directed=F)
set.vertex.attribute(nw7,"ses",rbinom(500,1,0.12))
# formation predicted by 3: (1)edges (2)mixing by ses (3)mean degree by ses
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodematch: let 25% edges nodematched, 250*0.25 = 62.5
# nodefactor: ses1 mean degree * size of ses1 population = 1*500*0.12 = 60 (number of nodes of ses1 that in an edge)
formation7 <- ~edges + nodematch("ses", diff=F) + nodefactor("ses")
target.stats7 <- c(250,62.5,60)
coef.diss7 <- dissolution_coefs(dissolution = ~offset(edges), duration=36)

# only very limited types of dissolution models are supported in EpiModel
# homogeneous, conditional on existence of the edge
est7 <- netest(nw7, formation7, target.stats7, coef.diss7)

# Diagnostics
dx7 <- netdx(est7, nsims=5, nsteps=500, 
             nwstats.formula = ~edges + nodematch("ses") + nodefactor("ses", base=0))
dx7
plot(dx7)
par(mfrow=c(1,2))
plot(dx7, type="duration")
plot(dx7, type="dissolution")
par(mfrow=c(1,1))

# SIS over Network
# # contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param7 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assign a baseline prevalence of 1%
status.vector7 <- rbinom(500, 1, 0.01)
status.vector7 <- ifelse(status.vector7 == 1, "i", "s")
init7 <- init.net(status.vector7 = status.vector7)
# simulation control
control7 <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "ses", verbose.int=0)
# simulation
sim7 <- netsim(est7, param7, init7, control7)
sim7
plot(sim7, main = "State Prevalences 12% low ses 25% nodematched")
plot(sim7, popfrac = FALSE, main = "State Sizes 12% low ses 25% nodematched", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)



#####################################
# 12% low ses; 50% edges nodematched
#####################################
nw8 <- network.initialize(500, directed=F)
set.vertex.attribute(nw8,"ses",rbinom(500,1,0.12))
# formation predicted by 3: (1)edges (2)mixing by ses (3)mean degree by ses
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodematch: let 50% edges nodematched, 250*0.5 = 125
# nodefactor: ses1 mean degree * size of ses1 population = 1*500*0.12 = 60 (number of nodes of ses1 that in an edge)
formation8 <- ~edges + nodematch("ses", diff=F) + nodefactor("ses")
target.stats8 <- c(250,125,60)
coef.diss8 <- dissolution_coefs(dissolution = ~offset(edges), duration=36)

# only very limited types of dissolution models are supported in EpiModel
# homogeneous, conditional on existence of the edge
est8 <- netest(nw8, formation8, target.stats8, coef.diss8)

# Diagnostics
dx8 <- netdx(est8, nsims=5, nsteps=500, 
             nwstats.formula = ~edges + nodematch("ses") + nodefactor("ses", base=0))
dx8
plot(dx8)
par(mfrow=c(1,2))
plot(dx8, type="duration")
plot(dx8, type="dissolution")
par(mfrow=c(1,1))

# SIS over Network
# # contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param8 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assign a baseline prevalence of 1%
status.vector8 <- rbinom(500, 1, 0.01)
status.vector8 <- ifelse(status.vector8 == 1, "i", "s")
init8 <- init.net(status.vector8 = status.vector8)
# simulation control
control8 <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "ses", verbose.int=0)
# simulation
sim8 <- netsim(est8, param8, init8, control8)
sim8
plot(sim8, main = "State Prevalences 12% low ses 50% nodematched")
plot(sim8, popfrac = FALSE, main = "State Sizes 12% low ses 50% nodematched", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)



#####################################
# 12% low ses; 75% edges nodematched
#####################################
nw9 <- network.initialize(500, directed=F)
set.vertex.attribute(nw9,"ses",rbinom(500,1,0.12))
# formation predicted by 3: (1)edges (2)mixing by ses (3)mean degree by ses
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodematch: let 75% edges nodematched, 250*0.75 = 187.5
# nodefactor: ses1 mean degree * size of ses1 population = 1*500*0.12 = 60 (number of nodes of ses1 that in an edge)
formation9 <- ~edges + nodematch("ses", diff=F) + nodefactor("ses")
target.stats9 <- c(250,187.5,60)
coef.diss9 <- dissolution_coefs(dissolution = ~offset(edges), duration=36)

# only very limited types of dissolution models are supported in EpiModel
# homogeneous, conditional on existence of the edge
est9 <- netest(nw9, formation9, target.stats9, coef.diss9)

# Diagnostics
dx9 <- netdx(est9, nsims=5, nsteps=500, 
             nwstats.formula = ~edges + nodematch("ses") + nodefactor("ses", base=0))
dx9
plot(dx9)
par(mfrow=c(1,2))
plot(dx9, type="duration")
plot(dx9, type="dissolution")
par(mfrow=c(1,1))

# SIS over Network
# # contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param9 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assign a baseline prevalence of 1%
status.vector9 <- rbinom(500, 1, 0.01)
status.vector9 <- ifelse(status.vector9 == 1, "i", "s")
init9 <- init.net(status.vector9 = status.vector9)
# simulation control
control9 <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "ses", verbose.int=0)
# simulation
sim9 <- netsim(est9, param9, init9, control9)
sim9
plot(sim9, main = "State Prevalences 12% low ses 75% nodematched")
plot(sim9, popfrac = FALSE, main = "State Sizes 12% low ses 75% nodematched", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)



#####################################
# 12% low ses; 90% edges nodematched
#####################################
nw10 <- network.initialize(500, directed=F)
set.vertex.attribute(nw10,"ses",rbinom(500,1,0.12))
# formation predicted by 3: (1)edges (2)mixing by ses (3)mean degree by ses
# edges: let pop mean degree = 1, translated to 250 edges for 500 nodes
# nodematch: let 90% edges nodematched, 250*0.9 = 225
# nodefactor: ses1 mean degree * size of ses1 population = 1*500*0.12 = 60 (number of nodes of ses1 that in an edge)
formation10 <- ~edges + nodematch("ses", diff=F) + nodefactor("ses")
target.stats10 <- c(250,187.5,60)
coef.diss10 <- dissolution_coefs(dissolution = ~offset(edges), duration=36)

# only very limited types of dissolution models are supported in EpiModel
# homogeneous, conditional on existence of the edge
est10 <- netest(nw10, formation10, target.stats10, coef.diss10)

# Diagnostics
dx10 <- netdx(est10, nsims=5, nsteps=500, 
             nwstats.formula = ~edges + nodematch("ses") + nodefactor("ses", base=0))
dx10
plot(dx10)
par(mfrow=c(1,2))
plot(dx10, type="duration")
plot(dx10, type="dissolution")
par(mfrow=c(1,1))

# SIS over Network
# # contacts/actions half transmissible
# 10 risky actions per edge per month 
# disease duration 6 months (like TB): recovery rate = 0.16
param10 <- param.net(inf.prob=0.5, act.rate=10, rec.rate=0.16)
# assign a baseline prevalence of 1%
status.vector10 <- rbinom(500, 1, 0.01)
status.vector10 <- ifelse(status.vector10 == 1, "i", "s")
init10 <- init.net(status.vector10 = status.vector10)
# simulation control
control10 <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "ses", verbose.int=0)
# simulation
sim10 <- netsim(est10, param10, init10, control10)
sim10
plot(sim10, main = "State Prevalences 12% low ses 90% nodematched")
plot(sim10, popfrac = FALSE, main = "State Sizes 12% low ses 90% nodematched", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)