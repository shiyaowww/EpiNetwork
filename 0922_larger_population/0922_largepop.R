library(EpiModel)
library(ggplot2)

##################################################
# 30% disadvantaged group (low ses in this case):
  # homogeneous mixing
  # 10% nodematched
  # 25% nodematched
  # 50% nodematched
  # 75% nodematched
  # 90% nodematched
  # totally segregated
# 12% disadvantaged group
  # same 7 scenarios
##################################################

# compared to runs before:
  # larger population: 1000 nodes
  # more overall contacts: mean degree = 2
  # lower beta value for a single infective over a single time step
  # SIS: personal contacts = 2, beta = 0.4, gamma = 0.16, personal reproduction number = 2*0.4/0.16 = 5

# formation predicted by 2 factors: (1)edges (2)nodematched level by ses
  # edges: let pop mean degree = 2, translated to 1000 edges for 1000 nodes
  # nodematch: let 0% edges nodematched, target stat = 1000 * nodematched level = counts of nodematched edges

# show proportion as percentage
percent <- function(x, digits=2, format="f", ...){
  paste(format(100*x, format=format, digits=digits, ...),"%")
}

# proportion of disadvantaged group
lowses <- rep(0.3,14)
lowses[8:14] <- rep(0.12,7) 
lowses_tag <- percent(lowses)

# assortative mixing (nodematched) level
assortlevel <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
assortlevel <- c(assortlevel,assortlevel)
assortlevel_tag <- percent(assortlevel)

# label all 14 scenarios
tag <- rep(0,14)
for(i in 1:14){
  tag[i] <- paste("scenario", i, ":",lowses_tag[i],"low ses;",assortlevel_tag[i],"nodematched")
}

# c = counts of scenarios, 14 in this case
# k = number of simulations for each scenario, 
c <- 14
k <- 10

# record incidence for plotting
IR_overall <- matrix(0, nr=480, nc=k+1)
IR_ses1 <- matrix(0, nr=480, nc=k+1)
IR_ses0 <- matrix(0, nr=480, nc=k+1)
RR <- matrix(0, nr=480, nc=k+1)

# iteration
for(i in 1:c){
  
  # network construction
  nw <- network.initialize(1000, directed=F)
  set.vertex.attribute(nw, "ses", rbinom(1000, 1, lowses[i]))
  formation <- ~edges + nodematch("ses", diff=F)
  target.stats <- c(1000, assortlevel[i]*1000)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration=36)
  est <- netest(nw, formation, target.stats, coef.diss)
  
  # diagnostics of the estimated dynamic network
  dx <- netdx(est, nsims=5, nsteps=480,
              nwstats.formula = ~edges + nodematch("ses"))
  plot(dx)
  par(mfrow=c(1,2))
  plot(dx, type="duration")
  plot(dx, type="dissolution")
  par(mfrow=c(1,1))
  
  # parametrize the transmission 
  param <- param.net(inf.prob=0.2, act.rate=2, rec.rate=0.16)
  
  # assign a baseline prevalence of 0.2%, 2 people infected to begin with
  status.vector <- rbinom(1000, 1, 0.002)
  status.vector <- ifelse(status.vector == 1, "i", "s")
  init <- init.net(status.vector = status.vector)
  
  # simulation control: run for 40 yrs by month, simulate 10 times
  control <- control.net(type = "SIS", nsteps=480, nsims=k, epi.by = "ses", verbose.int=0)
  
  # simulate
  sim <- netsim(est,param,init,control)
  sim
  
  # plot sizes of I state
  plot(sim, y=c("i.num","i.num.ses0", "i.num.ses1"), legend=T,
       popfrac=F, sim.lines=T, qnts=F, mean.smooth=F,
       main=paste("sizes of i state -",tag[i]))
  
  # calculate incidence
  for(j in 1:k){
    df <- head(as.data.frame(get_sims(sim, sim=j)),480)
    IR_overall[,j] <- df$i.num/df$num
    IR_ses1[,j] <- df$i.num.ses0/df$num.ses0
    IR_ses0[,j] <- df$i.num.ses1/df$num.ses1
  }
  
  # extract incidence and calculate the relative 
  IR_overall[,k+1] <- rowMeans(IR_overall[,1:k])
  IR_ses1[,k+1] <- rowMeans(IR_ses1[,1:k])
  IR_ses0[,k+1] <- rowMeans(IR_ses0[,1:k])
  RR <- IR_ses1/IR_ses0
  
  # as a data frame
  plot_ir_df <- data.frame(t=1:480,IR_overall=IR_overall[,k+1],IR_ses1=IR_ses1[,k+1],IR_ses0=IR_ses0[,k+1], RR=RR[,k+1])
  # plot_ir_df
  
  # plot incidence 
  ggplot(plot_ir_df) +
    geom_line(aes(t,IR_overall, color="overall")) +
    geom_line(aes(t,IR_ses1, color="ses1")) +
    geom_line(aes(t,IR_ses0, color="ses0")) +
    labs(colour="population") +
    ylim(0,1) +
    ggtitle(paste("incidence for sub-populations -",tag[i]))
 
   # plot relative risk
  ggplot(plot_ir_df) +
    geom_line(aes(t,RR)) +
    ggtitle(paste("relative risk ses1 to ses0 -",tag[i]))
}

