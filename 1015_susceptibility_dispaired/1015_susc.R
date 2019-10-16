library(EpiModel)
library(ggplot2)
library(gridExtra)
library(grid)

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

# susceptibility score
ses0_susc_score <- c(0.2,0.15,0.1,0.05)

# label all 14 scenarios
tag <- rep(0,14)
for(i in 1:14){
  tag[i] <- paste("scenario", i, ":",lowses_tag[i],"low ses;",assortlevel_tag[i],"nodematched")
}

# c = number of assortative mixing scenarios (14 preferred)
c <- 14
# s = number of different levels of disparity in susceptibility (4 preferred)
s <- 4
# k = number of simulations for each scenario (10 preferred)
k <- 10

# record incidence for plotting
IR_overall <- matrix(0, nr=480, nc=k+1)
IR_ses1 <- matrix(0, nr=480, nc=k+1)
IR_ses0 <- matrix(0, nr=480, nc=k+1)
RR <- matrix(0, nr=480, nc=k+1)

# record I size at the end of time series, for endemic comparison across scenarios
across_scenarios <- data.frame(scenario=1:c, 
                               overall_mean=rep(0,c), overall_sd=rep(0,c), 
                               ses1_mean=rep(0,c), ses1_sd=rep(0,c), 
                               ses0_mean=rep(0,c), ses0_sd=rep(0,c))

# a list for all plots
allplots <- list()

# iterations of differential susceptibility distributions
for(m in 1:s){
  
  # create list of plots for this susceptibility distribution
  plotlist <- list()
  
  # iterations of differential assortative mixing levels
  for(i in 1:c){
 
    # blank network; define 2 modes
    nw <- network.initialize(1000, bipartite=1000*lowses[i], directed=F)
    set.vertex.attribute(nw, "ses", rbinom(1000, 1, lowses[i]))
    bipvals(nw, mode=1, "ses")
  
    # dynamic linkage construction
    formation <- ~edges + nodematch("ses", diff=F)
    target.stats <- c(1000, assortlevel[i]*1000)
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration=36)
    est <- netest(nw, formation, target.stats, coef.diss)
  
    # diagnostics of the estimated dynamic network
    # dx <- netdx(est, nsims=5, nsteps=480,
                # nwstats.formula = ~edges + nodematch("ses"))
    # plot(dx)
    # par(mfrow=c(1,2))
    # plot(dx, type="duration")
    # plot(dx, type="dissolution")
    # par(mfrow=c(1,1))
  
    # parametrize the transmission
    # inf.prob (node feature): risk of transmission given contact with an infected person
    param <- param.net(inf.prob = (0.2-0.7*ses0_susc_score[m])/0.3,
                       inf.prob.m2 = ses0_susc_score[m],
                       act.rate=2, rec.rate=0.16, rec.rate.m2=0.16)
  
    # assign a baseline prevalence of 2%, 20 people infected to begin with
    status.vector <- rbinom(1000, 1, 0.02)
    status.vector <- ifelse(status.vector == 1, "i", "s")
    init <- init.net(status.vector = status.vector)
  
    # simulation control: run for 40 yrs by month, simulate k times (10 times prefered)
    control <- control.net(type = "SIS", nsteps=480, nsims=k, epi.by = "ses", verbose.int=0)
  
    # simulate
    sim <- netsim(est,param,init,control)
    sim
  
    # plot sizes of I state; saving as an object in the plotlist (a bit tricky)
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(sim, y=c("i.num","i.num.ses0", "i.num.ses1"), legend=T,
         popfrac=F, sim.lines=T, qnts=F, mean.smooth=F,
         main=paste("sizes of i state -",tag[i]))
    plotlist[[3*i-2]] <- recordPlot()
    invisible(dev.off())
  
    # calculate overall and by-group incidence over time
    for(j in 1:k){
      df <- head(as.data.frame(get_sims(sim, sim=j)),480)
      IR_overall[,j] <- df$i.num/df$num
      IR_ses1[,j] <- df$i.num.ses0/df$num.ses0
      IR_ses0[,j] <- df$i.num.ses1/df$num.ses1
    }
  
    # for overall and by-group: extract average incidence across k simulations; calculate incidence ratio 
    IR_overall[,k+1] <- rowMeans(IR_overall[,1:k])
    IR_ses1[,k+1] <- rowMeans(IR_ses1[,1:k])
    IR_ses0[,k+1] <- rowMeans(IR_ses0[,1:k])
    RR <- IR_ses1/IR_ses0
  
    # for this scenario, using average incidence across k simulations
    # calculate mean and std for I size for last 280 steps, in order to compare across scenarios
    across_scenarios$overall_mean[i] <- mean(IR_overall[201:480,k+1])
    across_scenarios$overall_sd[i] <- sd(IR_overall[201:480,k+1])
    across_scenarios$ses1_mean[i] <- mean(IR_ses1[201:480,k+1])
    across_scenarios$ses1_sd[i] <- sd(IR_ses1[201:480,k+1])
    across_scenarios$ses0_mean[i] <- mean(IR_ses0[201:480,k+1])
    across_scenarios$ses0_sd[i] <- sd(IR_ses0[201:480,k+1])
  
    # as a data frame
    plot_ir_df <- data.frame(t=1:480,IR_overall=IR_overall[,k+1],IR_ses1=IR_ses1[,k+1],IR_ses0=IR_ses0[,k+1], RR=RR[,k+1])
    # plot_ir_df
  
    # plot incidence 
    plotlist[[3*i-1]] <-
      ggplot(plot_ir_df) +
        geom_line(aes(t,IR_overall, color="overall")) +
        geom_line(aes(t,IR_ses1, color="ses1")) +
        geom_line(aes(t,IR_ses0, color="ses0")) +
        labs(colour="population") +
        ggtitle(paste("incidence for sub-populations -",tag[i]))
   
    # plot relative risk
    plotlist[[3*i]] <-
      ggplot(plot_ir_df) +
        geom_line(aes(t,RR)) +
        ylim(0,1.5) +
        ggtitle(paste("relative risk ses1 to ses0 -",tag[i]))
  
  }

  # plot table: endemic level across scenarios, for this susceptibility distribution
  across_scenarios[,-1] <- round(across_scenarios[,-1],3)
  pdf(NULL)
  dev.control(displaylist="enable")
  grid.table(across_scenarios)
  plotlist[[3*c+1]] <- recordPlot()
  invisible(dev.off())
  
  # save this list of plots as an element of the allplots list 
  allplots[[m]] <- plotlist

}

# check plotlists in allplots
length(allplots)
allplots

# open pdf device to output plots for susceptibility distribution 1
# pdf("/Users/shiyaow/Desktop/c2_s1st.pdf")
pdf("/home/shiyaow/c2_s1st.pdf")
allplots[[1]]
dev.off()

# open pdf device to output plots for susceptibility distribution 2
# pdf("/Users/shiyaow/Desktop/c2_s2nd.pdf")
pdf("/home/shiyaow/c2_s2nd.pdf")
allplots[[2]]
dev.off()

# open pdf device to output plots for susceptibility distribution 3
# pdf("/Users/shiyaow/Desktop/c2_s3rd.pdf")
pdf("/home/shiyaow/c2_s3rd.pdf")
allplots[[3]]
dev.off()

# open pdf device to output plots for susceptibility distribution 4
# pdf("/Users/shiyaow/Desktop/c2_s4th.pdf")
pdf("/home/shiyaow/c2_s4th.pdf")
allplots[[4]]
dev.off()



