library(EpiModel)
library(ggplot2)
library(gridExtra)
library(grid)
library(openxlsx)
library(reshape2)


########################################## BACKGROUND INFO #########################################
### proportion of disadvantaged group ses_1 ###
# 12%
# 30%

### assortative mixing level ###
# 10% 
# 25% 
# 50% 
# 75% 
# 90% 
# 100%

### susceptible disparity: proportions not protected by vaccination ###
### 1-91.1% (the overall coverage of mealses vaccine in the U.S.)
# susceptibility 0.091 equally shared by both ses groups
# 0.091*1.25 for ses1; rest attributed to ses0
# 0.091*1.5 for ses1; rest attributed to ses0
# 0.091*1.75 for ses1; rest attributed to ses0
# 0.091*2 for ses1; rest attributed to ses0

# population size: 5000
# mean degree: 10
# network formation predicted by 2 factors: (1)edges (2)assortative mixing level by ses
####################################################################################################



######################################### CONSTANT NUMBERS #########################################
# p = number of assortative mixing scenarios: 7
p <- 7
# q = number of different levels of disparity in susceptibility: 5 
q <- 5
# r = number of different proportions of lowses: 2
r <- 2
# n = size of total population
n <- 1000
# mean_degree = mean degree for each node: 10
mean_degree <- 10
# edge_duration = how long generally a social contact lasts: 36 steps
edge_duration <- 36
# overall susceptibility: overall inf.prob=0.091 (assuming dispaired when by ses group)
susc_overall <- 0.091
# act_rate = frequency of actual contact of people in an edge by week (assuming no difference by group)
act_rate <- 1 
# rec_rate = 1/duration of disease which is 2 weeks for measles (assuming no difference by group)
rec_rate <-0.5 
# init_prev = overall baseline prevalence: 5%
init_prev <- 0.05
# m = number of simulations for each netsim: 10 
nsim <- 2
# nsteps = how long each one simulation runs: 480 steps
nstep <- 480
####################################################################################################



########################################### INDEX VECTORS ##########################################
# (index variable i) for proportion of low ses 
ind_lowses <- c(0.12,0.3)

# (index variable k) for assortative mixing level
ind_nodematch <- c(0,0.1,0.25,0.5,0.75,0.9,1)

# (index variable j) for susceptibility disparity level
# susceptibility "shared" by well-being(ses0) group, holding overall susceptibiity constant
ind_susc <- c(susc_overall,1.25*susc_overall,1.5*susc_overall,1.75*susc_overall,2*susc_overall)

# define function to show proportion as percentage
percent <- function(x, digits=2, format="f", ...){
  paste(format(100*x, format=format, digits=digits, ...),"%")
}

# create tags for p scenarios of gradient assortative mixing levels
nodematch_tag <- percent(ind_nodematch)
tag <- rep(0,p)
for(k in 1:p){
  tag[k] <- paste("scenario", k, ":",nodematch_tag[k])
}
####################################################################################################



################################## INITIALIZE VECTORS & MATRICES ###################################
# record prevalence for each simulation among m simulations of a netsim
prev_overall <- matrix(0, nr=nstep, nc=nsim+1)
prev_ses1 <- matrix(0, nr=nstep, nc=nsim+1)
prev_ses0 <- matrix(0, nr=nstep, nc=nsim+1)
PR <- matrix(0, nr=nstep, nc=nsim+1)

# record incidence for each simulation among m simulations of a netsim
IR_overall <- matrix(0, nr=nstep, nc=nsim+1)
IR_ses1 <- matrix(0, nr=nstep, nc=nsim+1)
IR_ses0 <- matrix(0, nr=nstep, nc=nsim+1)
RR <- matrix(0, nr=nstep, nc=nsim+1)

# record averaged endemic levels across simulations
prev_across_scenarios <- data.frame(scenario=1:p, 
                               overall_prev_mean=rep(0,p), overall_prev_sd=rep(0,p), 
                               ses1_prev_mean=rep(0,p), ses1_prev_sd=rep(0,p), 
                               ses0_prev_mean=rep(0,p), ses0_prev_sd=rep(0,p),
                               PR_mean=rep(0,p), PR_sd=rep(0,p))

# record avereaged incidence curves across simulations
IR_across_scenarios <- data.frame(scenario=1:p, 
                                    overall_IR_mean=rep(0,p), overall_IR_sd=rep(0,p), 
                                    ses1_IR_mean=rep(0,p), ses1_IR_sd=rep(0,p), 
                                    ses0_IR_mean=rep(0,p), ses0_IR_sd=rep(0,p),
                                    RR_mean=rep(0,p), RR_sd=rep(0,p))

# record formation setting for 7 assortative mixing scenarios
formation <- list()
target.stats <- list()
est <- list()

# record prevalence_overall prevalence_ses1 prevalence_ses0 and prevalence_ratio into 4 matrices (7*5)  
prevalence_overall <- matrix(0,ncol=q,nrow=p)
prevalence_ses1 <- matrix(0,ncol=q,nrow=p)
prevalence_ses0 <- matrix(0,ncol=q,nrow=p)
prevalence_ratio <- matrix(0,ncol=q,nrow=p)

# record risk_overall risk_ses1 risk_ses0 and risk_ratio into 4 matrices (7*5)  
risk_overall <- matrix(0,ncol=q,nrow=p)
risk_ses1 <- matrix(0,ncol=q,nrow=p)
risk_ses0 <- matrix(0,ncol=q,nrow=p)
risk_ratio <- matrix(0,ncol=q,nrow=p)

# each list below has 2 matrices for 2 proportions of low ses
sum_prevalence_overall <- list()
sum_prevalence_ses1 <- list()
sum_prevalence_ses0 <- list()
sum_prevalence_ratio <- list()
sum_risk_overall <- list()
sum_risk_ses1 <- list()
sum_risk_ses0 <- list()
sum_risk_ratio <- list()

# create a list for all plots fof following 2 proportions of low ses group
allplots <- list()
####################################################################################################



########################################## ITERATIONS ##############################################
# iterate: proportions of lowses
for(i in 1:r){
  
  # assign low ses proportion
  lowses <- ind_lowses[i]
  
  # define blank network with 2 modes of ses
  nw <- network.initialize(n, bipartite=n*lowses, directed=F)
  set.vertex.attribute(nw, "ses", rep(c(1,0),times=c(n*lowses,n-n*lowses)))
  bipvals(nw, mode=1, "ses")
  
  # create a list for all plots for follwing 4 susc disparity levels 
  suscplots <- list()

  # iterate: susceptibility disparity levels
  for(j in 1:q){
    
    # assign susceptibility score shared by ses0 group
    susc_ses1 <- ind_susc[j]
    susc_ses0 <- (susc_overall-lowses*susc_ses1)/(1-lowses)
    
    # create a list for plots of following scenarios of assortative mixing
    assortplots <- list()
    
    # set up edge dynamics
    
      # k from 1 to p-1
      k <- 1
      while (k<p){
        formation[[k]] <- ~edges + nodematch("ses", diff=F)
        target.stats[[k]] <- c(n*mean_degree/2,ind_nodematch[[k]]*n*mean_degree/2)
        k <- k+1
      }
    
      # last scenario when k=p
      formation[[k]] <- ~edges + nodemix("ses", levels2=2)
      target.stats[[k]] <- c(n*mean_degree/2,0)
    
    # iterate: assortative mixing scenarios
    for(k in 1:p){
      
      # generate dynamic contacts over blank nerwork
      coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration=edge_duration)
      est[[k]] <- netest(nw,formation[[k]],target.stats[[k]],coef.diss)
      
      # parametrize transmission over network with contacts 
      
      # inf.prob (node feature): risk of transmission given contact with an infected person
      param <- param.net(inf.prob = susc_ses1,
                        inf.prob.m2 = susc_ses0,
                        act.rate=act_rate, rec.rate=rec_rate, rec.rate.m2=rec_rate)
      
      # assign baseline prevalence in ses0 and ses1 according to susc shared by group
      init <- init.net(i.num = n*lowses*init_prev*susc_ses1/susc_overall, 
                       i.num.m2 = n*(1-lowses)*init_prev*susc_ses0/susc_overall)
      
      # simulation controls
      control <- control.net(type="SIS", nsteps=nstep, nsims=nsim, epi.by="ses", verbose.int=0)
      
      # simulate
      sim <- netsim(est[[k]], param, init, control)
      sim
      
      # plot I curves and open a null device to record
      pdf(NULL)
      dev.control(displaylist="enable")
      plot(sim, y=c("i.num", "i.num.m2"), legend=T,
           popfrac=F, sim.lines=T, qnts=F, mean.smooth=F,
           main=paste("sizes of I state -",tag[k]))
      assortplots[[5*k-4]] <- recordPlot()
      invisible(dev.off())
      
      # record overall and by-group incidence for each time step
      for(s in 1:nsim){
        df <- head(as.data.frame(get_sims(sim,sim=s)),nstep)
        prev_overall[,s] <- (df$i.num + df$i.num.m2)/n
        prev_ses1[,s] <- df$i.num/df$num
        prev_ses0[,s] <- df$i.num.m2/df$num.m2
        IR_overall[,s] <- (df$si.flow + df$si.flow.m2)/n
        IR_ses1[,s] <- df$si.flow/df$num
        IR_ses0[,s] <- df$si.flow.m2/df$num.m2
      }
      
      # extract average incidence across nsim simulations and calculate the realtive
      prev_overall[,nsim+1] <- rowMeans(prev_overall[,1:nsim])
      prev_ses1[,nsim+1] <- rowMeans(prev_ses1[,1:nsim])
      prev_ses0[,nsim+1] <- rowMeans(prev_ses0[,1:nsim])
      PR <- prev_ses1/prev_ses0
      
      # extract average prevalence across nsim simulations and calculate the realtive
      IR_overall[,nsim+1] <- rowMeans(IR_overall[,1:nsim])
      IR_ses1[,nsim+1] <- rowMeans(IR_ses1[,1:nsim])
      IR_ses0[,nsim+1] <- rowMeans(IR_ses0[,1:nsim])
      RR <- IR_ses1/IR_ses0
      
      # calculate mean and std for endemic prevalence over last 100 steps of averaged simulation
      prev_across_scenarios$overall_prev_mean[k] <- mean(prev_overall[(nstep-100):nstep,nsim+1])
      prev_across_scenarios$overall_prev_sd[k] <- sd(prev_overall[(nstep-100):nstep,nsim+1])
      prev_across_scenarios$ses1_prev_mean[k] <- mean(prev_ses1[(nstep-100):nstep,nsim+1])
      prev_across_scenarios$ses1_prev_sd[k] <- sd(prev_ses1[(nstep-100):nstep,nsim+1])
      prev_across_scenarios$ses0_prev_mean[k] <- mean(prev_ses0[(nstep-100):nstep,nsim+1])
      prev_across_scenarios$ses0_prev_sd[k] <- sd(prev_ses0[(nstep-100):nstep,nsim+1])
      prev_across_scenarios$PR_mean[k] <- mean(PR[(nstep-100):nstep,nsim+1],na.rm=T)
      prev_across_scenarios$PR_sd[k] <- sd(PR[(nstep-100):nstep,nsim+1],na.rm=T)
      
      # calculate mean and std for incidence over the course of averaged simulation
      IR_across_scenarios$overall_IR_mean[k] <- mean(IR_overall[1:nstep,nsim+1])
      IR_across_scenarios$overall_IR_sd[k] <- sd(IR_overall[1:nstep,nsim+1])
      IR_across_scenarios$ses1_IR_mean[k] <- mean(IR_ses1[1:nstep,nsim+1])
      IR_across_scenarios$ses1_IR_sd[k] <- sd(IR_ses1[1:nstep,nsim+1])
      IR_across_scenarios$ses0_IR_mean[k] <- mean(IR_ses0[1:nstep,nsim+1])
      IR_across_scenarios$ses0_IR_sd[k] <- sd(IR_ses0[1:nstep,nsim+1])
      IR_across_scenarios$RR_mean[k] <- mean(RR[1:nstep,nsim+1],na.rm=T)
      IR_across_scenarios$RR_sd[k] <- sd(RR[1:nstep,nsim+1],na.rm=T)
      
      # write each across_scenarios as a data frame
      plot_prev_df <- data.frame(t=1:nstep,
                                 prev_overall=prev_overall[,nsim+1],
                                 prev_ses1=prev_ses1[,nsim+1],
                                 prev_ses0=prev_ses0[,nsim+1], 
                                 PR=PR[,nsim+1])
      
      plot_ir_df <- data.frame(t=1:nstep,
                               IR_overall=IR_overall[,nsim+1],
                               IR_ses1=IR_ses1[,nsim+1],
                               IR_ses0=IR_ses0[,nsim+1], 
                               RR=RR[,nsim+1])
      
      # plot prevalence and open a null device to record
      assortplots[[5*k-3]] <- 
        ggplot(plot_prev_df) +
          geom_line(aes(t,prev_overall, color="overall")) +
          geom_line(aes(t,prev_ses1, color="ses1")) +
          geom_line(aes(t,prev_ses0, color="ses0")) +
          labs(colour="population") +
          ggtitle(paste("prevalence for sub-populations -",tag[k]))
      
      # plot incidence and open a null device to record
      assortplots[[5*k-2]] <- 
        ggplot(plot_ir_df) +
          geom_line(aes(t,IR_overall, color="overall")) +
          geom_line(aes(t,IR_ses1, color="ses1")) +
          geom_line(aes(t,IR_ses0, color="ses0")) +
          labs(colour="population") +
          ggtitle(paste("incidence for sub-populations -",tag[k]))
      
      # plot PR and open a null device to record
      assortplots[[5*k-1]] <- 
        ggplot(plot_prev_df) +
          geom_line(aes(t,PR)) +
          ggtitle(paste("relative prevalence ses1 to ses0 -",tag[k]))
      
      # plot RR and open a null device to record
      assortplots[[5*k]] <- 
        ggplot(plot_ir_df) +
          geom_line(aes(t,RR)) +
          ggtitle(paste("relative risk ses1 to ses0 -",tag[k]))
    }
    
    # summarize the calculated means
    # overall prevalence
    prevalence_overall[,j] <- prev_across_scenarios$overall_prev_mean
    # ses1 prevalence
    prevalence_ses1[,j] <- prev_across_scenarios$ses1_prev_mean 
    # ses0 prevalence
    prevalence_ses0[,j] <- prev_across_scenarios$ses0_prev_mean 
    # PR ses1 to ses0
    prevalence_ratio[,j] <- prev_across_scenarios$PR_mean  
    # overall risk
    risk_overall[,j] <- IR_across_scenarios$overall_IR_mean
    # ses1 risk
    risk_ses1[,j] <- IR_across_scenarios$ses1_IR_mean 
    # ses0 risk
    risk_ses0[,j] <- IR_across_scenarios$ses0_IR_mean 
    # RR ses1 to ses0
    risk_ratio[,j] <- IR_across_scenarios$RR_mean  
    
    # output table as plot: endemic prevalence across scenarios, for this j
    prev_across_scenarios[,-1] <- round(prev_across_scenarios[,-1],3)
    pdf(NULL,width=12)
    dev.control(displaylist="enable")
    grid.table(prev_across_scenarios)
    assortplots[[5*p+1]] <- recordPlot()
    invisible(dev.off())
    
    # output table as plot: incidence across scenarios, for this j
    IR_across_scenarios[,-1] <- round(IR_across_scenarios[,-1],3)
    pdf(NULL,width=12)
    dev.control(displaylist="enable")
    grid.table(IR_across_scenarios)
    assortplots[[5*p+2]] <- recordPlot()
    invisible(dev.off())

    # assortplots of 7 scenarios as the jth element of suscplots
    suscplots[[j]] <- assortplots
  }
  
  # save the summary matrices (5*7 each), for this i
  sum_prevalence_overall[[i]] <- prevalence_overall
  sum_prevalence_ses1[[i]] <- prevalence_ses1
  sum_prevalence_ses0[[i]] <- prevalence_ses0
  sum_prevalence_ratio[[i]] <- prevalence_ratio
  sum_risk_overall[[i]] <- risk_overall
  sum_risk_ses1[[i]] <- risk_ses1
  sum_risk_ses0[[i]] <- risk_ses0
  sum_risk_ratio[[i]] <- risk_ratio
  
  # plot list of 5 levels as the ith element of allplots
  allplots[[i]] <- suscplots
}
####################################################################################################



############################################## OUTPUT ############################################## 
# check completion of plotting
length(allplots) # should be r=2
length(allplots[[1]]) # should be q=5
length(allplots[[1]][[1]]) # should be 5*p+1=37

# check completion of summary
length(sum_risk_ratio) # should be r=2
dim(sum_risk_ratio[[1]]) # should be p=7 q=5 

# display 16 matrices (5*7 each, to be contoured)
sum_prevalence_overall[[1]] 
sum_prevalence_ses1[[1]] 
sum_prevalence_ses0[[1]] 
sum_prevalence_ratio[[1]]

sum_risk_overall[[1]]
sum_risk_ses1[[1]]
sum_risk_ses0[[1]]
sum_risk_ratio[[1]]

sum_prevalence_overall[[2]] 
sum_prevalence_ses1[[2]] 
sum_prevalence_ses0[[2]] 
sum_prevalence_ratio[[2]]

sum_risk_overall[[2]]
sum_risk_ses1[[2]]
sum_risk_ses0[[2]]
sum_risk_ratio[[2]]

# output above 8 matrices into an excel workbook
key <- c("0.12lowses_prev_overall", "0.12lowses_prev_ses1", "0.12lowses_prev_ses0", "0.12lowses_PR",
         "0.12lowses_risk_overall", "0.12lowses_risk_ses1", "0.12lowses_risk_ses0", "0.12lowses_RR",
         "0.3lowses_prev_overall", "0.3lowses_prev_ses1", "0.3lowses_prev_ses0", "0.3lowses_PR",
         "0.3lowses_risk_overall", "0.3lowses_risk_ses1", "0.3lowses_risk_ses0", "0.3lowses_RR")

output_matrices <- list(sum_prevalence_overall[[1]], 
                        sum_prevalence_ses1[[1]], 
                        sum_prevalence_ses0[[1]],
                        sum_prevalence_ratio[[1]],
                        sum_risk_overall[[1]],
                        sum_risk_ses1[[1]],
                        sum_risk_ses0[[1]],
                        sum_risk_ratio[[1]],
                        sum_prevalence_overall[[2]], 
                        sum_prevalence_ses1[[2]], 
                        sum_prevalence_ses0[[2]], 
                        sum_prevalence_ratio[[2]],
                        sum_risk_overall[[2]],
                        sum_risk_ses1[[2]],
                        sum_risk_ses0[[2]],
                        sum_risk_ratio[[2]])

write.xlsx(output_matrices,"/home/shiyaow/matrices.xlsx", sheetName=key)

# open 10 pdf devices to record 10 pdf files
# pdf("/Users/shiyaow/Desktop/1024.pdf")
pdf("/home/shiyaow/1028_r1_q1_7p_sim2.pdf",width=12)
allplots[[1]][[1]]
dev.off()

pdf("/home/shiyaow/1028_r1_q2_7p_sim2.pdf",width=12)
allplots[[1]][[2]]
dev.off()

pdf("/home/shiyaow/1028_r1_q3_7p_sim2.pdf",width=12)
allplots[[1]][[3]]
dev.off()

pdf("/home/shiyaow/1028_r1_q4_7p_sim2.pdf",width=12)
allplots[[1]][[4]]
dev.off()

pdf("/home/shiyaow/1028_r1_q5_7p_sim2.pdf",width=12)
allplots[[1]][[5]]
dev.off()

pdf("/home/shiyaow/1028_r2_q1_7p_sim2.pdf",width=12)
allplots[[2]][[1]]
dev.off()

pdf("/home/shiyaow/1028_r2_q2_7p_sim2.pdf",width=12)
allplots[[2]][[2]]
dev.off()

pdf("/home/shiyaow/1028_r2_q3_7p_sim2.pdf",width=12)
allplots[[2]][[3]]
dev.off()

pdf("/home/shiyaow/1028_r2_q4_7p_sim2.pdf",width=12)
allplots[[2]][[4]]
dev.off()

pdf("/home/shiyaow/1028_r2_q5_7p_sim2.pdf",width=12)
allplots[[2]][[5]]
dev.off()
