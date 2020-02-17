library(statnet)
library(ggplot2)
library(ggnetwork)
library(openxlsx)
library(reshape2)
library(GGally)

#################################### PART1: CONST/VAR/CONTAINER/FUNC #######################################
# R0: reproductive number
R0 <- 16

# gamma:
gamma <- 0.5

# beta:
beta <- R0*gamma

# numRepeat: num ~10, repeat 10 times for each single scenario
numRepeat <- 1000

# timeRun: num ~52 if by week, how long we want the simulation run
timeRun <- 52

# stepLength: num ~1, dt
stepLength <- 1

# numNode: num ~10
numNode <- 10

# numPerNode: num ~1000, how many individuals per node
numPerNode <- 1000

# probLowSes: num or umeric vector ranged 0 to 1, ses0/all, increment by 1/numNode
#    part of the primary research question, thus as a vector here
probLowSes <- c(0.1, 0.3)

# suscScale: num or numeric vector, element > 1, relative susceptibility of low ses group 
#    part of the primary research question, thus as a vector here
suscScale <- seq(1, 1.5, by=0.05)

# suscVectortor: numeric vector length numNode, score of susceptibility
#    element = 1 if ses1, suscScale if ses0
suscVector <- rep(0, numNode)

# overallProbEdge: num or numeric vector ranged 0 to 1, edge prob
#    eg. overallProbEdge = 0.5 for ~20 edges in a 10-node network
#    could derive auxiliary research question, in which case it becomes a vector
overallProbEdge <- 0.5

# nodeMatch: num or numeric vector ranged 0 to 1, assortative mixing levels 
#    part of the primary research question, thus as a vector here
nodeMatch <- seq(0, 1, by=0.1)

# probEdgeWithinGroup: num ranged 0 to 1, later calculated with probLowSes
#    = overallProbEdge * nodeMatch +
#         [overallProbEdge * (1-nodeMatch)] *
#            { [probLowSes^2 + (1-probLowSes)^2] - [overallProbEdge * nodeMatch] } / 
#                { [2 * probLowSes * (1-probLowSes)] + [probLowSes^2 + (1-probLowSes)^2] - [overallProbEdge * nodeMatch] } 
probEdgeWithinGroup <- 0

# probEdgeAcrossGroup: num ranged 0 to 1, later calculated with probLowSes
#    =    [overallProbEdge * (1-nodeMatch)] *
#            [2 * probLowSes * (1-probLowSes)] / 
#                { [2 * probLowSes * (1-probLowSes)] + [probLowSes^2 + (1-probLowSes)^2] - [overallProbEdge * nodeMatch] } 
probEdgeAcrossGroup <- 0

# edgeWeightScale: num or numeric vector ranged 0 to 1, weaker edge across groups than within
#    could derive auxiliary research question, in which case it becomes a vector
edgeWeightScale <- 0.6

# probInterPerson: num or numeric vector ranged 0 to 1
#    probability of group-level "contactability"(edge weights) turned into a person-to-person contact
#    could derive auxiliary research question, in which case it becomes a vector
probInterPerson <- 0.015

# adjMatrix: adjacency matrix nrow = numNode ncol = numNode, rewritten each step of each scenario
adjMatrix <- matrix(0, nr=numNode, nc=numNode)

# nwBase: network object, rewritten each repeat of each scenario
nwBase <- network(numNode, directed=F)

# vAttr: text vector, length 2 
vAttr <- c("ID","SES")

# dfNwInfo: dataframe length numNode, width 2 at first and 4 at last 
#    rewritten each repeat of each scenario, store ID SES initInfectStatus and numPerNode

# numInitInfect: num ~1, how many individuals we want to infect in the first place
numInitInfect <- 0

# initNodeId: num or numeric vector, ranged 1 to numNode
#    sampled node id of the one(s) first infected
#    generated in loops, size may vary with diff value of numInitInfect

# dfIndexInitNode: 2-col numInitInfect-row index dataframe
#    1st col for node id: present only if initNodeId
#    2nd col named initInfectStatus: 1 for all row(s) in this dataframe
#    generated in loops, size may vary with diff value of numInitInfect

# nodeForceInfect = beta * adjMatrix %*% suscVector * I / numPerNode
#    [matrix %*% vector] for matrix multiplication, yield vector
#    [vector * vector] for hadamard(point-to-point) multiplication, yield vector
#    numeric vector length numNode, force of infection for each node 
nodeForceInfect <- rep(0,numNode)

# nodeForceRecov = gamma 
#    numeric vector length numNode, force of recovery, each element the same 
nodeForceRecov <- rep(0, numNode)

# nodeProbInfect = 1 - exp(- nodeForceInfect * stepLength), probability of infection
#    numeric vector length numNode, function of dt (stepLength)
nodeProbInfect <- rep(0, numNode)

# nodeProbRecov = 1 - exp(- nodeForceRecov * stepLength), probability of recovery
#    numeric vector length numNode, function of dt (stepLength)
nodeProbRecov <- rep(0, numNode)

# I: numeric vector, length numNode, rewritten each step of each repeat of each scenario
I <- rep(0, numNode)

# S: numeric vector, length numNode, rewritten each step of each repeat of each scenario
S <- rep(0, numNode)

# incidenceByNodeMatrix: numeric matrix, nrow = timeRun, ncol = numNode
#    rewritten each repeat of each scenario
incidenceByNodeMatrix <- matrix(0, nr=timeRun/stepLength, nc=numNode)

# newRecov: num, rewritten each step each repeat each scenario
#    not necessary to keep record of thus no storage for recovery each step
newRecov <- 0

# incidenceMatrix: numeric matrix, nrow = timeRun, ncol = numRepeat+1, rewritten each scenario
#    each incidenceByNodeMatrix collapse to 1 col of this incidenceMatrix
incidenceMatrix <- matrix(0, nr=timeRun/stepLength, nc=numRepeat+1)

# incidenceSes0Matrix: numeric matrix, nrow = timeRun, ncol = numRepeat+1, rewritten each scenario
#    each incidenceByNodeMatrix's ses0 columns collapse to 1 col of this incidenceSes0Matrix
incidenceSes0Matrix <- matrix(0, nr=timeRun/stepLength, nc=numRepeat+1)

# incidenceSes1Matrix: numeric matrix, nrow = timeRun, ncol = numRepeat+1, rewritten each scenario
#    each incidenceByNodeMatrix's ses1 columns collapse to 1 col of this incidenceSes1Matrix
incidenceSes1Matrix <- matrix(0, nr=timeRun/stepLength, nc=numRepeat+1)

# incidenceMatrices: list length the total counts of scenarios, store all incidenceMatrix
incidenceMatrices <- list()

# incidenceSes0Matrices: list length the total counts of scenarios, store all incidenceSes0Matrix
incidenceSes0Matrices <- list()

# incidenceSes1Matrices: list length the total counts of scenarios, store all incidenceSes1Matrix
incidenceSes1Matrices <- list()

# prevByNodeMatrix: numeric matrix, nrow = timeRun, ncol = numNode
#    rewritten each repeat of each scenario
#    recording not the incidence but the I state
prevByNodeMatrix <- matrix(0, nr=timeRun/stepLength, nc=numNode)

# prevMatrix: numeric matrix, nrow = timeRun, ncol = numRepeat+1, rewritten each scenario
#    each prevByNodeMatrix collapse to 1 col of this prevMatrix
prevMatrix <- matrix(0, nr=timeRun/stepLength, nc=numRepeat+1)

# prevSes0Matrix: numeric matrix, nrow = timeRun, ncol = numRepeat+1, rewritten each scenario
#    each prevByNodeMatrix's ses0 columns collapse to 1 col of this prevSes0Matrix
prevSes0Matrix <- matrix(0, nr=timeRun/stepLength, nc=numRepeat+1)

# prevSes1Matrix: numeric matrix, nrow = timeRun, ncol = numRepeat+1, rewritten each scenario
#    each prevByNodeMatrix's ses1 columns collapse to 1 col of this prevSes1Matrix
prevSes1Matrix <- matrix(0, nr=timeRun/stepLength, nc=numRepeat+1)

# prevMatrices: list length the total counts of scenarios, store all prevMatrix
prevMatrices <- list()

# prevSes0Matrices: list length the total counts of scenarios, store all prevSes0Matrix
prevSes0Matrices <- list()

# prevSes1Matrices: list length the total counts of scenarios, store all prevSes1Matrix
prevSes1Matrices <- list()

# The primary research question explores:
#    length(probLowSes) * length(nodeMatch) * length(suscScale) scenarios.
#    Each scenario has 6 columns of summary (from the last col of the 6 matrices).
#    dataframe length = timeRun/stepLength, ncol = 6 at first, 8 at last with incidenceRatio and prevRatio
#    dfSumScenario rewritten every scenario
dfSumScenario <- data.frame(incidence = rep(0,timeRun/stepLength),
                            incidenceSes0 = rep(0,timeRun/stepLength),
                            incidenceSes1 = rep(0,timeRun/stepLength),
                            prev = rep(0,timeRun/stepLength),
                            prevSes0 = rep(0,timeRun/stepLength),
                            prevSes1 = rep(0,timeRun/stepLength),
                            incidenceRatio = rep(0,timeRun/stepLength),
                            prevRatio = rep(0,timeRun/stepLength))

# avgEndemicPrev: num, averaged endemic prevalence from last steps, rewritten each scenario
#    eg. avgEndemicPrev = dfSumScenario$prev[(timeRun/stepLength-20):timeRun], averaged last 20 steps
avgEndemicPrev <- 0

# avgIncidenceRatio: num, averaged dfSumScenario$incidenceRatio over the timeRun/stepLength steps
#    rewritten each scenario
avgIncidenceRatio <- 0

# avgPrevRatio: num, averaged dfSumScenario$prevRatio over the timeRun/stepLength steps
#    rewritten each scenario
avgPrevRatio <- 0

# sdEndemicPrev: num, sd of endemic prevalence from last steps, rewritten each scenario
#    eg. avgEndemicPrev = dfSumScenario$prev[(timeRun/stepLength-20):timeRun], averaged last 20 steps
sdEndemicPrev <- 0

# sdIncidenceRatio: num, sd of dfSumScenario$incidenceRatio over the timeRun/stepLength steps
#    rewritten each scenario
sdIncidenceRatio <- 0

# sdPrevRatio: num, sd of dfSumScenario$prevRatio over the timeRun/stepLength steps
#    rewritten each scenario
sdPrevRatio <- 0

# tableAvgEP: dataframe nrow = length(nodeMatch), ncol = length(suscScale)
#    store the calculated avgEndemicPrev value of multiple scenarios
#    rewritten each iteration of probLowSes
tableAvgEP <- data.frame(matrix(NA, nrow=length(nodeMatch), ncol=length(suscScale)))
rownames(tableAvgEP) <- nodeMatch
colnames(tableAvgEP) <- suscScale

# tableAvgIR: dataframe nrow = length(nodeMatch), ncol = length(suscScale)
#    store the calculated avgIncidenceRatio value of multiple scenarios
#    rewritten each iteration of probLowSes
tableAvgIR <- data.frame(matrix(NA, nrow=length(nodeMatch), ncol=length(suscScale)))
rownames(tableAvgIR) <- nodeMatch
colnames(tableAvgIR) <- suscScale

# tableAvgPR: dataframe nrow = length(nodeMatch), ncol = length(suscScale)
#    store the calculated avgPrevRatio value of multiple scenarios
#    rewritten each iteration of probLowSes
tableAvgPR <- data.frame(matrix(NA, nrow=length(nodeMatch), ncol=length(suscScale)))
rownames(tableAvgPR) <- nodeMatch
colnames(tableAvgPR) <- suscScale

# tableSdEP: dataframe nrow = length(nodeMatch), ncol = length(suscScale)
#    store the calculated sdEndemicPrev value of multiple scenarios
#    rewritten each iteration of probLowSes
tableSdEP <- data.frame(matrix(NA, nrow=length(nodeMatch), ncol=length(suscScale)))
rownames(tableSdEP) <- nodeMatch
colnames(tableSdEP) <- suscScale

# tableSdIR: dataframe nrow = length(nodeMatch), ncol = length(suscScale)
#    store the calculated sdIncidenceRatio value of multiple scenarios
#    rewritten each iteration of probLowSes
tableSdIR <- data.frame(matrix(NA, nrow=length(nodeMatch), ncol=length(suscScale)))
rownames(tableSdIR) <- nodeMatch
colnames(tableSdIR) <- suscScale

# tableSdPR: dataframe nrow = length(nodeMatch), ncol = length(suscScale)
#    store the calculated sdPrevRatio value of multiple scenarios
#    rewritten each iteration of probLowSes
tableSdPR <- data.frame(matrix(NA, nrow=length(nodeMatch), ncol=length(suscScale)))
rownames(tableSdPR) <- nodeMatch
colnames(tableSdPR) <- suscScale

# avgTables: list length = 3 * length(probLowSes), store all "3 avg tables" for contour plots
avgTables <- list()

# sdTables: list length = 3 * length(probLowSes), store all "3 sd tables" for trouble shooting
sdTables <- list()

# avgTablesMelted：all avgTables melted for contour plots
avgTablesMelted <- list()

# scenarioKey: vector length 3 * length(probLowSes) 
#    key for output xlsx, avgTables and sdTables
scenarioKey <- rep(0, 3*length(probLowSes))

# i: loop var, for probLowSes
i <- 0

# j: loop var, for suscScale
j <- 0

# k: loop var, for nodeMatch
k <- 0

# l: loop var, for repeats 
l <- 0

# p: loop var, when building up adjacency matrix 
p <- 0

# q: loop var, when building up adjacency matrix 
q <- 0

# m: loop var, when assigning vertex attributes
m <- 0

# t: loop var, over the course of timeRun, by stepLength
t <- 0

# c: loop var, when creating contour plots
c <- 0



##################################### PART2: ITERATIONS ########################################
# iterate over diff proportions of low ses
for (i in 1:length(probLowSes))
{
  
  # iterate over diff levels of nodematching(assortative mixing)
  for (j in 1:length(nodeMatch))
  {
    # calculate probEdgeWithinGroup
    probEdgeWithinGroup <-
      overallProbEdge * nodeMatch[j] +
      overallProbEdge * (1-nodeMatch[j]) *
      ((probLowSes[i])^2 + (1-probLowSes[i]^2) - overallProbEdge * nodeMatch[j]) /
      (2 * probLowSes[i] * (1-probLowSes[i]) + 
         (probLowSes[i])^2 + (1-probLowSes[i])^2 - overallProbEdge * nodeMatch[j])
    
    # calculate probEdgeAcrossGroup
    probEdgeAcrossGroup <-
      overallProbEdge * (1-nodeMatch[j]) *
      2 * probLowSes[i] * (1-probLowSes[i]) /
      (2 * probLowSes[i] * (1-probLowSes[i]) + 
         (probLowSes[i])^2 + (1-probLowSes[i])^2 - overallProbEdge * nodeMatch[j])
    
    # iterate over diff relative susceptibility low vs. high ses
    for (k in 1:length(suscScale))
    {
      # write up suscVector for this value of suscScale
      suscVector <- rep(c(suscScale[k],1), c(numNode*probLowSes[i],numNode*(1-probLowSes[i])))
      
      # iterate over repeated runs of a signle scenario
      for (l in 1:numRepeat)
      {
        # initialize empty matrix
        adjMatrix <- matrix(0, nc=numNode, nr=numNode)
        
        # build up ses0 nodematched edges
        for (p in 1:numNode*probLowSes[i])
        {
          for (q in 1:numNode*probLowSes[i])
          {
            adjMatrix[p,q] <- rbinom(1,1,probEdgeWithinGroup)
          }
        }
        
        # build up ses1 nodematched edges
        for (p in (numNode*probLowSes[i] + 1):numNode)
        {
          for (q in(numNode*probLowSes[i] + 1):numNode)
          {
            adjMatrix[p,q] <- rbinom(1,1,probEdgeWithinGroup)
          }
        }
        
        # build up acorss-group edges
        for (p in (numNode*probLowSes[i] + 1):numNode)
        {
          for (q in 1:numNode*probLowSes[i])
          {
            adjMatrix[p,q] <- rbinom(1,1,probEdgeAcrossGroup)
          }
        }
        
        # let diagnal = 1 which ensures "always true" within node itself
        diag(adjMatrix) <- 1
        
        # add edge weights
        adjMatrix[1:numNode*probLowSes[i],(numNode*probLowSes[i]+1):numNode] <-
          edgeWeightScale * adjMatrix[1:numNode*probLowSes[i],(numNode*probLowSes[i]+1):numNode]
        
        adjMatrix[(numNode*probLowSes[i]+1):numNode,1:numNode*probLowSes[i]] <-
          edgeWeightScale * adjMatrix[(numNode*probLowSes[i]+1):numNode,1:numNode*probLowSes[i]]
        
        # make adj symmetric using only lower triangle 
        adjMatrix[upper.tri(adjMatrix)] <- t(adjMatrix)[upper.tri(adjMatrix)]
        
        # generate network with statnet-network 
        nwBase <- network(adjMatrix, directed=F)
        
        # network info with this adjacency matrix
        dfNwInfo <- data.frame(ID=1:numNode, SES=rep(c(0,1),c(numNode*probLowSes[i],numNode*(1-probLowSes[i]))))
        
        # attach 2 vertex attributes: node ID and node SES
        for(m in 1:2)
        {
          set.vertex.attribute(nwBase, vAttr[m], dfNwInfo[m])
        }
        
        # designate the first infected one(s)
        numInitInfect <- 1
        initNodeId <- sample(1:numNode, numInitInfect) # num or vector, depends on numInitInfect
        
        # create the dataframe that index the init node(s)
        dfIndexInitNode <- data.frame(ID=initNodeId, initInfectStatus=1)
        
        # left-join dfIndexInitNode to dfNwInfo
        dfNwInfo <- dplyr::left_join(dfNwInfo, dfIndexInitNode, by="ID", all=T)
        
        # create the 4th col in dfNwInfo: numPerNode, the node population
        dfNwInfo$numPerNode <- numPerNode
        
        # initialize states I and S for all numNode nodes
        I <- rep(0, numNode)
        S <- dfNwInfo$numPerNode
        
        # modify the states I and S considering initial infection seeded into the populaton
        I[initNodeId] <- I[initNodeId] + dfNwInfo$initInfectStatus[initNodeId]
        S[initNodeId] <- S[initNodeId] - dfNwInfo$initInfectStatus[initNodeId]
        
        # SIS over the course of timeRun
        for (t in 1:(timeRun/stepLength))
        {
          # nodewise force of infection & recovery
          nodeForceInfect <- probInterPerson * beta * suscVector * (adjMatrix %*% I) / numPerNode
          nodeForceRecov <- gamma 
          
          # nodewise probiliry of infection & recovery
          nodeProbInfect <- 1 - exp(- nodeForceInfect * stepLength)
          nodeProbRecov <- 1 - exp(- nodeForceRecov * stepLength)
          
          # sample incidence and new recovery using probability above
          incidenceByNodeMatrix[t,] <- rbinom(length(S), S, nodeProbInfect)
          newRecov <- rbinom(length(I), I, nodeProbRecov)
          
          # update states I and S
          S <- S - incidenceByNodeMatrix[t,] + newRecov
          I <- I + incidenceByNodeMatrix[t,] - newRecov
          
          # update prevalence
          prevByNodeMatrix[t,] <- I
        }
        
        # update 1 col of incidenceMatrix for this repeat
        incidenceMatrix[,l] <- rowSums(incidenceByNodeMatrix[drop=F])
        
        # update 1 col of incidenceSes0Matrix for this repeat
        incidenceSes0Matrix[,l] <- rowSums(incidenceByNodeMatrix[,1:(numNode*probLowSes[i]), drop=F])
        
        # update 1 col of incidenceSes1Matrix for this repeat
        incidenceSes1Matrix[,l] <- rowSums(incidenceByNodeMatrix[,(numNode*probLowSes[i]+1):numNode, drop=F])
        
        # update 1 col of prevMatrix for this repeat
        prevMatrix[,l] <- rowSums(prevByNodeMatrix[drop=F])
        
        # update 1 col of prevSes0Matrix for this repeat
        prevSes0Matrix[,l] <- rowSums(prevByNodeMatrix[,1:(numNode*probLowSes[i]), drop=F])
        
        # update 1 col of prevSes1Matrix for this repeat
        prevSes1Matrix[,l] <- rowSums(prevByNodeMatrix[,(numNode*probLowSes[i]+1):numNode, drop=F])
        
        # REPEAT ENDED.
      }
      
      # calculate last col of incidenceMatrix and store incidenceMatrix for this scenario
      incidenceMatrix[,numRepeat+1] <- rowMeans(incidenceMatrix[,1:numRepeat, drop=F])
      incidenceMatrices[[i*j*k]] <- incidenceMatrix
      
      # calculate last col of incidenceSes0Matrix and store incidenceSes0Matrix for this scenario
      incidenceSes0Matrix[,numRepeat+1] <- rowMeans(incidenceSes0Matrix[,1:numRepeat, drop=F])
      incidenceSes0Matrices[[i*j*k]] <- incidenceSes0Matrix
      
      # calculate last col of incidenceSes1Matrix and store incidenceSes1Matrix for this scenario
      incidenceSes1Matrix[,numRepeat+1] <- rowMeans(incidenceSes1Matrix[,1:numRepeat, drop=F])
      incidenceSes1Matrices[[i*j*k]] <- incidenceSes1Matrix
      
      # calculate last col of prevMatrix and store prevMatrix for this scenario
      prevMatrix[,numRepeat+1] <- rowMeans(prevMatrix[,1:numRepeat, drop=F])
      prevMatrices[[i*j*k]] <- prevMatrix
      
      # calculate last col of prevSes0Matrix and store prevSes0Matrix for this scenario
      prevSes0Matrix[,numRepeat+1] <- rowMeans(prevSes0Matrix[,1:numRepeat, drop=F])
      prevSes0Matrices[[i*j*k]] <- prevSes0Matrix
      
      # calculate last col of prevSes1Matrix and store prevSes1Matrix for this scenario
      prevSes1Matrix[,numRepeat+1] <- rowMeans(prevSes1Matrix[,1:numRepeat, drop=F])
      prevSes1Matrices[[i*j*k]] <- prevSes1Matrix
      
      # sum this scenario
      dfSumScenario <- data.frame(steps = 1:timeRun/stepLength,
                                  incidence = incidenceMatrix[,numRepeat+1],
                                  incidenceSes0 = incidenceSes0Matrix[,numRepeat+1],
                                  incidenceSes1 = incidenceSes1Matrix[,numRepeat+1],
                                  prev = prevMatrix[,numRepeat+1],
                                  prevSes0 = prevSes0Matrix[,numRepeat+1],
                                  prevSes1 = prevSes1Matrix[,numRepeat+1],
                                  incidenceRatio = (incidenceSes0Matrix[,numRepeat+1]/(probLowSes*numNode*numPerNode))/
                                    (incidenceSes1Matrix[,numRepeat+1]/((1-probLowSes)*numNode*numPerNode)),
                                  prevRatio = (prevSes0Matrix[,numRepeat+1]/(probLowSes*numNode*numPerNode))/
                                    (prevSes1Matrix[,numRepeat+1]/((1-probLowSes)*numNode*numPerNode)))
      
      # calculate and store the mean (and sd) of endemic prevalence over last 10 steps 
      tableAvgEP[j,k] <- 
        avgEndemicPrev <- mean(dfSumScenario$prev[(timeRun/stepLength-9):timeRun/stepLength], na.rm=T)
      tableSdEP[j,k] <-
        sdEndemicPrev <- sd(dfSumScenario$prev[(timeRun/stepLength-9):timeRun/stepLength])
      
      # calculate and store the mean (and sd) of incidence ratio over all steps  
      tableAvgIR[j,k] <-
        avgIncidenceRatio <- mean(dfSumScenario$incidenceRatio[1:timeRun/stepLength], na.rm=T)
      tableSdIR[j,k] <- 
        sdIncidenceRatio <- sd(dfSumScenario$incidenceRatio[1:timeRun/stepLength])
      
      # calculate and store the mean (and sd) of prevalence ratio over all steps   
      tableAvgPR[j,k] <-
        avgPrevRatio <- mean(dfSumScenario$prevRatio[1:timeRun/stepLength], na.rm=T)
      tableSdPR[j,k] <-
        sdPrevRatio <- sd(dfSumScenario$prevRatio[1:timeRun/stepLength])
      
      # SCENARIO ENDED. 
      # This value of relative susceptibility (this k) ENDED.
    }
    # This value of nodeMatch (this j) ENDED.
  }
  # This value of probLowSes ENDED.
  # store the EP IR PR tables into the lists
  avgTables[[3*i-2]] <- tableAvgEP
  avgTables[[3*i-1]] <- tableAvgIR
  avgTables[[3*i]] <- tableAvgPR
  sdTables[[3*i-2]] <- tableSdEP
  sdTables[[3*i-1]] <- tableSdIR
  sdTables[[3*i]] <- tableSdPR
  
  # update sceanrioKey for this value of probLowSes
  scenarioKey[3*i-2] <- paste("probLowSes",probLowSes[i],"EP")
  scenarioKey[3*i-1] <- paste("probLowSes",probLowSes[i],"IR")
  scenarioKey[3*i] <- paste("probLowSes",probLowSes[i],"PR")
}

## output lists of tables to desktop as xlsx 
write.xlsx(avgTables, "/Users/shiyaow/Desktop/0131_avgTables.xlsx", colNames=T, rowNames=T, sheetName=scenarioKey)
write.xlsx(sdTables, "/Users/shiyaow/Desktop/0131_sdTables.xlsx", colNames=T, rowNames=T, sheetName=scenarioKey)
## output lists of tables to greatlakes as xlsx
# write.xlsx(avgTables, "/home/shiyaow/0131_avgTables.xlsx", sheet=scenarioKey)
# write.xlsx(sdTables, "/home/shiyaow/0131_sdTables.xlsx", sheet=scenarioKey)

# dataframes to matrices
for(c in 1:length(scenarioKey))
{
  avgTables[[c]] <- data.matrix(avgTables[[c]])
  rownames(avgTables[[c]], nodeMatch)
  colnames(avgTables[[c]], suscScale)
  avgTablesMelted[[c]] <- melt(avgTables[[c]])
}
# melt matrices for contour plots

# create contour plots
contour <- list()
for (c in 1:length(scenarioKey))
{
  contour[[c]] <-
    ggplot(avgTablesMelted[[c]], aes(x=Var1, y=Var2, z=value)) +
    geom_point(aes(color=value)) +
    stat_contour(na.rm=T) +
    labs(x="assortative mixing level", y="relative susceptibility ses0 vs ses1", title=scenarioKey[c])
}

# store "length(scenarioKey)" plots with pdf device
pdf("/Users/shiyaow/Desktop/0131_contours.pdf")
# pdf("/home/shiyaow/0131_contour.pdf")
contour
dev.off()




