#### Single locus allele frequency simulation

## This function simulates constant selection (a != 0), fluctuating selection (b != 0) and neutral evolution (a and b = 0)
## populations = number of populations
## generations = total number of generations, including initial
## p0 = initial allele frequency
## popsize = population size
## samplesize = sample size for observed allele freq, should be <= popsize
## h = het. effect
## a = constant selection parameter
## b = fluctuating selection parameter
## model is s = a + b*env
Sim <- function(populations = 10, generations = 11, p0 = 0.5, popsize = 200,
                                  samplesize = 30, h = 0.5, a = 0, b = 0){
    # population number
    nPol <- populations
    # number of generations
    nGen <- generations

    cat("a =",a,";b = ",b,"\n")

    # Initialize the matrix
    evoSim <- matrix(NA, nrow = nGen, ncol = 2*nPol)
    # create a vector of odd numbers (for populations) and even numbers
    # (environmental data)
    popCol <- seq(from = 1, to = 2*nPol - 1, by = 2)
    envCol <- seq(from = 2, to = 2*nPol, by = 2)
    # fill in row one with starting allele frequency and environmental data
    evoSim[1, popCol] <- 0.5

    ##Simulation

    ## sample the environment
    for (env in envCol) {
        for (gen in 1:nGen) {
            evoSim[gen, env] <- runif(1, 0, 1)
        }
    }

    for (pop in popCol){
        for (gen in 2:nGen) { #set selection equal to a function of environmental data from the previous generation
            s <- a + (b*((evoSim[gen - 1, pop + 1]) -0.5 ))## center environment for simulations
            freq <-  evoSim[gen - 1 , pop] #save the allele frequency for the previous generation
            expectedFreq <- freq + (s*freq*(1 - freq)*(freq + h*(1-2*freq))) #calculate expected frequency
            newFreq <- rbinom(1, (2*popsize), expectedFreq)/(2*popsize) #genetic drift based on population size
            evoSim[gen, pop] <- newFreq
        }
        ## sample size perturbs allele frequences from true values
        evoSim[,pop]<-rbinom(nGen,(2*samplesize),evoSim[,pop])/(2*samplesize)
        ## arcsine square root transformation
        ## evoSim[,pop]<-2 * asin(sqrt(evoSim[,pop]))
    }
    

    ## recode as change in p, but scaled to be between 0, and 1
    evoSimNew <- matrix(NA, nrow = nGen - 1, ncol = 2*nPol)

    for (env in envCol) {
        for (gen in 1:nGen - 1) {
            evoSimNew[gen, env] <- evoSim[gen, env]
        }
    }

    for (pop in popCol){
        for (gen in 1:nGen-1){ 
            ## pi is range of p, so 2 pi is range of delta p
            ##evoSimNew[gen, pop] <- (evoSim[gen + 1, pop] - evoSim[gen, pop])/(2*pi) + .5 
            ## SD standardize and scale to 0 to 1 
            evoSimNew[gen, pop] <- (evoSim[gen + 1, pop] - evoSim[gen, pop])/sqrt(evoSim[gen, pop] * (1-evoSim[gen, pop]))
            evoSimNew[gen, pop] <- (evoSimNew[gen, pop] + 1)/2
        }
    }

    return(evoSimNew)
}
#generate 5000 data sets for strong, constant selection, samples a ~ +- .2
for (i in 1:5000){
  data <- Sim(10,11,0.5,popsize=200,samplesize=200,h=0.5,a=rbeta(1,10,40)*sample(c(-1,1),1),b=0)
  write.table(data, paste0("Type1DataSet_", i,".csv"),
            col.names = F, row.names = F)
}
#generate 5000 data sets for strong, fluctuating selection, a = 0 (no directional), b ~ +- .3
for (i in 1:5000){
  data <- Sim(10,11,0.5,popsize=200,samplesize=200,h=0.5,a=0,b=rbeta(1,15,35)*sample(c(-1,1),1))
  write.table(data, paste0("NewType2DataSet_", i,".csv"),
            col.names = F, row.names = F)
}
#generate 5000 data sets for pure drift
for (i in 1:5000){
  data <- Sim(10,11,0.5,popsize=200,samplesize=200,h=0.5,a=0,b=0)
  write.table(data, paste0("NewType3DataSet_", i,".csv"),
            col.names = F, row.names = F)
}

