#################################################
####      Group management v23 blank        #####
#################################################
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

##########################
#overview of the workflow#
##########################

## 0. Set parameters
## 1. Create parameter space (ps)
## 2. Run model
## 3. Generate overall dataset and plots
## 4. Statistical analysis of the simulation outputs (linear mixed-effect model)

## 0. Set parameters

# biological parameters
survival <-            # overarching survival rate
sex_ratio_N <-         # sex ratio of the founders
age_max <-             # maximum longevity of the species
first_repro <-         # age at which an individual of this species reproduces for the first time
last_repro <-          # age at which an individual of this species reproduces for the first time

mating_system <-       # "M": monogamous; "P": polygynous; "Ha": harem
harem_size <-          # number of females per male
stdev_harem_size <-    # standard deviation harem size
litter_size <-         # average litter size
sd_litter_size <-      # standard deviation litter size

# run parameters
generations <-         # number of generations you want to run the simulation
simulations <-         # number of simulations you want to run

# population structure parameters population
NG_BC <-                      # number of groups in a breeding circle
NBC_I <-                      # number of breeding circles per institute
NI <-                         # number of institutes

NG <- NG_BC*NBC_I*NI          # total number of groups

holding_populations <-        # 1 =  present; 0 = absent
restocking <-                 # 1 =  yes; 0 =  no
n_restock <-                  # number of individuals taken from holding pop to restock crashed pops
years_in_holding <-           # maximum number of years an individual can stay in a holding population, before it is removed from the system
level_holding_pop <-          # holding populations can be absent or installed at the breeding circle, institute or global level: "non","BC","I","global"
K_HP <-                       # maximum size ("carrying capacity") of a holding population

# molecular data parameters
mol_data <-          # 'I': random data each founder unique; 'II' : real data manual input
no_markers <- 10
store_history <- 0   # 1 = keep track of all historic individuals; 0 = only keep track of living individuals

# transfer parameters
transfer_age <-         # list the ages at which animals may be transferred (e.g., if c(1), only individuals with an age of one year will be transferred)
transfer_gender <-      # list the sexes that may be transferred (e.g., if c(1), only males will be transferred)

# kinship parameters
kinships <-          # 1 if you want to calculate a kinship matrix
kin_method <-        # method of calculating kinship: 1 = pedigree based; 2 = based on molecular data

# define the combination of metapopulation designs and transfer strategies to be tested
GroupSize <-         # specify the minimum and maximum group size to be simulated, and the step
Transfer_Strategy <- c()  # specify the breeding strategies you want to simulate (options: "LINE_BREEDING": line breeding;"RANDOM": random breeding;"CB": circular breeding; "FALCONERS": falconers' breeding scheme; "MAI": maximum avoidance of inbreeding;"MK": mean kinship breeding)
TFG <- seq()    # specify transfer frequency between groups
TFHL <- seq()   # specify the transfer frequency ratio at the higher level
FounderSize <- c()   # specify the founder sizes (= number of individuals in each group at the start) you want to simulate

scenario <- paste(NG_BC,"*",NBC_I,"*",NI,sep="")   # this line generates the scenario a name, used later when plotting

## 1. Create parameter space

parameter_space <- expand.grid(GroupSize,TFG,TFHL,FounderSize)    # holds all combinations to be simulated
colnames(parameter_space) <- c("K_G","TFG","TFHL","FounderSize")

## 2. Run model

for (ps in c(1:nrow(parameter_space)))        # each combination of parameters will be simulated
{print(paste("ps analysis nr",ps,sep=" "))    # to keep track of progress while running the simulation
 TransFrequencyGroups<-parameter_space$TFG[ps]
 TransRatioHigherLevel<-parameter_space$TFHL[ps]

  if(NG_BC>1) {TFG<-seq(TransFrequencyGroups,generations,TransFrequencyGroups)}else {TFG<-generations+1}
  if(NBC_I>1 & generations > TransFrequencyGroups*TransRatioHigherLevel) {TFBC<-seq(TransFrequencyGroups*TransRatioHigherLevel,generations,TransFrequencyGroups*TransRatioHigherLevel)} else {TFBC<-generations+2}
  if(NI>1 & generations > TransFrequencyGroups*TransRatioHigherLevel*TransRatioHigherLevel ) {TFI<-seq(TransFrequencyGroups*TransRatioHigherLevel^2,generations,TransFrequencyGroups*TransRatioHigherLevel^2)} else {TFI<-generations+3}
  TFG<-TFG[!TFG%in%c(TFBC,TFI)] 
  TFBC<-TFBC[!TFBC%in%TFI]

K_G <- parameter_space$K_G[ps]            # maximum size ("carrying capacity") of a group
K <- c(rep(K_G,NG),K_HP)                  # a vector listing the carrying capacity of all groups in order, including the holding population at the end
NFG <- parameter_space$FounderSize[ps]    # number of founders per group
# end defining run specific parameters

################
#molecular data#
################

if(mol_data=="I")  # each founder unique
{no_alleles_marker <- 2*NFG*NG
  loci <- vector(mode='list', length=no_markers)
  probs <- vector(mode='list', length=no_markers)
  for (marker in c(1:no_markers))
       {loci[[marker]] <- seq(100,100+no_alleles_marker,1)
        probs[[marker]] <- rep(1/no_alleles_marker,no_alleles_marker)}
 }
 
if(mol_data=="II") # founder genotypes based on predefined allele frequencies
{m_1 <- c(219,221,223,225,227,229,233)
 prob_m_1 <- c(0.368733,0.003646,0.19918,0.391978,0.034184,0.001367,0.000912)
 m_2 <- c(100,110)
 prob_m_2 <- c(0.35,0.65)
 m_3 <- c(100,100)
 prob_m_3 <- c(0.5,0.5)

loci <- list(m_1,m_2,m_3)
probs <- list(prob_m_1,prob_m_2,prob_m_3)
no_markers <- length(loci)}

################################
##generate overarching dataset##
################################

out <- vector(mode='list', length=6)
names(out) <- c("pop_str","K_G","TransferStrategy","Simulation","TFG","TFHL")
out_GD <- vector(mode='list')
out_N <- vector(mode='list')
out_Fis <- vector(mode='list')
out_Cum_crash <- vector(mode='list')
out_Cum_transfer <- vector(mode='list')
out_Cum_pop_rescue <- vector(mode='list')
out_MMK <- vector(mode='list')

##############################
##define founding population## -> to be saved for later analysis
##############################

# define founder population structure     # alternatively use sample(c(1,2),N,replace=T,prob=c(sex_ratio_N,1-sex_ratio_N))
NF <- NFG*NG     # total number of founders (all groups combined)

output_parameter <- 0

for(TS in Transfer_Strategy)
{print(TS)
  ########################
  ##OUTPUT TABLES per TS##
  ########################
  GD <- vector(mode='list', length=generations)        # initialize a vector that will keep track of the genetic diversity in each generation
  Fis <- vector(mode='list', length=generations)       # initialize a vector that will keep track of the inbreeding coefficient in each generation
  N <- vector(mode='list', length=generations)         # initialize a vector that will keep track of the size of a group in each generation
  Cum_transfers <- rep(0,generations)     # initialize a vector that will hold the number of transfers in each generation
  Cum_pop_crash <- rep(0,generations)     # initialize a vector that will hold the number of groups crashing in each generation
  Cum_pop_rescue <- rep(0,generations)    # initialize a vector that will hold the number of groups rescued by re-establishment from the holding population in each generation
  MMK_gen <- rep(NA,generations)
  
  for (sim in 1:simulations)              # this loops make sure that you run the amount of replicates per combination that you specified
  { print(sim)                            # so you can keep track of where you are while running the script
    output_parameter <- output_parameter+1
    
    nr_migrants <- 0     # initialize a count that will keep track of the number of individuals transferred during a simulation
    nr_rescue <- 0       # initialize a count that will keep track of the number of rescues during a simulation
    nr_crash <- 0        # initialize a count that will keep track of the number of (group) crashes during a simulation
    
    step_G <- 0          # do not change: to be used in the falconers' breeding scheme and maximum avoidance of inbreeding
    step_BC <- 0         # do not change: to be used in the falconers' breeding scheme and maximum avoidance of inbreeding
    step_I <- 0          # do not change: to be used in the falconers' breeding scheme and maximum avoidance of inbreeding
    
    founder <- data.frame(id=c(1:NF),
                        gender=rbinom(NF,1,prob=sex_ratio_N),
                        age=rpois(NF,2),
                        sire=0,
                        dam=0,
                        alive=1,
                        birthyear=0,
                        deathyear=NA,
                        MK=0,
                        G_Unique=sort(rep(c(1:NG),NFG))-1,
                        G=rep(rep(sort(rep(c(1:NG_BC),NFG)),NBC_I),NI)-1,
                        BC=rep(sort(rep(c(1:NBC_I),NFG*NG_BC)),NI)-1,
                        I=sort(rep(c(1:NI),NFG*NG_BC*NBC_I))-1,
                        GBirth=rep(rep(sort(rep(c(1:NG_BC),NFG)),NBC_I),NI)-1,
                        Homo=0,
                        t_in_holding=0
                        )    # make a dataframe holding the founder population, with one row holding the information for every animal in this population
    
    # add molecular data to the founder population dataframe
    if(mol_data=="II")
    {for (marker in 1:no_markers)
    {founder[,paste("m",marker,"1",sep="_")]<-sample(loci[[marker]],NF,replace=T,prob=probs[[marker]])
     founder[,paste("m",marker,"2",sep="_")]<-sample(loci[[marker]],NF,replace=T,prob=probs[[marker]])}}
    
    if(mol_data=="I")
    {for (marker in 1:no_markers)
    {founder[,paste("m",marker,"1",sep="_")]<-seq(100,100+no_alleles_marker,2)[1:NF]
     founder[,paste("m",marker,"2",sep="_")]<-seq(101,101+no_alleles_marker,2)[1:NF]}}
    
    # create kinship matrix
    if (kinships==1){kinship<-diag(NF)/2;colnames(kinship) <- founder$id;rownames(kinship)<-founder$id}
    
    # define dynamic population structure
    current <- founder    # defines the table 'current' that will keep track of the current population throughout the simulation (it starts of course from the founder population)
    
    first_m <- min(grep("m_",colnames(current))) 
    last_m <- max(grep("m_",colnames(current)))-1
    
    death <- current[current$alive==0,]    # makes a table that will keep track of all dead individuals, if you chose to store history
    id_nr <- max(founder$id)
    
    # simulate population dynamics by looping over the specified number of generations to simulate
    for(j in (1:generations))
    {
     # breeding in a Group loop
     for (gr in c(1:NG)-1)        # loop through all groups
      {current_pop <- filter(current, G_Unique==gr & alive ==1)    # take into account only living animals in the specific group we are now working with
       males <- filter(current_pop,age>=first_repro & age<=last_repro & gender==1 & alive==1)$id    # all males that can reproduce during this generation
       males <- males[sample(length(males))]          # number of males that can reproduce during this generation
       females <- filter(current_pop,age>=first_repro & age<=last_repro & gender==0 & alive==1)$id  # all females that can reproduce during this generation
       females <- females[sample(length(females))]    # number of males that can reproduce during this generation
       
       c2 <- max(0,round((K[gr+1]-nrow(current_pop))/litter_size))
        couples <- min(length(males),length(females),c2)    # number of couples that will reproduce during this generation

       # generate random harem contributions 'hs' listing males for all couples being created (multiple females with the same male)
       if (couples>0)
       { hs <- rnorm(couples,harem_size,stdev_harem_size)
         hs <- round(hs,0)
         hs[hs<1] <- 1
         harem_males <- NULL
           for (m2 in c(1:length(hs))){harem_males <- c(harem_males,rep(males[m2],hs[m2]))}}
       
       # let the couples produce offspring
       while (couples>0)
       {ls1 <- max(1,round(rnorm(1,litter_size,sd_litter_size),0))  #litter size
        dam <- females[couples]
        col_dame_i <- which(colnames(kinship)%in%dam)
        # which males are coupled
         if(mating_system=="M")  {sire<-rep(sire<-males[couples],ls1)} 
         if(mating_system=="P")  {sire<-as.numeric(sample(as.character(males),ls1,replace=T))}
         if(mating_system=="Ha") {sire<-rep(harem_males[couples],ls1)}
       for (k in (1:ls1)) 
       {b <- nrow(current)+1
        id_nr <- id_nr+1
        col_sire_i <- which(colnames(kinship)%in%sire[k])  
        
        sire2 <- current_pop[current_pop$id==sire[k],]
        dam2 <- current_pop[current_pop$id==dam,]
        
        current[b,]<-c(id_nr,rbinom(1,1,sex_ratio_N),0,sire[k],dam,1,j,NA,0,gr,sire2$G,sire2$BC,sire2$I,gr,0,0,rep(0,2*no_markers))
        current[b,seq(first_m,last_m,2)]<-sire2[,seq(first_m,last_m,2)+sample(c(0,1),replace=T,no_markers)]
        current[b,seq(first_m+1,last_m+1,2)]<-dam2[,seq(first_m,last_m,2)+sample(c(0,1),replace=T,no_markers)]
       
       # update kinship matrix (after breeding)
       if (kinships==1)
       {kinship_i <- colSums(kinship[c(col_dame_i,col_sire_i),])*0.5
        kinship <- cbind(kinship,kinship_i)
        kinship_i <- c(kinship_i,0.5*(1+kinship[col_dame_i,col_sire_i]))
        kinship <- rbind(kinship,kinship_i)
        colnames(kinship)[ncol(kinship)] <- id_nr
        rownames(kinship)[nrow(kinship)] <- id_nr
        current$Homo[b] <- kinship[col_dame_i,col_sire_i]}
       } 
       couples <- couples-1
       } # END COUPLES
      } # END Group loop

      ##############################################
      #realized kinship calculations between groups#
      ##############################################
      if(2 %in%  kin_method)    # this means: when approaching kinship calculation by looking at the molecular data (instead of pedigree)
      {CNM <- grep("^m_.*\\_1$",colnames(current))
       current_G <- sort(unique(current$G_Unique))
       if(length(current_G)>1)
       {MMK <- matrix(0,length(current_G),length(current_G))  # Molecular Mean Kinship
         for(cnm in CNM)
           {AF<-data.frame(alleles=unique(c(current[,cnm],current[,cnm+1])))
           # create allele frequencies for group
           for (g2 in current_G)
             {a <- current %>% filter(G_Unique == g2) %>% select(cnm,cnm+1) %>% unlist() %>% table()
              a <- data.frame(a/sum(a))
              names(a) <- c('alleles',g2)
              AF <- merge(AF,a,by="alleles",all=T)}
      AF[is.na(AF)]<-0
      
      for (k1 in c(1:length(current_G)))
        for (k2 in c(k1:length(current_G)))
        {MMK[k1,k2]<-MMK[k1,k2]+sum(AF[,k1+1]*AF[,k2+1])}
      
      }
      MMK <- MMK/no_markers
      MMK[lower.tri(MMK)] <- t(MMK)[lower.tri(t(MMK))]
      colnames(MMK) <- current_G
      rownames(MMK) <- current_G  
       } else {MMK <- 'na'}
      }

      #####################################
      #transfers###########################
      #####################################
      migrants <- which(current$age%in%transfer_age & current$gender%in%transfer_gender & current$G_Unique%in%c(0:(NG-1)))    # specifies which individuals will be transferred
      
      if(TS=="CB")          # transferring according to the circular breeding strategy
      {if(j %in%TFG)        # once every TFG generations...
      {current$G[migrants] <- (current$G[migrants]+1)%%NG_BC}      # ... transfer migrants to the next group in the breeding circle
        if(j %in%TFBC)      # once every TFBC generations...
        {current$BC[migrants] <- (current$BC[migrants]+1)%%NBC_I}  # ... transfer migrants to the next breeding circle/life support system
        if(j %in%TFI)       # once every TFI generations...
        {current$I[migrants] <- (current$I[migrants]+1)%%NI}       # ... transfer migrants to the next institute
      }
      if(TS=="RANDOM")       # transferring according to the random breeding strategy
      {if (j %in%TFG)        # once every TFG generations...
      {if(NG_BC>2) {current$G[migrants] <- (current$G[migrants]+sample(c(1:(NG_BC-1)),length(migrants),replace=T))%%NG_BC}             # .... if there are more than two groups per breeding circle/life support system, move migrants to a random group
        if(NG_BC==2){current$G[migrants] <- (current$G[migrants]+1)%%2}}       # .... if there are two groups per breeding circle/life support system only, move migrants to the other group
        if (j %in%TFBC)      # once every TFBC generations...
        {if(NBC_I>2) {current$BC[migrants] <- (current$BC[migrants]+sample(c(1:(NBC_I-1)),length(migrants),replace=T))%%NBC_I}         # .... if there are more than two breeding circles/life support systems per institute, move migrants to a random breeding circle/life support system
          if(NBC_I==2){current$BC[migrants] <- (current$BC[migrants]+1)%%2}}   # .... if there are two breeding circles per institute only, move migrants to the other breeding circle
        if (j %in%TFI)       # once every TFI generations...
        {if(NI>2) {current$I[migrants] <- (current$I[migrants]+sample(c(1:(NI-1)),length(migrants),replace=T))%%NI}                    # .... if there are more than two institutes, move migrants to a random institute
          if(NI==2){current$I[migrants] <- (current$I[migrants]+1)%%2}}        # .... if there are two institutes only, move migrants to the other institute
      }
      if(TS=="FALCONERS")      # transferring according to the falconers' breeding scheme
      {if (j %in%TFG)          # once every TFG generations...
      {step_G <- step_G+1      # ... in the falconers' breeding scheme, one is added to the step in each round (see explanation and graphical representation in the READ.ME file)
       step_G <-(step_G%%NG_BC>0)*step_G+(step_G%%NG_BC==0)        # ... when "the circle is round", the step must reset
       current$G[migrants] <- (current$G[migrants]+step_G)%%NG_BC}         # ... transfer individuals within a breeding circle/life support system according to the appropriate step
        if (j %in%TFBC)        # once every TFBC generations...
        {step_BC <- step_BC+1  # ... in the falconers' breeding scheme, one is added to the step in each round (see explanation and graphical representation in the READ.ME file)
         step_BC <-(step_BC%%NBC_I>0)*step_BC+(step_BC%%NBC_I==0)  # ... when "the circle is round", the step must reset
        current$BC[migrants] <- (current$BC[migrants]+step_BC)%%NBC_I}     # ... transfer individuals within an institute according to the appropriate step
        if (j %in%TFI)         # once every TFI generations...
        {step_I <- step_I+1    # ... in the falconers' breeding scheme, one is added to the step in each round (see explanation and graphical representation in the READ.ME file)
         step_I <-(step_I%%NI>0)*step_I+(step_I%%NI==0)            # ... when "the circle is round", the step must reset
         current$I[migrants] <- (current$I[migrants]+step_I)%%NBC_I}       # ... transfer individuals between institutes according to the appropriate step
      }
      if(TS=="MAI")            # transferring according to the maximum avoidance of inbreeding strategy
      {if (j %in%TFG)          # once every TFG generations...
      {step_G <- (2^step_G%%NG_BC>0)*step_G+(2^step_G%%NG_BC==0)*0         # ... calculate the step according to the formula 2^[t-1]
       current$G[migrants] <- (current$G[migrants]+2^step_G)%%NG_BC        # ... transfer individuals within a breeding circle/life support system according to the appropriate step
       step_G <- step_G+1}     # ... add one to t, so the next step will be calculated properly
        if (j %in%TFBC)        # once every TFBC generations...
       {step_BC <- (2^step_BC%%NBC_I>0)*step_BC+(2^step_BC%%NBC_I==0)*0    # ... calculate the step according to the formula 2^[t-1]
        current$BC[migrants] <- (current$BC[migrants]+2^step_BC)%%NBC_I    # ... transfer individuals within an institute according to the appropriate step
        step_BC <- step_BC+1}  # ... add one to t, so the next step will be calculated properly
        if (j %in%TFI)         # once every TFI generations...
        {step_I <- (2^step_I%%NBC_I>0)*step_I+(2^step_I%%NBC_I==0)*0       # ... calculate the step according to the formula 2^[t-1]
        current$I[migrants] <- (current$I[migrants]+2^step_I)%%NBC_I       # ... transfer individuals between institutes according to the appropriate step
        step_I <- step_I+1}    # ... add one to t, so the next step will be calculated properly
      }
      if(TS=="MK" & j %in%c(TFG,TFBC,TFI))           # transferring according to the mean kinship breeding strategy
      {G_unique_kinships <- current$G_Unique[current$id%in%colnames(kinship)]     # unique groups for all individuals listed in kinship matrix
       for(m in migrants)  # m is line number migrant
        {migrant <- current$id[m]
         G_unique <- current$G_Unique[m]   # unique group number migrant
         Inst <- G_unique%/%(NG_BC*NBC_I)
         z <- which(G_unique_kinships%in%G_unique)  # line numbers of individuals group in kinship
         MKm<-kinship[which(rownames(kinship)%in%migrant),]   
         MKsource <- mean(kinship[z,z])
      
         kinship_selector <- data.frame(Group=c(1:NG)-1,mk_m_g=NA)
           for (x in kinship_selector$Group) {kinship_selector$mk_m_g[x+1]<-mean(MKm[G_unique_kinships%in%x])}
           if(kinship_selector$mk_m_g[kinship_selector$Group%in%G_unique]>mean(kinship[z,z]))   #only move when mean kinship is above average MK
    
      {
        if(j%in%TFG)     # alternative groups in breeding circle
        {first<-(G_unique%/%NG_BC)*NG_BC                                                  # first group id of breeding circle
         alternative_groups <- c(first:(first+NG_BC-1))                                   # all groups within breeding circle
         alternative_groups <- alternative_groups[-which(alternative_groups%in%G_unique)] # remove current group from list 
         ks2 <- kinship_selector[alternative_groups+1,]                                   # reduce matrix to selected groups, +1 since population number ranges from 0 to x as such differing 
         ks2 <- ks2[sample(1:nrow(ks2)),]                                                 # shuffle to avoid repeated transfers to the same group 
         current$G_Unique[m] <- ks2$Group[order(ks2$mk_m_g)][1]}                          # take first element of sorted matrix
        
        if(j%in%TFBC)    # alternative groups in institutes
        {alternative_groups <- c(0:((NG_BC*NBC_I)-1))+Inst*NG_BC*NBC_I
         alternative_groups <- alternative_groups[-which(alternative_groups%in%G_unique)] # remove current group from list 
         ks2 <- kinship_selector[alternative_groups+1,]                                   # reduce matrix to selected groups, +1 since population number ranges from 0 to x as such differing
         ks2 <- ks2[sample(1:nrow(ks2)),]                                                 # shuffle to avoid repeated transfers
         current$G_Unique[m] <- ks2$Group[order(ks2$mk_m_g)][1]}                          # take first element of sorted matrix
             
        if(j%in%TFI)
        {alternative_groups <- c(0:(NG-1))
        alternative_groups <- alternative_groups[!alternative_groups%in%(c(0:((NG_BC*NBC_I)-1))+Inst*NG_BC*NBC_I)]            # remove current group from list 
        ks2<-kinship_selector[alternative_groups+1,]                                      # reduce matrix to selected groups, +1 since population number ranges from 0 to x as such differing
        ks2<-ks2[sample(1:nrow(ks2)),]                                                    # shuffle to avoid repeated transfers
        current$G_Unique[m] <- ks2$Group[order(ks2$mk_m_g)][1]}                           # take first element of sorted matrix
        
      }
      }
      current$G <- current$G_Unique%%NG_BC 
      current$BC <- (current$G_Unique%%(NG_BC*NBC_I))%/%NG_BC
      current$I <- current$G_Unique%/%(NG_BC*NBC_I)
      }
      
      current$G_Unique <- current$G+(current$BC)*NG_BC+NG_BC*NBC_I*(current$I)
      current[migrants,] %>% filter(GBirth!=G_Unique) %>% nrow()-> new_migrants
      nr_migrants <- nr_migrants+new_migrants    # update 'nr_migrants' that keeps track of how many transfers have taken place
      
      ########
      #ageing#
      ########
      current$age <- current$age+1    # all individuals have grown one year older
      current$t_in_holding[current$G_Unique>=1000] <- current$t_in_holding[current$G_Unique>=1000]+1    # all individuals in a holding population have spent one more year in a holding population

      # list all culled individuals in the new 'death_gen' table
      death_gen <- current[current$age>age_max | current$t_in_holding>years_in_holding,] # specifies that all individuals older than the maximum age, or longer in a holding population than desirable, are to be culled
      current <- current[!current$id%in%death_gen$id,] # these individuals die (removed from current population)
      current$alive <- rbinom(nrow(current),1,survival)   # kill or cull individuals
      death_gen <- rbind(death_gen,current[current$alive==0,])
      current <- current[current$alive==1,]
        if(nrow(death_gen)>0){death_gen$deathyear<-j}
      # A. pinpoint surplus individuals and move them to holding if relevant
      for (hold_inst in c(0:(NI-1)))
        {sources <- table(current$G_Unique[current$I%in%(hold_inst)])
         sources <- sources[which(sources>K_G)]
         cull <- NULL
          while(length(sources)>0)
           {cull <- c(cull,sample(current$id[current$G_Unique==names(sources)[1]],sources[1]-K_G))
            sources <- sources[-1]}
          if(holding_populations==1 & length(cull)>0 & sum(current$G_Unique==NG)<K_HP)
          # put ids in holding population based on available "space"
           {cull_to_holding<-sample(cull,min(length(cull),K_HP-sum(current$G_Unique==NG)))
            cull <- cull[!cull%in%cull_to_holding] # remove individuals that will go to a holding population from the table of individuals to be culled
            current$G_Unique[current$id%in%cull_to_holding] <- (1000+hold_inst)} # specify that these individuals are now in the appropriate holding population
          death_gen <- rbind(death_gen,current[current$id%in%cull,])
          current <- current[!current$id%in%cull,] # cull individuals by updating the table of current population, so all culled individuals are removed

      # B. restock Unique groups from specific holding
      to_restock <- c(0:(NG-1))[!c(0:(NG-1)) %in% unique(current$G_Unique)]
      size_holding <- sum(current$G_Unique==1000+hold_inst)
      while (size_holding>0 & restocking == 1 & length(to_restock)>0)
       {nr_rescue <- nr_rescue + 1    # keeps track of the number of times rescue (restocking a group from the holding population) takes place
        current$G_Unique[sample(which(current$G_Unique==(1000+hold_inst)),min(size_holding,n_restock),replace=F)] <- as.numeric(sample(as.character(to_restock),1))
        to_restock <- c(0:(NG-1))[!c(0:(NG-1)) %in% unique(current$G_Unique)]
        size_holding <- sum(current$G_Unique==(1000+hold_inst))}
      
    }
      
      # C. clean kinship matrix based on deceased individuals
      if(kinships==1 & nrow(death_gen)>0 & length(kinship)>1)
      {kinship <- kinship[-which(colnames(kinship)%in%death_gen$id),-which(colnames(kinship)%in%death_gen$id)]}
      
      
      if(store_history==1){death <- rbind(death,death_gen)}    # if a history of dead individuals is stored, add the animals that died during this generation to that inventory
      death_gen <- death_gen[numeric(0), ]    # empty the dataframe 'death_gen'
      
      Cum_pop_rescue[j] <- nr_rescue          # for the current generation, update the vector 'Cum_pop_rescue' that holds the number of groups rescued by re-establishment from the holding population in each generation
      current$t_in_holding[current$G_Unique<1000] <- 0
      
      #################################################
      ##create output file#############################
      #################################################
      GD_temp <- rep(NA,NG)        # initialize a vector that holds genetic diversity of all groups
      Fis_temp <- rep(NA,NG)       # initialize a vector that holds inbreeding coefficient of all groups
      N_temp <- rep(NA,NG)         # initialize a vector that holds size of all groups
      
      groups_to_include <- sort(unique(current$G_Unique))
      groups_to_include <- groups_to_include[groups_to_include<1000] # groups to include here excluded holding population (that have a number over 1000)
      for (ds in groups_to_include)   # loop through all groups except the holding population(s)
      {ids_g <- current$id[current$G_Unique%in%(ds)]  # here -1 because groups range from 0 to NG-1
       N_temp[ds+1] <- length(ids_g)
      if(length(kinship)>1)
        {k_all <- kinship[which(rownames(kinship)%in%ids_g),]
         k_gr <- kinship[which(rownames(kinship)%in%ids_g),which(colnames(kinship)%in%ids_g)]
         GD_temp[ds+1] <- mean(k_all)                            # genetic diversity (GD)  based on kinship values
         Fis_temp[ds+1] <- 2*(mean(diag(k_gr))-0.5)}             # inbreeding coefficient (Fis) based on kinship values
      }
      
      N[[j]] <- c(N_temp,nrow(current))                   # last element is population level average
      GD[[j]] <- c(GD_temp,mean(kinship))                 # last element is population level average
      Fis[[j]] <- c(Fis_temp,2*(mean(diag(kinship))-0.5)) # last element is population level average
      Cum_transfers[j] <- nr_migrants                     # update the vector 'Cum_transfers' that holds the number of transfers in each generation
      MMK_gen[j] <- mean(MMK,na.rm=T)                     # genetic diversity (GD) based on molecular data
      Cum_pop_crash[j] <- nr_crash                        # for the current generation, update the vector 'Cum_pop_Crash that holds the number of groups crashing in each generation
      

    } # END generations
    
    # generate output figures while running script
    d1<-data.frame(gen_div=1-sapply(GD,"[[",NG+1),Gen=1:generations)
    d1<-cbind(d1,GD_mk=1-MMK_gen)
    
    # the following figure shows genetic diversity (GD) throughout the simulated generations
    fig1 <- ggplot(data=d1,aes(x=Gen))+
      geom_line(aes(y=gen_div),colour="red")+
      geom_line(aes(y=GD_mk),colour="blue")+
      xlab("Generations")+
      ylab("GD")+
      lims(y = c(min(d1$gen_div,d1$GD_mk),1))+ 
      theme(legend.position = "none")
    
    d1<-cbind(d1,Fis=sapply(Fis,"[[",NG+1))
    # the following figure shows inbreeding coefficient (Fis) throughout the simulated generations
    fig2 <- ggplot(d1,aes(x=Gen,y=Fis,color=which(Transfer_Strategy%in%TS)))+
      geom_line()+
      xlab("Generations")+
      ylab("Fis")+
      lims(y = c(0,max(d1$Fis)))+ 
      theme(legend.position = "none")
    
    # if you want, you can also include more plots that show the cumulative number of transfers and reszcues over the simulated generations (by removing the hashtags)
    #fig3 <- plot(Cum_transfers~c(1:generations),type="l",ylim=c(0,max(200,max(Cum_transfers))),col=which(Transfer_Strategy%in%TS),xlab="generation",ylab="Cum. N. transfers")
    #fig4 <- plot(Cum_pop_rescue~c(1:generations),type="l",ylim=c(0,max(200,max(Cum_pop_rescue))),col=which(Transfer_Strategy%in%TS),xlab="generation",ylab="Cum. N. rescues")
    
    N2 <- t(array(unlist(N), dim = c(NG+1,generations)))
    N2 <- data.frame(rep(1:generations,NG+1),sort(rep(1:(NG+1),generations)),c(N2))
    colnames(N2) <-c ('Gen','group','N')
    N2$group <- as.factor(N2$group)
    # the following figure shows the group sizes of the different simulated groups throughout the simulated generations
    fig5 <- ggplot(data=N2[N2$group %in% c(1:NG),],aes(x=Gen,y=N,group=group,color=group))+ #NG+1 for holding population
      geom_line()+
      xlab("Generations")+
      theme(legend.position = "none")
   grid.arrange(fig1,fig5,fig2, ncol = 1, nrow = 3, top = paste(TS," SCENARIO ",NG_BC,"-",NBC_I,"-",NI,"(",paste(parameter_space[ps,], collapse = '_'),")",sep=""))
   
    ################################
    ##complete overarching dataset##   OUT dataset
    ################################
    out$pop_str[output_parameter] <- paste(NG_BC,NBC_I,NI,sep="_")
    out$K_G[output_parameter] <- K_G
    out$TransferStrategy[output_parameter] <- TS
    out$Simulation[output_parameter] <- sim     
    out$TFG[output_parameter] <- TransFrequencyGroups
    out$TFHL[output_parameter] <- TransRatioHigherLevel    
    out_GD[[output_parameter]] <- t(array(unlist(GD), dim = c(NG+1,generations)))
    out_N[[output_parameter]] <- t(array(unlist(N), dim = c(NG+1,generations)))
    out_Fis[[output_parameter]] <- t(array(unlist(Fis), dim = c(NG+1,generations)))
    out_Cum_transfer[[output_parameter]] <- unlist(Cum_transfers)
    out_MMK[[output_parameter]] <- unlist(MMK_gen)
    out_Cum_pop_rescue[[output_parameter]] <- unlist(Cum_pop_rescue)
    } # END simulations
}# END Transfer strategies

## 3. Generate overall dataset and plots

#####################
##overarching plots##
#####################

plot_data <- expand.grid(ps,paste(NG_BC,"_",NBC_I,"_",NI,sep=""),Transfer_Strategy,seq(10,generations,10),c(1:simulations),c(1:(NG+1)),1,0,0,0)
colnames(plot_data) <- c('ps_nr','Scenario','TS','Gen','sim','Group','GD','Fis','N','Transfers')

for (i in c(1:nrow(plot_data)))
{el <- which(out$TransferStrategy%in%plot_data$TS[i] # el element
        & out$Simulation%in%plot_data$sim[i])
  plot_data$GD[i] <- out_GD[[el]][plot_data$Gen[i],plot_data$Group[i]]
  plot_data$Fis[i] <- 2*(out_Fis[[el]][plot_data$Gen[i],plot_data$Group[i]]-0.5)
  plot_data$N[i] <- out_N[[el]][plot_data$Gen[i],plot_data$Group[i]]
  plot_data$Transfers[i]<-out_Cum_transfer[[el]][plot_data$Gen[i]]
    }
plot_data$ps_nr <- ps
plot_data$Gen <- as.factor(plot_data$Gen)

if(ps==1){plot_data_all<-plot_data}else{plot_data_all<-rbind(plot_data_all,plot_data)}
} # end specific ps (combination of paramters) analysis

##########################
#generate overall dataset#
##########################

write.table(parameter_space, file = paste("PS",NG_BC,NBC_I,NI,".csv",sep="_"), sep = ",", quote = FALSE, row.names = F)                # make a CSV file in the working directory that holds the parameter space
write.table(plot_data_all, file = paste("OUT",NG_BC,NBC_I,NI,mating_system,".csv",sep="_"), sep = ",", quote = FALSE, row.names = F)   # make a CSV file in the working directory that holds the eventual output data

parameter_space <- read.csv(paste(c('PS',NG_BC,NBC_I,NI,'.csv'), collapse = '_'),header=T,sep=",")
plot_data_all <- read.csv(paste(c('OUT',NG_BC, NBC_I, NI, mating_system, '.csv'), collapse = '_'),header=T,sep=",")
plot_data_all$Gen <- as.factor(plot_data_all$Gen)

#########################
##generate plots per ps##
#########################

# specify the path to the location where you want to save a PDF containing the plots of the different combinations of parameters (ps)
destination <- paste("Boxplots",
                     NG_BC,"_",NBC_I,"_",NI,"_",mating_system,"pdf_files.pdf",sep="")

pdf(file = destination)

for (ps in unique(plot_data_all$ps_nr))
{plot_data<-plot_data_all %>% filter(ps_nr ==ps)
 GD_plot<-ggplot(plot_data[plot_data$Group==NG+1,], aes(x=Gen, y=1-GD, fill=TS)) +
 geom_boxplot(position=position_dodge(1))
 GD_plot    # this is a boxplot comparing the maintained genetic diversity for the different

 Fis_plot<-ggplot(plot_data[plot_data$Group==NG+1,], aes(x=Gen, y=Fis, fill=TS)) +
 geom_boxplot(position=position_dodge(1))
 Fis_plot   # this is a boxplot comparing the eventual inbreeding coefficient for the different

 Transfer_plot<-ggplot(plot_data, aes(x=Gen, y=Transfers, fill=TS)) +
 geom_boxplot(position=position_dodge(1))
 Transfer_plot    # this is a boxplot comparing the eventual number of transfers for the different

 Popsize_plot<-ggplot(plot_data[plot_data$Group==NG+1,], aes(x=Gen, y=N, fill=TS)) +
 geom_boxplot(position=position_dodge(1))
 Popsize_plot    # this is a boxplot comparing
 
 # generate the title to accompany the plots for each ps
 tgrob <- text_grob(paste("SCENARIO ",NG_BC,"-",NBC_I,"-",NI," PARAMETERS (K,TFG,TFHL,F_GR) = ",paste(parameter_space[ps,], collapse = '_')),color = "red", face = "bold", size = 10)
 plot_0 <- as_ggplot(tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))

p <- ggarrange(plot_0,GD_plot,ggarrange(Fis_plot, Transfer_plot, ncol = 2, labels = c("B", "C")), 
  nrow = 3,labels = c("","A"),heights=c(1,4,4))
plot(p)
}
dev.off()

## 4. Statistical analysis of the simulation outputs (linear mixed-effect model)

# finally, we will make a linear mixed-effect model to statistically analyse the simulation outputs
library(lme4)
library(lmerTest)
head(parameter_space)
head(plot_data_all)
parameter_space$ps_nr<-rownames(parameter_space)

data <- merge(parameter_space,plot_data_all,by="ps_nr")

m1 <- lmer((1-GD)~TS*TFG+(1|ps_nr), data = data[data$Gen==40,])
summary(m1) # interpret this summary as explained in the READ.ME file