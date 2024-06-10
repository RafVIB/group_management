# group_management
Code of my Biodiversity Integrated Research Project on group management at the Centre for Research and Conservation.

In my Master of Science in Biology at Ghent University, I took the course Biodiversity Integrated Research Project. I went to the Centre for Research and Conservation at Antwerp Zoo to work on a project on group management of captive breeding populations. In this repository, you can access the R code we built for simulation of captive populations managed under different scenarios. It simulates different strategies of group management to unravel what strategies are most suitable for managing certain captive populations in zoos.

# Introduction to the topic
Managing captive populations is important for the conservation of certain threatened species. European zoos cooperate closely to do so. There are two types of breeding programmes: the more intensive European Endangered Species Program or EEP and the less intensive European Studbook or ESB. Both approaches rely on collecting all data on births, deaths, transfers etc. of specific individuals in so-called studbooks. Based on this information, the coordinator of the breeding program can analyse the population and can make recommendations on what animals should move to where and what animals should mate with what animals in order to keep the population demographically and genetically sustainable. Each row in a studbook holds the known data about one specific individual. In these traditional breeding programs, it is crucial that you can recognize individual animals and that you need to collect detailed information about it. For bonobos, elephants or okapis, this works fine. But what to think of small animals that are kept in large groups, such as flamingos, sometimes kept in groups of hundreds of individuals? Zookeepers can’t tell the different individuals apart, so the only thing we can rely on to keep these population genetically diverse is exchange of ‘non-specific’ individuals. There has been little scientific attention for how to manage these populations properly, so what I did in my research project was investigating what possible management strategies work best for keeping these populations sustainable.

I investigated this via modelling in R. Here you can see a simple representation of the situation I modelled. You can see that the animals, in this case seahorses, are subdivided at three different levels. At the lowest level, they live together in small groups, so in this case let’s say one fish tank. Groups are clustered in so called ‘life support systems’. A life support system can for example represent different fish tanks coupled to a single filtering installation. Breeding takes place at the group level. In the model, animals are transferred between groups within a life support system at a certain frequency (TFG), for example once every two years. They can also be transferred between life support systems at a fixed ratio (TRHL), e.g. three times two years so six years. At the highest level, animals can also be transferred between institutes at a rate of three times three times two years, so eighteen years.
Every fish tank has a certain carrying capacity. If there are more animals in a group than possible, the surplus is transferred to the holding populations of that institute. When a group goes extinct, it then can be restocked with individuals from this holding population. That is the context of the model.
![Graphical representation of the modelled group management situation](https://github.com/RafVIB/group_management/blob/main/overview%20group%20management.png)
The main goal of the research project was to investigate the performance of different group management strategies. I included six of them in my model.
A very simple one is circular breeding. Within a life support system, the different fish tanks are arranged in a circle and animals are always transferred to the next tank in a circular movement.
A variation on this is the falconers’ breeding scheme, in which animals are first transferred to the next tank, then to the second next tank, and so on.
![Graphical representation of circular breeding and the falconers' breeding scheme](https://github.com/RafVIB/group_management/blob/main/CB%20%26%20FBS.png)
A more complicated one is called maximum avoidance of inbreeding. In this strategy, the destination is calculated with this formula. At t = 1, the formula equates 1, so animals are transferred to the next tank. At t = 2, the formula equates 2 (the second next tank), but at t = 3, the formula equates 5, so they will go to the fifth tank, and so on.
Another very simple strategy is line breeding. In this case, there are no transfers, and breeding always proceeds within the same group.
![Graphical representation of maximum avoidance of inbreeding and line breeding](https://github.com/RafVIB/group_management/blob/main/MAI%20%26%20LB.png)
The last two strategies are random breeding and mean kinship breeding. In random breeding, a random number of individuals is transferred.
Mean kinship breeding is the most complicated one. This strategy relies on calculation of the kinship between all individuals. Kinship is here defined as the probability that a pair of randomly sampled homologous alleles are identical by descent. So for us, as diploid organisms, the kinship with ourselves is 50%, because if I draw an allele and then draw an allele at the same locus again, there is a probability of 50% that it is the same allele again. The kinship between me and my sister or my parents is then, logically, 25%.
The mean kinship strategy uses kinship values to decide whether an individual has to be transferred to another population or not. More specifically, an animal is only transferred if the kinship of this animal with the animals in the destination group is, on average, lower than its kinship with the animals in its own group. This has some similarities with classic breeding programmes using studbooks, because it implies that you keep track of the descent of all individuals.
![Graphical representation of random breeding and mean kinship breeding](https://github.com/RafVIB/group_management/blob/main/RB%20%26%20MK.png)

# About the code
This is a graphic simulation of how the script is built up. First, in an input block, you specify a series of input variables: demographic variables, details about how the groups are subdivided among life support systems and institutes and so on. Next, a founder population is made as a dataframe that holds a line of information about every individual.
After that, the actual simulation begins. A series of steps is repeated for every generation: first, breeding takes place, so new individuals are added to the dataframe. Then, animals die and are removed from the dataframe. Then, if applicable, kinships between individuals are calculated. Animals are then transferred to other groups according to the selected transfer strategy and frequency, and data is stored for later analysis. In the end, the script also gives graphical representation of genetic diversity, inbreeding, number of individuals and number of transfers, crashes and rescues as an output.
![Overview of how the script is built up](https://github.com/RafVIB/group_management/blob/main/script%20overview.png)

When the script is simulating different scenarios, it renders the following graph for each simulation. On this graph, for example, the following scenario is simulated for a certain species: one institute thas has two life support systems, each with four breeding groups that have a founder size of 10 and carrying capacity of 50 animals. Animals are transferred between groups every year and the transfer ratio at the higher levels is 2.
These graphs display how the overall genetic diversity and inbreeding coefficient evolves for this simulation, and also how big the in this case eight different groups are throughout the 100 simulated generations.
![Example of the plots rendered for each simulation](https://github.com/RafVIB/group_management/blob/main/example%20graph.png)

When you let the script run multiple times for every possible combination, you can graphically inspect the simulated data for a given breeding design. Here, you can see an example of the summarized data for a very similar scenario, but in this case with two instutes, with a carrying capacity of 20 animals per groups and different transfer frequencies.
Each simulation simulates one hundred generations. Graph A summarizes how genetic diversity evolves across the generations for the different transfer strategies. Graph B does the same for inbreeding, while graph C keeps track of the number of transfers.
![Example of an overarching boxplot for a certain scenario]([https://github.com/RafVIB/group_management/assets/169991371/ee8fb61c-f61d-439e-b28e-f2c0ce0b7b19](https://github.com/RafVIB/group_management/blob/main/example%20overarching%20boxplot.png))
(What stands out here is that the mean kinship strategy comes with way less transfers than other strategies, except of course line breeding which doesn’t have any transfers at all. But this is logical: strategies such as circular breeding, Falconer’s breeding scheme and maximum avoidance of inbreeding will always transfer animals when the time to transfer has come, while mean kinship breeding will only transfer an animal if its mean kinship toe the animals in the destination group is lower than its mean kinship to the animals in its current group. This reduced number of transfers also explains why inbreeding is higher in mean kinship breeding. On the other hand, the mean kinship strategy is quite good at maintaining genetic diversity.)

For a more formal consideration of the different transfer strategies, a statistical test is also built in in the script. It fits a linear mixed-effect model on the simulation outputs. If you ask R to give a summary of the fitted model, you get an output like this, that shows the effects of the different transfer strategies and their interactions with transfer frequency. The intercept represents the circular breeding strategy. The estimates of the effects than show how much better or worse these different transfer strategies are for maintaining genetic diversity, and how significant this effect is.

# Case study: common marmoset
As an example to test the script, I here applied it to simulate some breeding strategies for the common marmoset. This is not an endangered species, but it has been studied a lot and it breeds well in captivity, so it can be used as a sort of ‘model species’ for a breeding program of an endangered but less well-studied marmoset species. I estimated the following species specific parameters from litterature and I simulated a breeding design with one institute thas has two life support systems, each with four breeding groups.
After I ran the simulations, the linear mixed-effect model rendered the following effects on genetic diversity. In this case, line breeding and mean kinship breeding are – very significantly – better at maintaining genetic diversity than circular breeding, but the effect sizes are very small. Furthermore, line breeding also comes with a very high inbreeding and mean kinship is technically a very complex strategy, so this doesn’t mean that they per se are the two best strategies to follow. For this specific case study, the interpretation is thus not that rich, but that’s alo because I simulated only one breeding design. More than one is at the moment computationally still very heavy. I’ll come back to that in the conclusion.
![Overview of the common marmoset case study](https://github.com/RafVIB/group_management/blob/main/overview%20case%20study.png)

# Conclusion
This research project has delivered a ‘protoype’ for a program that captive population managers can use to simulate how different breeding designs and transfer strategies perform in maintaining genetic diversity and avoiding inbreeding. They simply have to fill in the demographic parameters of their species, pick a series of breeding designs they want to test and interpret the output of the simulations.
This is only a first version of the program. Numerous upgrades and expansions can still be added to it. For example, the program currently works best for K species, a next version can be adapted to also simulate populations of r species in a better way. Another thing that can be upgraded is, for example, making survival age-dependent, which it isn’t now in the program, but in real life, it very often is. But the most important thing to work on further is the speed of the program. Ironically, we first worked in the program Vortex, which is a program specifically made to perform population viability analyses such as this one, but this worked too slow, so then we switched to R and started writing the whole program ourselves. Not only did this take a very long time to write, leaving little time to run the program for specific species, the program that we wrote also turned out to take a very, very long time to run.
For example, to come back to the common marmoset case, you can see in the calculation table below that the program had to run 2,700 simulations of one hundred generations. Because of course, there are already the six management strategies, the ones that I explained you with the clownfish examples, to start with.
We tested all this for three different group carrying capacities, for five transfer frequencies between groups, for three transfer ratios at the higher level and two founder group sizes. Finally, we did five different runs for every combination. In total, this means the script had to simulate 2,700 runs.
But remember that this is only for one design, the 4-2-1 breeding design, so a breeding design with one institute thas has two life support systems, each with four breeding groups.
![Calculation table for the common marmoset case study](https://github.com/RafVIB/group_management/blob/main/calculation%20table.png)
If you want to test multiple breeding designs to see which one is best, the program will take such a long time to run that it becomes practically impossible to compare these multiple breeding designs, which is a shame, because that would give us so much more usable information that just testing the different transfer strategies. So the main bottleneck that has to be fixed in the future, is running time.

# What you can find in this repository
The programme script, documented with notes throughout for readability:
- A blank programme where you have to fill in all paramaters yourself ('GMI_V2024_0023_blank');
- A programme filled in for the common marmoset example ('GMI_V2024_0023_example').

The simulation output data from the common marmoset case study:
- The CSV-file 'PS_4_2_1_' describes the parameter space (all simulated combination). Each line represent one scenario, with successively carrying capacity of a group (K_G), transfer frequency between groups (TFG), transfer ratio at the higher level (TFHL), and number of founder individuals per group (FounderSize);
- The CSV-file 'OUT_4_2_1_M_' hols the simulation output. Each line specifies a certain combination of parameters, a generation, the simulation (out of five per combination of parameters), the group, and the genetic diversity, inbreeding coefficient, group size and cumulative number of transfers in this generation.
