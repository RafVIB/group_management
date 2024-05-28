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

I had to run this script a huge number of times, because we wanted to test a large number of scenarios. Of course, there are already the six management strategies that I explained with the clownfish examples to start with. But we wanted to test each different strategy for different breeding designs. So with breeding design I mean a combination of a certain number of institutes, each with a certain number of life support systems, each with a certain number of groups. In total, we tested 68 combinations. We tested all this for three different group carrying capacities, for five transfer frequencies between groups, 2 transfer ratios at the higher level and two founder group sizes. Finally, we did three different runs for every combination. In total, this means the script had to simulate 122,400 runs.
If you want to use the code yourself, you can of course change the code to test a specific scenario of your own choice.

In the end, you will get three summarizing graphs per scenario. Graph A summarizes how genetic diversity evolves across the generations for the different transfer strategies. Graph B does the same for inbreeding, while graph C keeps track of the number of transfers.
