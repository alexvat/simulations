from pylab import *
from core import *
import numpy
import simuPOP as sim
from simuPOP.utils import *
import sys
import time


class Model(Simulation):

    """ Class that provides facilities to run simulations in simuPOP with an island model 
  
    :param Gen: number of generations over 1 simulation
    :param loci: number of loci
    :param alleles: number of alleles
    :param dist: distance between loci (it can be either a constant number or a list)
    :param numPop: number of populations we want to simulate
    :param sizePop: population size of each subpopulation

    """


    def __init__(self, Gen,loci,alleles,dist,numPop):
        self.loci=loci
        self.alleles = alleles
        self.dist=dist
        self.numPop = numPop
        Simulation.__init__(self, Gen, loci=self.loci, alleles=self.alleles, dist=self.dist, numPop=self.numPop)
        
         
       
    
    def reset(self):
        """ important to reset the parameters of the simulation for each run"""
        Simulation.reset(self)
        
    
    
    def run(self,step=2,
                 sizePop=100,
                 infoFields=['migrate_to','fitness'],
                 recombination_rate = 0.00375,
                 migration_rate = 0.01,
                 mutation_rate = [0.00000001],
                 subPopNames = ['x','y','z','w'],
                 alleleNames = ['A','B'],
                 s1 = 0.1,
                 burnin=50,
                 **kargs):

        """ This method will perform the simulation in simuPOP.
            The above parameters are the default ones but they can easily be changed.

            :param step: the step in the generations
            :param sizePop: define the population size. If you have more than one subpopulation, this is the size of the one subpopulation.
                            if you want your subpopulations to have different size, you need to change a little the script.

            :param numPop: number of subpopulations
            :param infoFields: fields needed for simupop
            :param recombination_rate
            :param migration_rate
            :param mutation rate
            :param initFreq: initialize the frequencies for all loci
            :param subPopNames: names of the subpopulations
            :param alleleNames
            :param s1: value of the selection coefficient
            :param burnin:
        """

        self.reset()
        pop=sim.Population(size=[sizePop]*self.numPop, loci=self.loci, lociPos=list(range(self.dist, (self.dist*self.loci)+1,self.dist)), subPopNames=subPopNames, infoFields=infoFields)
    
        simu = sim.Simulator(pop)
        print("The simulation has started")
        t1 = time.time()


        mutate_snps=range(0,50)+range(51,101)

        # define the initialization of each loci based the beta distribution where a and b parameters are allele frequencies from noncoding human regions
        snps=[0.14, 0.11, 0.17, 0.11, 0.32, 0.33, 0.21, 0.11, 0.11, 0.28, 0.11, 0.12, 0.8, 0.66, 0.74, 0.68, 0.66, 0.77, 0.77, 0.76, 0.77, 0.74, 0.72, 0.11, 0.73, 0.72, 0.72, 0.72, 0.54, 0.17, 0.78, 0.64, 0.78, 0.2, 0.24, 0.25, 0.78, 0.66, 0.2, 0.14, 0.75, 0.16, 0.72, 0.18, 0.77, 0.42, 0.34, 0.7, 0.17, 0.14, 0.2, 0.46, 0.13, 0.26, 0.16, 0.13, 0.14, 0.24, 0.18, 0.36, 0.71, 0.27, 0.28, 0.25, 0.25, 0.3, 0.19, 0.14, 0.16, 0.3, 0.39, 0.16, 0.24, 0.32, 0.11, 0.18, 0.48, 0.31, 0.21, 0.15, 0.34, 0.71, 0.33, 0.18, 0.71, 0.13, 0.23, 0.2, 0.22, 0.23, 0.16, 0.23, 0.23, 0.22, 0.24, 0.82, 0.36, 0.37, 0.72, 0.16, 0.14]
        self.initFreq=[]

        
        for i in range(len(snps)):
           alpha=float(4*sizePop*migration_rate*snps[i])
           bhta=float(4*sizePop*migration_rate*(1-snps[i]))                  
           p=numpy.random.beta(alpha,bhta)
           while (p>=0.9 or p<=0.1):
                   p=numpy.random.beta(alpha,bhta)
               
           print " SNP {snp} with alpha {alpha}, bhta {bhta} and frequency {p}".format(snp=i, alpha=alpha, bhta=bhta, p=p)
           self.initFreq.append(p)

        simu.evolve(
            
            initOps=[sim.InitGenotype(freq=[self.initFreq[i], 1-self.initFreq[i]], loci=i) for i in range(len(snps))],
 

            # initialize the sex and select the 50 loci (parents)
            preOps = [sim.InitSex(maleProp=0.5,at=[0]),

                      # initialize the genotype of locus 50 at generation 0 (in the beginning of the simulation)
                      sim.PyOperator(self.genotypeBegin,at=[0]),
                      
                      # Wait 50 generations for the system to reach equilibrium
                      # Then, change the the genotype of locus 50 at generation 50 by inserting a single copy of allele 0 in one individual 
                      sim.PyOperator(self.genotypeAfter,at=[50]),

                      # function that carries out the selection proccess
                      sim.MaSelector(loci=50,wildtype=0,fitness=[1+s1, 1+s1/2, 1],begin=50, end=-1,subPops=1)],

            # recombination
            matingScheme=sim.RandomMating(ops=[
                sim.Recombinator(rates=recombination_rate)]),
                                                                                                              
            # mutation and migration of offsprings
            postOps = [

      
                sim.SNPMutator(u=mutation_rate,loci=mutate_snps),
                          
                # call function to calculate Fst and check for equilibrium state
                sim.PyOperator(self.calcFst,step=step),

                #migration
                # Here we define an island model, but this can easily be changed.
                # For more information about the migration models, please look in the documentation of SimuPOP here http://simupop.sourceforge.net/manual_svn/build/userGuide_ch7_sec3.html
                sim.Migrator(sim.utils.migrIslandRates(migration_rate,self.numPop)),
               
                # call function to save the allele frequencies
                sim.PyOperator(self.checkAlleles, step=step, param = subPopNames),
 
                                    
                # check if locus 50 is lost due to genetic drift. If yes, we terminate the simulation
                sim.Stat(alleleFreq=50,step=step,subPops=1,begin=50,end=-1),
                sim.TerminateIf('alleleFreq[50][0] == 0',step=step,begin=50,end=-1),
                
                # check the progress of the simulation
                sim.PyEval('"Gen: %d" % gen',step=step),
                sim.PyOutput('\n',step=step),
      
            ],
            gen=self.Gen
              
        )
        
 
        t2 = time.time()
        print "simulation took", t2-t1, "seconds."
       
 #       sys.stderr = old_stderr
 #      debugOutput.close()



class MultiModel(MultiSimulation):
    
    """ Class that provides facilities to run simulations in simuPOP mutiple times 

    >>> g = MultiModel(Gen=2000,Nruns=100)
    >>> g.run()
  
    :param Gen: number of generations over 1 simulation
    :param Nruns: number of runs

    """

    # the init arguments must be those of the optim_func
    def __init__(self, Gen=1000, loci=100, dist=4, alleles=2, Nruns=10, numPop=4):
        super(MultiModel, self).__init__(Gen,  loci, dist, alleles, Nruns, numPop,
            optim_func=Model)

        
    # the parameters are the same as the ones as running one simulation
    def run(self, step=10, sizePop=100, loci=100, infoFields=['migrate_to','fitness'],recombination_rate = 0.00001,
            migration_rate = 0.05, mutation_rate = [0.00000001], subPopNames = ['x','y','w','z'], alleleNames=['A','B'],
            s1 = 0.015, burnin=100, **kargs):
        super(MultiModel, self).run( step=step,
                                     sizePop=sizePop,
                                     infoFields=infoFields,
                                     recombination_rate = recombination_rate,
                                     migration_rate = migration_rate,
                                     mutation_rate = mutation_rate,
                                     subPopNames = subPopNames,
                                     alleleNames = alleleNames,
                                     s1 = s1,
                                     burnin=burnin, **kargs)
       



