from pylab import *
import simuPOP as sim
from simuPOP.utils import *
import numpy
import random
import math
import time
import csv



__all__ = ["Results", "MultiResults", "Simulation", "MultiSimulation"]



class MultiResults(object):
    """A class to store results from simulation using SimuPOP.


    The results of each simulation are: allele frequencies, haplotypes, fst over time.
    So, it is convenient to store those results in well defined objects.
    This class is useful to store multiple instance of Results 
    

    """
    def __init__(self):
        self._results = []


    def _get_Gen(self):
        """ Return the number of generations"""
        return self._results[0].Gen
    Gen = property(_get_Gen)

    def _get_step(self):
        """ step to save the results"""
        return self._results[0].step
    step = property(_get_step)

    def _get_xGen(self):
        return range(self.step, self.Gen+self.step, self.step)
    xGen = property(_get_xGen)
    
    def add_results(self, res):
        self._results.append(res)

    def reset(self):
        self._results = []
    
    def _get_all_Sim_alleleFreq(self):
        """ Return all the allele frequencies from all the simulations.
            Each simulation is a dictionary which contains all the alleles frequencies from all loci.
        """
        return ([dict(x['all_allelesFreq']) for x in self._results])
    all_Sim_alleleFreq = property(_get_all_Sim_alleleFreq)

    def _get_all_fst(self):
        """Return the fst values from each simulation"""
        return [list(x['fst']) for x in self._results]
    all_fst = property(_get_all_fst)


    def _get_all_YSelectedLoci(self):
        """ return the allele frequency of the selected locus, which is locus 50 in our case, just to make it easier to produce the plots.
            The same is procedure follow the XSelectedLoci,WSelectedLoci and ZSelectedLoci.
            This information can also be extracted from the all_Sim_alleleFreq object"""
        return [list(x['YSelectedLoci']) for x in self._results]
    all_YSelectedLoci = property(_get_all_YSelectedLoci)


    def _get_all_XSelectedLoci(self):
        return [list(x['XSelectedLoci']) for x in self._results]
    all_XSelectedLoci = property(_get_all_XSelectedLoci)

    def _get_all_WSelectedLoci(self):
        return [list(x['WSelectedLoci']) for x in self._results]
    all_WSelectedLoci = property(_get_all_WSelectedLoci)

    def _get_all_ZSelectedLoci(self):
        return [list(x['ZSelectedLoci']) for x in self._results]
    all_ZSelectedLoci = property(_get_all_ZSelectedLoci)
    
    def _get_all_Sim_haplotypes(self):
        """ Return all the haplotypes from all the subpopulation from all the simulations.
        """
        return ([dict(x['all_haplotypes']) for x in self._results])
    all_Sim_haplotypes = property(_get_all_Sim_haplotypes)

    def _get_all_haplo(self):
        """ Return all the haplotypes from all the subpopulation from all the simulations.
        """
        return ([dict(x['haplo']) for x in self._results])
    all_haplo = property(_get_all_haplo) 
        
    def _get_results(self):
        return self._results
    results = property(_get_results)

    
        

class Results(dict):
    """A class that stores results from one simlation

       The results from a simulation are all similar: allele frequencies, fst values to examine
       if it has reached an equilibrium and haplotypes values for both individuals.

       """
    
    def __init__(self,Gen=2000,loci=100,alleles=2,step=10,dist=8,numPop=4):
        super(Results, self).__init__()
        self.step = step
        self._Gen = Gen
 #       self.selected_loci=0
        self._loci = loci
        self._alleles = alleles
        self.dist=8
        self._numPop=numPop

        # lists to store the allele frequencies for each allele and for each loci
        for  locus in range(self.loci):
            for allele in range(self.alleles):
                self['alleleFr{locus}{allele}' .format(locus=locus,allele=allele)]=[]
        

       
        self['fst'] = []

        # dictionary that will save all the allele frequencies from all loci for all the alleles over the generations
        self['all_allelesFreq']={}
        self['all_haplotypes']={}
        self['haplo']={}
     
        # define the list of names of the subpopulations 
        sub=['x','y','z','w']
        for  i in sub:
            # initialization of the lists where we wiil save the haplotypes of each subpopulations when the allele frequency of the selected locus is close to 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 and 1 

            self['all_haplotypes']['{pop}01'.format(pop=i)]=[]
            self['all_haplotypes']['{pop}02'.format(pop=i)]=[]
            self['all_haplotypes']['{pop}03'.format(pop=i)]=[]
            self['all_haplotypes']['{pop}04'.format(pop=i)]=[]
            self['all_haplotypes']['{pop}05'.format(pop=i)]=[]
            self['all_haplotypes']['{pop}06'.format(pop=i)]=[]
            self['all_haplotypes']['{pop}07'.format(pop=i)]=[]
            self['all_haplotypes']['{pop}08'.format(pop=i)]=[]
            self['all_haplotypes']['{pop}09'.format(pop=i)]=[]
            self['all_haplotypes']['{pop}1'.format(pop=i)]=[]
            self['haplo']['{pop}'.format(pop=i)]=[]
       
        self['YSelectedLoci']=[]
        self['XSelectedLoci']=[]
        self['WSelectedLoci']=[]
        self['ZSelectedLoci']=[]
  
    
    def _getGen(self):    
        return self._Gen
    Gen = property(_getGen)

    def _getLoci(self):    
        return self._loci
    loci = property(_getLoci)

    def _getAlleles(self):    
        return self._alleles
    alleles = property(_getAlleles)
    
    def _get_xGen(self):
        return range(self.step, self.Gen+self.step, self.step)
    xGen = property(_get_xGen)

    def _get_fst(self):
        return self['fst']
    fst = property(_get_fst)

    def _get_YSelectedLoci(self):
        return self['YSelectedLoci']
    YSelectedLoci = property(_get_YSelectedLoci)

    def _get_XSelectedLoci(self):
        return self['XSelectedLoci']
    XSelectedLoci = property(_get_XSelectedLoci)

    def _get_WSelectedLoci(self):
        return self['WSelectedLoci']
    WSelectedLoci = property(_get_WSelectedLoci)

    def _get_ZSelectedLoci(self):
        return self['ZSelectedLoci']
    ZSelectedLoci = property(_get_ZSelectedLoci)
    
    def _get_all_allelesFreq(self):
        return self['all_allelesFreq']
    all_allelesFreq = property(_get_all_allelesFreq)

    def _get_all_haplotypes(self):
        return self['all_haplotypes']
    all_haplotypes = property(_get_all_haplotypes)

    def _get_haplo(self):
        return self['haplo']
    haplo = property(_get_haplo)

       

class Simulation(object):
    """ This class provides facilities to run a simulation with simupop."""

    def __init__(self, Gen=None, loci = None, alleles = None, dist= None,numPop=None):

        self.dist = dist
        self.Gen = Gen
        self.loci=loci
        self.alleles=alleles
        self.numPop=numPop
        

        self.results = Results(Gen=self.Gen, loci = self.loci, alleles = self.alleles, step=10, dist = self.dist, numPop = self.numPop)
        
        
    def _setGen(self, value):
        if value<0:
            raise ValueError("N must be positive.")
        else:
            self._Gen = value

    def _getGen(self):   
        return self._Gen
    Gen = property(_getGen, _setGen)


    def _setLoci(self, value):
        if value<0:
            raise ValueError("Loci Number must be positive.")
        else:
            self._loci = value

    def _getLoci(self):   
        return self._loci
    loci = property(_getLoci, _setLoci)

    def _setAlleles(self, value):
        if value<0:
            raise ValueError("Number of alleles must be positive.")
        else:
            self._alleles = value

    def _getAlleles(self):   
        return self._alleles
    alleles = property(_getAlleles, _setAlleles)


    def reset(self):
        self.results = Results(step=10, Gen=self.Gen, loci = self.loci, alleles = self.alleles)

    def genotypeBegin(self,pop):
        """ initialize the genotype of locus 50 (selected locus) to be monomorphic which means that all individuals will have the the allele 1 in the beginning of the simulation. """

        for ind in pop.individuals():
             geno = ind.genotype(0)
             geno[50]=1
             geno = ind.genotype(1)
             geno[50]=1
  
        return True
        
     
    def genotypeAfter(self,pop):
        """ model the hard sweep by inserting manually a single copy (allele 0) in the individual 2 in subpopulation 1(Y) at locus 50 """
 
        ind=pop.individual(2,1)
        geno = ind.genotype(0)
        geno[50]=0
        return True
 
     
    def checkAlleles(self,pop,param):
        """ save all allele frequencies of all loci over all generations in self.results.alleleFr.
            All data are saved in a dictionary self.results.
            To acquire the data from a specific loci and allele, here is an example:
            self.results['alleleFreq01'] where 0 corresponds to the loci and 1 to the allele.

            It also saves the allele frequencies from the last generation from all loci
            in the following format: [allele0 of loci 0, allele 1 of loci 0, ... allele N of loci N]

        """
        # store the allele frequencies from subpopulation 1 (Y) from all loci 
        
        sim.stat(pop,alleleFreq=range(self.loci),subPops=[1])                                                
        for loci in range(self.loci):
            for allele in range(self.alleles):
                self.results['alleleFr{loci}{allele}'.format(loci=loci,allele=allele)].append(pop.dvars().alleleFreq[loci][allele])
                    
                # add all of them in the dictionary that I can extract from the object
                lociAllele = 'alleleFr{loci}{allele}'.format(loci=loci,allele=allele)
                self.results['all_allelesFreq'][str(lociAllele)]=list((self.results['alleleFr{loci}{allele}'.format(loci=loci,allele=allele)]))

        # save the allele frequency of allele 0 (which is selected) from locus 50 from all SubPopulations
        for  i in range(0,len(param)):
            sim.stat(pop,alleleFreq=range(self.loci),subPops=param[i])
            if i ==0:
                self.results['XSelectedLoci'].append(pop.dvars().alleleFreq[50][0])
            if i ==1:
                self.results['YSelectedLoci'].append(pop.dvars().alleleFreq[50][0])
            if i ==2:
                self.results['ZSelectedLoci'].append(pop.dvars().alleleFreq[50][0])
            if i ==3:
                self.results['WSelectedLoci'].append(pop.dvars().alleleFreq[50][0])


        # save the haplotypes for each suppopulation 
        sim.stat(pop,alleleFreq=range(self.loci),subPops=[1])
       
        for  i in param:
            if i=='x':
                popInd=0
            if i=='y':
                popInd=1
            if i=='z':
                popInd=2
            if i=='w':
                popInd=3

            # save a sample of the haplotypes from 100 individuals from all populations when the allele frequency of locus 50 is close to 0.1         
            if float(pop.dvars().alleleFreq[50][0])>0.08 and (pop.dvars().alleleFreq[50][0])<0.11: 
                if self.results['all_haplotypes']['{pop}01'.format(pop=i)]==[]:
                     for ind in random.sample(range(0,pop.subPopSize(i)), 100):
                        self.results['all_haplotypes']['{pop}01'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(0)))
                        self.results['all_haplotypes']['{pop}01'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(1)))

            # save a sample of the haplotypes from 100 individuals from all pops when the allele frequency of locus 50 is close to 0.2                        
            if (pop.dvars().alleleFreq[50][0])>0.18 and (pop.dvars().alleleFreq[50][0])<0.22:
                if self.results['all_haplotypes']['{pop}02'.format(pop=i)]==[]:
                    for ind in random.sample(range(0,pop.subPopSize(i)), 100):
                        self.results['all_haplotypes']['{pop}02'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(0)))
                        self.results['all_haplotypes']['{pop}02'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(1)))

            # save a sample of the haplotypes from 100 individuals from all pops when the allele frequency of locus 50 is close to 0.3             
            if (pop.dvars().alleleFreq[50][0])>0.28 and (pop.dvars().alleleFreq[50][0])<0.32:
                if self.results['all_haplotypes']['{pop}03'.format(pop=i)]==[]:
                    for ind in random.sample(range(0,pop.subPopSize(i)), 100):
                        self.results['all_haplotypes']['{pop}03'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(0)))
                        self.results['all_haplotypes']['{pop}03'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(1)))

            # save a sample of the haplotypes from 100 individuals from all pops when the allele frequency of locus 50 is close to 0.4  
            if (pop.dvars().alleleFreq[50][0])>0.38 and (pop.dvars().alleleFreq[50][0])<0.42:
                if self.results['all_haplotypes']['{pop}04'.format(pop=i)]==[]:
                    for ind in random.sample(range(0,pop.subPopSize(i)), 100):
                        self.results['all_haplotypes']['{pop}04'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(0)))
                        self.results['all_haplotypes']['{pop}04'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(1)))

            # save a sample of the haplotypes from 100 individuals from all pops when the allele frequency of locus 50 is close to 0.5  
            if (pop.dvars().alleleFreq[50][0])>0.48 and (pop.dvars().alleleFreq[50][0])<0.52:
                if self.results['all_haplotypes']['{pop}05'.format(pop=i)]==[]:
                    for ind in random.sample(range(0,pop.subPopSize(i)), 100):
                        self.results['all_haplotypes']['{pop}05'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(0)))
                        self.results['all_haplotypes']['{pop}05'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(1)))

            # save a sample of the haplotypes from 100 individuals from all pops when the allele frequency of locus 50 is close to 0.6 
            if (pop.dvars().alleleFreq[50][0])>0.58 and (pop.dvars().alleleFreq[50][0])<0.62:
                if self.results['all_haplotypes']['{pop}06'.format(pop=i)]==[]:
                        
                    for ind in random.sample(range(0,pop.subPopSize(i)), 100):                      
                        self.results['all_haplotypes']['{pop}06'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(0)))
                        self.results['all_haplotypes']['{pop}06'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(1)))
                        
            # save a sample of the haplotypes from 100 individuals from all pops when the allele frequency of locus 50 is close to 0.7 
            if (pop.dvars().alleleFreq[50][0])>0.68 and (pop.dvars().alleleFreq[50][0])<0.72:
                 if self.results['all_haplotypes']['{pop}07'.format(pop=i)]==[]:
                      for ind in random.sample(range(0,pop.subPopSize(i)), 100):
                         self.results['all_haplotypes']['{pop}07'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(0)))
                         self.results['all_haplotypes']['{pop}07'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(1)))
                         
            # save a sample of the haplotypes from 100 individuals from all pops when the allele frequency of locus 50 is close to 0.8 
            if (pop.dvars().alleleFreq[50][0])>0.78 and (pop.dvars().alleleFreq[50][0])<0.82:
                 if self.results['all_haplotypes']['{pop}08'.format(pop=i)]==[]:
                    for ind in random.sample(range(0,pop.subPopSize(i)), 100):
                        self.results['all_haplotypes']['{pop}08'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(0)))
                        self.results['all_haplotypes']['{pop}08'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(1)))
                        
            # save a sample of the haplotypes from 100 individuals from all pops when the allele frequency of locus 50 is close to 0.9 
            if (pop.dvars().alleleFreq[50][0])>0.88 and (pop.dvars().alleleFreq[50][0])<0.92:
                if self.results['all_haplotypes']['{pop}09'.format(pop=i)]==[]:
                    for ind in random.sample(range(0,pop.subPopSize(i)), 100):
                        self.results['all_haplotypes']['{pop}09'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(0)))
                        self.results['all_haplotypes']['{pop}09'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(1)))

            # save a sample of the haplotypes from 100 individuals from all pops when the allele frequency of locus 50 is close to 1 
            if (pop.dvars().alleleFreq[50][0])>0.97 and (pop.dvars().alleleFreq[50][0])<0.99:
                if self.results['all_haplotypes']['{pop}1'.format(pop=i)]==[]:
                    for ind in random.sample(range(0,pop.subPopSize(i)), 100):
                        self.results['all_haplotypes']['{pop}1'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(0)))
                        self.results['all_haplotypes']['{pop}1'.format(pop=i)].append(list(pop.individual(ind,popInd).genotype(1)))
         
        return True


    
    def calcFst(self,pop):
        """ Calculate the Fst values for 1 simulation based on all loci

        """
 
        sim.stat(pop, structure=range(self.loci),vars=['F_st'])
        self.results['fst'].append(pop.dvars().F_st)
        
                
        return True

    
        

    def run(self, *args, **kargs):
        """This method should perform the simulation and will save all the results in the object results. """
        raise NotImplementedError
                

class MultiSimulation(object):
    """
    :param optim_func: a class performing the simulations that we want multiple times.
    called run(). func has a results dictionary to be stored Resulst.

    """

    def __init__(self, Gen=1000,loci=100, dist=8, alleles=2, Nruns=10, numPop=4, optim_func=None,
        optim_params={}):
        self.Gen = Gen
        self.loci = loci
        self.dist = dist
        self.alleles = alleles
        self.numPop = numPop
        self.Nruns = Nruns
        if optim_func == None:
            raise ValueError("optim_func should be provided.")
        try:
            self.func = optim_func(self.Gen, self.loci, self.alleles, self.dist, self.numPop, **optim_params)
        except:
            try:
                self.func = optim_func(**optim_params)
            except Exception:
                raise Exception

        self.multi = MultiResults()

    def reset(self):
        self.multi.reset()

    def run(self, *args, **kargs):
        import time
        for i in range(0, self.Nruns):

            filename = ('Sim %s',i)
            t1 = time.time()
            print "Running simulation %s over %s " % (i, self.Nruns)
            self.func.reset()
            self.func.run(*args, **kargs)
            
            # add the results of each simulation
            self.multi.add_results(self.func.results)

            t2 = time.time()
            print t2-t1 
