""" call automatically the class to run multiple simulations.
    all the data are saved in pickle objects.
    Pickle objects, as tested here can not save data for more than 100 simulations (approximately)
    It is suggested to repeat this procedure many time.
    In the following example we run in total 4 simulations (we repeat twice the MultiModel where it runs 2 simulations each time) and we save the results.
    For each loop, a directory s is produced where it contains all the data saved from each MultiModel call """

import models as s
from pylab import *
import pickle
import numpy
import sys,os


for i in range(0,2):
    # create a new file to save the data
    file =file="s{r}".format(r=i)
    os.makedirs(file)
    os.chdir(file)

    # call the MultiModel with the parameters that we want
    g = s.MultiModel(Gen=200,loci=101,dist=4,alleles=2,Nruns=2,numPop=4)
    g.run(sizePop=200,step=2,recombination_rate=0.00375,migration_rate=0.01,s1=0.1,mutation_rate=0.00000001,subPopNames = ['x','y','z','w'],burnin=10)

    # save the data in pickle objects
    pickle.dump((g.multi.all_YSelectedLoci),open("all_YSelectedLoci.dat","w"))
    pickle.dump((g.multi.all_XSelectedLoci),open("all_XSelectedLoci.dat","w"))
    pickle.dump((g.multi.all_WSelectedLoci),open("all_WSelectedLoci.dat","w"))
    pickle.dump((g.multi.all_ZSelectedLoci),open("all_ZSelectedLoci.dat","w"))

    pickle.dump((g.multi.all_fst),open("all_fst.dat","w"))

    pickle.dump((g.multi.all_Sim_alleleFreq),open("all_Sim_alleleFreq.dat","w"))
    pickle.dump((g.multi.all_Sim_haplotypes),open("all_Sim_haplotypes.dat","w"))
    os.chdir('..')
