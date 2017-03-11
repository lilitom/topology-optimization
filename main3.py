from getmax_position import getmax_po
from findone import Formatpoprandom
from continues_format import continues_format
import time
import csv
import os
import deap
import random
from deap import base
from deap import creator
from deap import tools


creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)
toolbox = base.Toolbox()
toolbox.register("attr_bool", random.randint, 0, 1)
toolbox.register("individual", tools.initRepeat, creator.Individual,
                 toolbox.attr_bool, 400)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.08)
toolbox.register("select", tools.selTournament, tournsize=3)


def main():
    random.seed(64)
    pop = toolbox.population(n=150)
    CXPB, MUTPB, NGEN = 0.5, 0.6, 135
    #Evaluate the entire population-----------------------------------------------------------
    number=160
    itear=5
    for i in range(len(pop)):
        continues_format(pop[i], number, itear)       
    f=open('abaqusrun_input.txt','w')
    for i in range(len(pop)):
        for j in range(len(pop[i])):
            f.write(str(pop[i][j]))
        f.write('\n')
    f.close()
    os.system('abaqus cae noGUI=abaqusrun.py')
    fitnesses = []
    with open('abaqusrun_output.txt') as f:
        for i in f.readlines():
            c = i.split('\n')[0]
            fitnesses.append((float(c),))
    os.remove('abaqusrun_input.txt')
    os.remove('abaqusrun_output.txt')
    #------------------------------------------------
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
        
    # Begin the evolution
    for g in range(NGEN):
        offspring = toolbox.select(pop, len(pop))
        offspring = list(map(toolbox.clone, offspring))
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values
        for mutant in offspring:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        #-----------------------------------------------------------------------
        fitnesses=[]
        for i in range(len(invalid_ind)):
            continues_format(invalid_ind[i], number, itear)     
        f = open('abaqusrun_input.txt', 'w')
        for i in range(len(invalid_ind)):
            for j in range(len(invalid_ind[i])):
                f.write(str(invalid_ind[i][j]))
            f.write('\n')
        f.close()
        os.system('abaqus cae noGUI=abaqusrun.py')
        with open('abaqusrun_output.txt') as f:
            for i in f.readlines():
                c = i.split('\n')[0]
                fitnesses.append((float(c),))
        os.remove('abaqusrun_input.txt')
        os.remove('abaqusrun_output.txt')
        #---------------------------------------------------------------------------
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        pop[:] = offspring

        fits = [ind.fitness.values[0] for ind in pop]
        max_position = getmax_po(fits)


        length = len(pop)
        mean = sum(fits) / length
        sum2 = 0
        for i in fits:
            sum2 = sum2 + i
        std = abs(sum2 / length - mean ** 2) ** 0.5

        print("  Min %s" % min(fits))
        print("  Max %s" % max(fits))
        print("  Avg %s" % mean)
        print("  Std %s" % std)

        #-----------------------------------------------------------------------
        f = open('abaqusrun2_input.txt', 'w')
        tt=pop[max_position]
        for j in range(len(tt)):
            f.write(str(tt[j]))
        f.write('\n')
        f.close()
        os.system('abaqus cae noGUI=abaqusrun2.py')
        new_fitnesses=[]
        with open('abaqusrun2_output.txt') as f:
            for i in f.readlines():
                c = i.split('\n')[0]
                new_fitnesses.append((float(c),))
        os.remove('abaqusrun2_input.txt')
        os.remove('abaqusrun2_output.txt')
        #---------------------------------------------------------------------------
        f = open('fitness.txt', 'a')
        f.write(str(max(fits)) + '\r\n')
        f.close()
        #---------------------------------------------------------------------

    best_ind = tools.selBest(pop, 1)[0]
    print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))



if __name__ == "__main__":
    main()